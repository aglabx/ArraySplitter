/*
 * minhash_main.c - Main entry point for MinHash-based classification
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "minhash.h"

/* Initialize context */
MHContext* mh_init(void) {
    MHContext *ctx = calloc(1, sizeof(MHContext));

    ctx->monomer_capacity = 1024;
    ctx->monomers = malloc(ctx->monomer_capacity * sizeof(MHMonomer));

    ctx->level1_threshold = LEVEL1_THRESHOLD;
    ctx->level2_threshold = LEVEL2_THRESHOLD;
    ctx->level3_threshold = LEVEL3_THRESHOLD;
    ctx->num_threads = 1;
    ctx->verbose = true;

    pthread_mutex_init(&ctx->uf_mutex, NULL);
    lsh_init(&ctx->lsh);

    return ctx;
}

/* Free context */
void mh_free(MHContext *ctx) {
    if (!ctx) return;

    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        free(ctx->monomers[i].sequence);
        free(ctx->monomers[i].sequence_id);
        free(ctx->monomers[i].orientation);
    }
    free(ctx->monomers);

    pthread_mutex_destroy(&ctx->uf_mutex);
    lsh_free(&ctx->lsh);

    if (ctx->subfamilies) {
        for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
            free(ctx->subfamilies[i].members);
        }
        free(ctx->subfamilies);
    }

    if (ctx->families) {
        for (uint32_t i = 0; i < ctx->family_count; i++) {
            free(ctx->families[i].members);
        }
        free(ctx->families);
    }

    if (ctx->superfamilies) {
        for (uint32_t i = 0; i < ctx->superfamily_count; i++) {
            free(ctx->superfamilies[i].members);
        }
        free(ctx->superfamilies);
    }

    free(ctx);
}

/* Parse monomer type */
static int parse_type(const char *s) {
    if (strcmp(s, "LEFT_FLANK") == 0) return 0;
    if (strcmp(s, "MONOMER") == 0) return 1;
    if (strcmp(s, "RIGHT_FLANK") == 0) return 2;
    return 1;
}

/* Load monomers from TSV file */
int mh_load_monomers(MHContext *ctx, const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open %s\n", filename);
        return -1;
    }

    char *line = NULL;
    size_t len = 0;
    int line_num = 0;

    while (getline(&line, &len, f) != -1) {
        line_num++;

        /* Skip header */
        if (line_num == 1) continue;

        /* Remove newline */
        char *nl = strchr(line, '\n');
        if (nl) *nl = '\0';

        /* Parse TSV: sequence_id, orientation, index, type, length, is_flank, sequence */
        char *seq_id = strtok(line, "\t");
        char *orient = strtok(NULL, "\t");
        char *idx_s = strtok(NULL, "\t");
        char *type_s = strtok(NULL, "\t");
        char *len_s = strtok(NULL, "\t");
        char *flank_s = strtok(NULL, "\t");
        char *seq = strtok(NULL, "\t");

        if (!seq_id || !orient || !idx_s || !type_s || !len_s || !flank_s || !seq) {
            continue;  /* Skip malformed lines */
        }

        /* Grow array if needed */
        if (ctx->monomer_count >= ctx->monomer_capacity) {
            ctx->monomer_capacity *= 2;
            ctx->monomers = realloc(ctx->monomers,
                                    ctx->monomer_capacity * sizeof(MHMonomer));
        }

        MHMonomer *m = &ctx->monomers[ctx->monomer_count];
        memset(m, 0, sizeof(MHMonomer));

        m->id = ctx->monomer_count;
        m->sequence_id = strdup(seq_id);
        m->orientation = strdup(orient);
        m->index = atoi(idx_s);
        m->type = parse_type(type_s);
        m->length = strlen(seq);
        m->is_flank = (strcmp(flank_s, "TRUE") == 0);
        m->sequence = strdup(seq);

        ctx->monomer_count++;
    }

    free(line);
    fclose(f);

    if (ctx->verbose) {
        fprintf(stderr, "Loaded %u monomers from %s\n", ctx->monomer_count, filename);
    }

    return 0;
}

/* Main classification pipeline */
int mh_classify(MHContext *ctx) {
    clock_t start = clock();

    /* Step 1: Compute MinHash signatures */
    compute_all_signatures(ctx);

    /* Step 2: Level 1 - Subfamilies */
    find_subfamilies(ctx);

    /* Step 3: Level 2 - Families */
    find_families(ctx);

    /* Step 4: Level 3 - Superfamilies */
    find_superfamilies(ctx);

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    if (ctx->verbose) {
        fprintf(stderr, "\n=== Summary ===\n");
        fprintf(stderr, "  Input monomers: %u\n", ctx->monomer_count);
        fprintf(stderr, "  Subfamilies: %u\n", ctx->subfamily_count);
        fprintf(stderr, "  Families: %u\n", ctx->family_count);
        fprintf(stderr, "  Superfamilies: %u\n", ctx->superfamily_count);
        fprintf(stderr, "  Time: %.2f seconds\n", elapsed);
    }

    return 0;
}

/* Write results */
int mh_write_results(MHContext *ctx, const char *prefix) {
    char filename[512];

    /* Main output: monomer assignments */
    snprintf(filename, sizeof(filename), "%s.hierarchy.tsv", prefix);
    FILE *f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Error: cannot create %s\n", filename);
        return -1;
    }

    fprintf(f, "monomer_id\tsequence_id\torientation\tindex\ttype\tlength\t"
               "subfamily_id\tfamily_id\tsuperfamily_id\tsequence\n");

    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        MHMonomer *m = &ctx->monomers[i];
        const char *type_str = m->type == 0 ? "LEFT_FLANK" :
                               m->type == 1 ? "MONOMER" : "RIGHT_FLANK";
        fprintf(f, "%u\t%s\t%s\t%u\t%s\t%u\t%u\t%u\t%u\t%s\n",
                m->id, m->sequence_id, m->orientation, m->index,
                type_str, m->length,
                m->subfamily_id, m->family_id, m->superfamily_id,
                m->sequence);
    }
    fclose(f);

    /* Stats output */
    snprintf(filename, sizeof(filename), "%s.stats.tsv", prefix);
    f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Error: cannot create %s\n", filename);
        return -1;
    }

    fprintf(f, "level\tcluster_id\tn_members\trepresentative_id\n");

    for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
        fprintf(f, "1\t%u\t%u\t%u\n",
                ctx->subfamilies[i].id,
                ctx->subfamilies[i].member_count,
                ctx->subfamilies[i].representative_id);
    }

    for (uint32_t i = 0; i < ctx->family_count; i++) {
        fprintf(f, "2\t%u\t%u\t%u\n",
                ctx->families[i].id,
                ctx->families[i].member_count,
                ctx->families[i].representative_id);
    }

    for (uint32_t i = 0; i < ctx->superfamily_count; i++) {
        fprintf(f, "3\t%u\t%u\t%u\n",
                ctx->superfamilies[i].id,
                ctx->superfamilies[i].member_count,
                ctx->superfamilies[i].representative_id);
    }

    fclose(f);

    if (ctx->verbose) {
        fprintf(stderr, "\nOutput written to:\n");
        fprintf(stderr, "  %s.hierarchy.tsv\n", prefix);
        fprintf(stderr, "  %s.stats.tsv\n", prefix);
    }

    return 0;
}

/* Print usage */
static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s -i input.tsv -o output_prefix [options]\n\n", prog);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -i FILE    Input TSV file with monomers\n");
    fprintf(stderr, "  -o PREFIX  Output prefix\n");
    fprintf(stderr, "  -t N       Number of threads (default: 1)\n");
    fprintf(stderr, "  -1 FLOAT   Level 1 threshold (default: %.2f)\n", LEVEL1_THRESHOLD);
    fprintf(stderr, "  -2 FLOAT   Level 2 threshold (default: %.2f)\n", LEVEL2_THRESHOLD);
    fprintf(stderr, "  -3 FLOAT   Level 3 threshold (default: %.2f)\n", LEVEL3_THRESHOLD);
    fprintf(stderr, "  -q         Quiet mode\n");
    fprintf(stderr, "  -h         Show this help\n");
}

int main(int argc, char *argv[]) {
    char *input = NULL;
    char *output = NULL;
    float t1 = LEVEL1_THRESHOLD;
    float t2 = LEVEL2_THRESHOLD;
    float t3 = LEVEL3_THRESHOLD;
    int num_threads = 1;
    bool quiet = false;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:t:1:2:3:qh")) != -1) {
        switch (opt) {
            case 'i': input = optarg; break;
            case 'o': output = optarg; break;
            case 't': num_threads = atoi(optarg); break;
            case '1': t1 = atof(optarg); break;
            case '2': t2 = atof(optarg); break;
            case '3': t3 = atof(optarg); break;
            case 'q': quiet = true; break;
            case 'h':
            default:
                print_usage(argv[0]);
                return opt == 'h' ? 0 : 1;
        }
    }

    if (!input || !output) {
        print_usage(argv[0]);
        return 1;
    }

    MHContext *ctx = mh_init();
    ctx->level1_threshold = t1;
    ctx->level2_threshold = t2;
    ctx->level3_threshold = t3;
    ctx->num_threads = num_threads > 0 ? num_threads : 1;
    ctx->verbose = !quiet;

    if (mh_load_monomers(ctx, input) != 0) {
        mh_free(ctx);
        return 1;
    }

    if (mh_classify(ctx) != 0) {
        mh_free(ctx);
        return 1;
    }

    if (mh_write_results(ctx, output) != 0) {
        mh_free(ctx);
        return 1;
    }

    mh_free(ctx);
    return 0;
}
