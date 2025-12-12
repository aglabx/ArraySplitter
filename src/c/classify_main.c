/*
 * classify_main.c - CLI, I/O, and main classification pipeline
 *
 * Usage: arraysplitter_classify_c -i monomers.tsv -o output_prefix [options]
 */

#include "classify.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

/* Initialize classification context */
ClassifyContext* classify_init(void) {
    ClassifyContext *ctx = calloc(1, sizeof(ClassifyContext));
    if (!ctx) return NULL;

    ctx->min_shared_kmers = DEFAULT_MIN_SHARED_KMERS;
    ctx->jaccard_threshold = DEFAULT_JACCARD_THRESHOLD;
    ctx->skip_flanks = true;
    ctx->num_threads = 4;
    ctx->verbose = false;

    return ctx;
}

/* Free classification context */
void classify_free(ClassifyContext *ctx) {
    if (!ctx) return;

    /* Free monomers */
    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        free(ctx->monomers[i].sequence_id);
        free(ctx->monomers[i].orientation);
        free(ctx->monomers[i].sequence);
    }
    free(ctx->monomers);

    /* Free k-mer index */
    for (int i = 0; i < KMER_BUCKETS; i++) {
        free(ctx->buckets[i].hits);
    }

    /* Free components */
    for (uint32_t i = 0; i < ctx->component_count; i++) {
        free(ctx->components[i].member_ids);
    }
    free(ctx->components);

    /* Free subcomponents */
    for (uint32_t i = 0; i < ctx->subcomponent_count; i++) {
        free(ctx->subcomponents[i].member_ids);
    }
    free(ctx->subcomponents);

    free(ctx);
}

/* Parse monomer type from string */
static MonomerType parse_type(const char *type_str) {
    if (strcmp(type_str, "LEFT_FLANK") == 0) return TYPE_LEFT_FLANK;
    if (strcmp(type_str, "RIGHT_FLANK") == 0) return TYPE_RIGHT_FLANK;
    return TYPE_MONOMER;
}

/* Load monomers from TSV file */
int load_monomers_tsv(ClassifyContext *ctx, const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return -1;
    }

    char line[1048576];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -1;
    }

    ctx->monomer_capacity = 10000;
    ctx->monomers = malloc(ctx->monomer_capacity * sizeof(Monomer));
    ctx->monomer_count = 0;

    uint32_t line_num = 1;
    while (fgets(line, sizeof(line), fp)) {
        line_num++;

        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') line[--len] = '\0';
        if (len > 0 && line[len-1] == '\r') line[--len] = '\0';
        if (len == 0) continue;

        char *seq_id = strtok(line, "\t");
        char *orientation = strtok(NULL, "\t");
        char *index_str = strtok(NULL, "\t");
        char *type_str = strtok(NULL, "\t");
        char *length_str = strtok(NULL, "\t");
        char *is_flank_str = strtok(NULL, "\t");
        char *sequence = strtok(NULL, "\t");

        if (!seq_id || !orientation || !index_str || !type_str ||
            !length_str || !is_flank_str || !sequence) {
            fprintf(stderr, "Warning: Malformed line %u, skipping\n", line_num);
            continue;
        }

        if (ctx->monomer_count >= ctx->monomer_capacity) {
            ctx->monomer_capacity *= 2;
            ctx->monomers = realloc(ctx->monomers,
                                    ctx->monomer_capacity * sizeof(Monomer));
        }

        Monomer *mono = &ctx->monomers[ctx->monomer_count];
        memset(mono, 0, sizeof(Monomer));

        mono->sequence_id = strdup(seq_id);
        mono->orientation = strdup(orientation);
        mono->index = (uint32_t)atoi(index_str);
        mono->type = parse_type(type_str);
        mono->length = (uint32_t)atoi(length_str);
        mono->is_flank = (strcmp(is_flank_str, "TRUE") == 0 ||
                          strcmp(is_flank_str, "true") == 0 ||
                          strcmp(is_flank_str, "1") == 0);
        mono->sequence = strdup(sequence);

        mono->component_id = UINT32_MAX;
        mono->subcomponent_id = UINT32_MAX;
        mono->jaccard_to_rep = 0.0f;

        ctx->monomer_count++;
    }

    fclose(fp);

    if (ctx->verbose) {
        fprintf(stderr, "Loaded %u monomers from %s\n", ctx->monomer_count, filename);
    }

    return 0;
}

/* Main classification pipeline */
int classify_monomers(ClassifyContext *ctx) {
    clock_t start = clock();

    if (ctx->verbose) {
        fprintf(stderr, "Starting classification of %u monomers\n", ctx->monomer_count);
        fprintf(stderr, "Parameters: min_shared_kmers=%d, jaccard_threshold=%.2f\n",
                ctx->min_shared_kmers, ctx->jaccard_threshold);
    }

    /* Step 1: Build k-mer bitvectors and index (combined for efficiency) */
    if (ctx->verbose) {
        fprintf(stderr, "Step 1: Building k-mer bitvectors and index...\n");
    }
    build_kmer_index(ctx);

    /* Step 3: Clustering */
    if (ctx->verbose) {
        fprintf(stderr, "Step 3: Clustering...\n");
    }
    cluster_monomers(ctx);

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    if (ctx->verbose) {
        fprintf(stderr, "Classification complete in %.2f seconds\n", elapsed);
        fprintf(stderr, "  Components: %u\n", ctx->component_count);
        fprintf(stderr, "  Subcomponents: %u\n", ctx->subcomponent_count);
    }

    return 0;
}

/* Write components TSV output */
int write_components_tsv(ClassifyContext *ctx, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write to %s\n", filename);
        return -1;
    }

    fprintf(fp, "sequence_id\torientation\tindex\ttype\tlength\tcomponent_id\tsubcomponent_id\tjaccard_to_rep\tn_distinct_kmers\n");

    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        Monomer *mono = &ctx->monomers[i];

        const char *type_str;
        switch (mono->type) {
            case TYPE_LEFT_FLANK: type_str = "LEFT_FLANK"; break;
            case TYPE_RIGHT_FLANK: type_str = "RIGHT_FLANK"; break;
            default: type_str = "MONOMER"; break;
        }

        int comp_id = (mono->component_id == UINT32_MAX) ? -1 : (int)mono->component_id;
        int subcomp_id = (mono->subcomponent_id == UINT32_MAX) ? -1 : (int)mono->subcomponent_id;

        fprintf(fp, "%s\t%s\t%u\t%s\t%u\t%d\t%d\t%.4f\t%u\n",
                mono->sequence_id,
                mono->orientation,
                mono->index,
                type_str,
                mono->length,
                comp_id,
                subcomp_id,
                mono->jaccard_to_rep,
                mono->bitvec.popcount);
    }

    fclose(fp);

    if (ctx->verbose) {
        fprintf(stderr, "Wrote %s\n", filename);
    }

    return 0;
}

/* Write component statistics */
int write_component_stats(ClassifyContext *ctx, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot write to %s\n", filename);
        return -1;
    }

    fprintf(fp, "component_id\tn_monomers\tn_subcomponents\tavg_distinct_kmers\n");

    for (uint32_t c = 0; c < ctx->component_count; c++) {
        Component *comp = &ctx->components[c];

        /* Count subcomponents */
        uint32_t n_subs = 0;
        for (uint32_t s = 0; s < ctx->subcomponent_count; s++) {
            if (ctx->subcomponents[s].component_id == c) {
                n_subs++;
            }
        }

        /* Average distinct k-mers */
        uint64_t total_kmers = 0;
        for (uint32_t i = 0; i < comp->member_count; i++) {
            total_kmers += ctx->monomers[comp->member_ids[i]].bitvec.popcount;
        }
        float avg_kmers = (float)total_kmers / comp->member_count;

        fprintf(fp, "%u\t%u\t%u\t%.1f\n",
                c, comp->member_count, n_subs, avg_kmers);
    }

    fclose(fp);

    if (ctx->verbose) {
        fprintf(stderr, "Wrote %s\n", filename);
    }

    return 0;
}

/* Print usage */
static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s -i monomers.tsv -o output_prefix [options]\n\n", prog);
    fprintf(stderr, "Classify monomers into components using shared 8-mer content.\n\n");
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -i, --input FILE      Input monomers TSV file\n");
    fprintf(stderr, "  -o, --output PREFIX   Output file prefix\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -m, --min-shared N    Min shared k-mers for connectivity (default: %d)\n",
            DEFAULT_MIN_SHARED_KMERS);
    fprintf(stderr, "  -j, --jaccard F       Jaccard threshold for subcomponents (default: %.2f)\n",
            DEFAULT_JACCARD_THRESHOLD);
    fprintf(stderr, "  -t, --threads N       Number of threads (default: 4)\n");
    fprintf(stderr, "  --include-flanks      Include LEFT_FLANK and RIGHT_FLANK\n");
    fprintf(stderr, "  -v, --verbose         Verbose output\n");
    fprintf(stderr, "  -h, --help            Show this help\n\n");
    fprintf(stderr, "Output files:\n");
    fprintf(stderr, "  PREFIX.components.tsv   Monomer classification results\n");
    fprintf(stderr, "  PREFIX.component_stats.tsv  Component statistics\n");
}

int main(int argc, char *argv[]) {
    char *input_file = NULL;
    char *output_prefix = NULL;

    ClassifyContext *ctx = classify_init();
    if (!ctx) {
        fprintf(stderr, "Error: Failed to initialize context\n");
        return 1;
    }

    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"min-shared", required_argument, 0, 'm'},
        {"jaccard", required_argument, 0, 'j'},
        {"threads", required_argument, 0, 't'},
        {"include-flanks", no_argument, 0, 'F'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:m:j:t:vh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'm':
                ctx->min_shared_kmers = atoi(optarg);
                break;
            case 'j':
                ctx->jaccard_threshold = (float)atof(optarg);
                break;
            case 't':
                ctx->num_threads = atoi(optarg);
                break;
            case 'F':
                ctx->skip_flanks = false;
                break;
            case 'v':
                ctx->verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                classify_free(ctx);
                return 0;
            default:
                print_usage(argv[0]);
                classify_free(ctx);
                return 1;
        }
    }

    if (!input_file || !output_prefix) {
        fprintf(stderr, "Error: Input file and output prefix are required\n\n");
        print_usage(argv[0]);
        classify_free(ctx);
        return 1;
    }

    if (load_monomers_tsv(ctx, input_file) != 0) {
        classify_free(ctx);
        return 1;
    }

    if (classify_monomers(ctx) != 0) {
        classify_free(ctx);
        return 1;
    }

    char filename[1024];

    snprintf(filename, sizeof(filename), "%s.components.tsv", output_prefix);
    write_components_tsv(ctx, filename);

    snprintf(filename, sizeof(filename), "%s.component_stats.tsv", output_prefix);
    write_component_stats(ctx, filename);

    printf("Classification complete:\n");
    printf("  Input: %u monomers\n", ctx->monomer_count);
    printf("  Components: %u\n", ctx->component_count);
    printf("  Subcomponents: %u\n", ctx->subcomponent_count);
    printf("  Output: %s.components.tsv\n", output_prefix);

    classify_free(ctx);
    return 0;
}
