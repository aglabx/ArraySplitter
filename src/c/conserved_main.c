/*
 * Conserved Region Finder (Multi-threaded)
 *
 * Finds highly conserved regions in edit distance output.
 * A conserved region is a stretch of consecutive monomer pairs
 * with similarity above threshold.
 *
 * Usage: arraysplitter_conserved -i file.ed -o output [-t threshold] [-m min_length]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <pthread.h>

#define MAX_LINE_LENGTH 4096
#define MAX_ID_LENGTH 2048

/* Edit distance pair from .ed file */
typedef struct {
    int idx1;
    int idx2;
    size_t len1;
    size_t len2;
    size_t edit_dist;
    double similarity;
    char type1[32];
    char type2[32];
} EdPair;

/* Array of pairs for one sequence */
typedef struct {
    char *sequence_id;
    char *orientation;
    size_t n_monomers;
    EdPair *pairs;
    size_t pair_count;
    size_t pair_capacity;

    /* Results (filled by worker) */
    char *result_regions;
    size_t result_regions_len;
    char *result_stats;
    size_t result_stats_len;
} EdArray;

/* Work queue for threads */
typedef struct {
    EdArray *arrays;
    size_t count;
    size_t next_item;
    pthread_mutex_t mutex;
    size_t completed;
    double threshold;
    int min_region_length;
} WorkQueue;

/* Conserved region */
typedef struct {
    int start_idx;
    int end_idx;
    int length;
    double avg_similarity;
    double min_similarity;
    size_t total_monomer_length;  /* Sum of monomer lengths in region */
    size_t min_monomer_length;
    size_t max_monomer_length;
} ConservedRegion;

/* Free EdArray */
static void free_ed_array(EdArray *arr) {
    if (arr) {
        free(arr->sequence_id);
        free(arr->orientation);
        free(arr->pairs);
        free(arr->result_regions);
        free(arr->result_stats);
    }
}

/* Parse header line: >seq_id\torientation=X\tn_monomers=N */
static int parse_header(const char *line, EdArray *arr) {
    if (line[0] != '>') return -1;

    char *copy = strdup(line + 1);  /* Skip '>' */
    size_t len = strlen(copy);
    while (len > 0 && (copy[len-1] == '\n' || copy[len-1] == '\r')) {
        copy[--len] = '\0';
    }

    char *saveptr;
    char *token;

    /* sequence_id */
    token = strtok_r(copy, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    arr->sequence_id = strdup(token);

    /* orientation=X */
    token = strtok_r(NULL, "\t", &saveptr);
    if (token && strncmp(token, "orientation=", 12) == 0) {
        arr->orientation = strdup(token + 12);
    } else {
        arr->orientation = strdup("unknown");
    }

    /* n_monomers=N */
    token = strtok_r(NULL, "\t", &saveptr);
    if (token && strncmp(token, "n_monomers=", 11) == 0) {
        arr->n_monomers = (size_t)atol(token + 11);
    }

    free(copy);
    return 0;
}

/* Parse data line: idx1->idx2\tlen1\tlen2\tedit_dist\tsimilarity%\ttype1\ttype2 */
static int parse_pair(const char *line, EdPair *pair) {
    char *copy = strdup(line);
    size_t len = strlen(copy);
    while (len > 0 && (copy[len-1] == '\n' || copy[len-1] == '\r')) {
        copy[--len] = '\0';
    }

    /* Parse idx1->idx2 */
    char *arrow = strstr(copy, "->");
    if (!arrow) { free(copy); return -1; }

    *arrow = '\0';
    pair->idx1 = atoi(copy);

    char *rest = arrow + 2;
    char *saveptr;
    char *token;

    token = strtok_r(rest, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    pair->idx2 = atoi(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    pair->len1 = (size_t)atol(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    pair->len2 = (size_t)atol(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    pair->edit_dist = (size_t)atol(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(copy); return -1; }
    pair->similarity = atof(token);  /* Parses "95.23%" as 95.23 */

    token = strtok_r(NULL, "\t", &saveptr);
    if (token) strncpy(pair->type1, token, 31);

    token = strtok_r(NULL, "\t", &saveptr);
    if (token) strncpy(pair->type2, token, 31);

    free(copy);
    return 0;
}

/* Add pair to array */
static int add_pair(EdArray *arr, EdPair *pair) {
    if (arr->pair_count >= arr->pair_capacity) {
        size_t new_cap = arr->pair_capacity == 0 ? 100 : arr->pair_capacity * 2;
        EdPair *new_pairs = realloc(arr->pairs, new_cap * sizeof(EdPair));
        if (!new_pairs) return -1;
        arr->pairs = new_pairs;
        arr->pair_capacity = new_cap;
    }
    arr->pairs[arr->pair_count++] = *pair;
    return 0;
}

/* Find conserved regions in array */
static void find_conserved_regions(EdArray *arr, double threshold, int min_length) {
    if (arr->pair_count == 0) {
        arr->result_regions = NULL;
        arr->result_regions_len = 0;
        arr->result_stats = NULL;
        arr->result_stats_len = 0;
        return;
    }

    /* Find regions */
    ConservedRegion *regions = NULL;
    size_t region_count = 0;
    size_t region_capacity = 0;

    int in_region = 0;
    int region_start = -1;
    double region_sum = 0;
    double region_min = 100.0;
    int region_len = 0;
    size_t region_mono_len = 0;
    size_t region_min_mono = SIZE_MAX;
    size_t region_max_mono = 0;

    /* Count conserved pairs for stats */
    size_t conserved_pairs = 0;
    size_t monomer_pairs = 0;  /* Only MONOMER-MONOMER pairs */
    size_t conserved_monomer_pairs = 0;

    for (size_t i = 0; i < arr->pair_count; i++) {
        EdPair *p = &arr->pairs[i];
        int is_conserved = (p->similarity >= threshold);
        int is_monomer_pair = (strcmp(p->type1, "MONOMER") == 0 &&
                               strcmp(p->type2, "MONOMER") == 0);

        if (is_conserved) conserved_pairs++;
        if (is_monomer_pair) {
            monomer_pairs++;
            if (is_conserved) conserved_monomer_pairs++;
        }

        if (is_conserved && is_monomer_pair) {
            if (!in_region) {
                /* Start new region */
                in_region = 1;
                region_start = p->idx1;
                region_sum = p->similarity;
                region_min = p->similarity;
                region_len = 1;
                region_mono_len = p->len1 + p->len2;
                region_min_mono = p->len1 < p->len2 ? p->len1 : p->len2;
                region_max_mono = p->len1 > p->len2 ? p->len1 : p->len2;
            } else {
                /* Extend region */
                region_sum += p->similarity;
                if (p->similarity < region_min) region_min = p->similarity;
                region_len++;
                region_mono_len += p->len2;  /* Add only second monomer (first already counted) */
                if (p->len2 < region_min_mono) region_min_mono = p->len2;
                if (p->len2 > region_max_mono) region_max_mono = p->len2;
            }
        } else {
            if (in_region && region_len >= min_length) {
                /* Save region */
                if (region_count >= region_capacity) {
                    size_t new_cap = region_capacity == 0 ? 100 : region_capacity * 2;
                    ConservedRegion *new_regions = realloc(regions, new_cap * sizeof(ConservedRegion));
                    if (new_regions) {
                        regions = new_regions;
                        region_capacity = new_cap;
                    }
                }
                if (region_count < region_capacity) {
                    regions[region_count].start_idx = region_start;
                    regions[region_count].end_idx = arr->pairs[i-1].idx2;
                    regions[region_count].length = region_len;
                    regions[region_count].avg_similarity = region_sum / region_len;
                    regions[region_count].min_similarity = region_min;
                    regions[region_count].total_monomer_length = region_mono_len;
                    regions[region_count].min_monomer_length = region_min_mono;
                    regions[region_count].max_monomer_length = region_max_mono;
                    region_count++;
                }
            }
            in_region = 0;
            region_len = 0;
            region_mono_len = 0;
            region_min_mono = SIZE_MAX;
            region_max_mono = 0;
        }
    }

    /* Handle last region */
    if (in_region && region_len >= min_length) {
        if (region_count >= region_capacity) {
            size_t new_cap = region_capacity == 0 ? 100 : region_capacity * 2;
            ConservedRegion *new_regions = realloc(regions, new_cap * sizeof(ConservedRegion));
            if (new_regions) {
                regions = new_regions;
                region_capacity = new_cap;
            }
        }
        if (region_count < region_capacity) {
            regions[region_count].start_idx = region_start;
            regions[region_count].end_idx = arr->pairs[arr->pair_count-1].idx2;
            regions[region_count].length = region_len;
            regions[region_count].avg_similarity = region_sum / region_len;
            regions[region_count].min_similarity = region_min;
            regions[region_count].total_monomer_length = region_mono_len;
            regions[region_count].min_monomer_length = region_min_mono;
            regions[region_count].max_monomer_length = region_max_mono;
            region_count++;
        }
    }

    /* Build regions output */
    size_t buf_size = 256 + region_count * 300;
    arr->result_regions = malloc(buf_size);
    if (arr->result_regions) {
        char *ptr = arr->result_regions;
        size_t remaining = buf_size;

        for (size_t i = 0; i < region_count; i++) {
            ConservedRegion *r = &regions[i];
            int n_monomers = r->end_idx - r->start_idx + 1;
            double avg_mono_len = n_monomers > 0 ? (double)r->total_monomer_length / n_monomers : 0;
            int written = snprintf(ptr, remaining,
                "%s\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%zu\t%zu\t%zu\t%.1f\n",
                arr->sequence_id, r->start_idx, r->end_idx,
                r->length, n_monomers,
                r->avg_similarity, r->min_similarity,
                r->total_monomer_length, r->min_monomer_length, r->max_monomer_length,
                avg_mono_len);
            ptr += written;
            remaining -= written;
        }
        arr->result_regions_len = ptr - arr->result_regions;
    }

    /* Build stats output */
    arr->result_stats = malloc(512);
    if (arr->result_stats) {
        double conserved_pct = arr->pair_count > 0 ?
            100.0 * conserved_pairs / arr->pair_count : 0;
        double monomer_conserved_pct = monomer_pairs > 0 ?
            100.0 * conserved_monomer_pairs / monomer_pairs : 0;

        /* Calculate total length of conserved regions */
        int total_conserved_len = 0;
        for (size_t i = 0; i < region_count; i++) {
            total_conserved_len += regions[i].length;
        }

        arr->result_stats_len = snprintf(arr->result_stats, 512,
            "%s\t%zu\t%zu\t%.2f%%\t%zu\t%zu\t%.2f%%\t%zu\t%d\n",
            arr->sequence_id,
            arr->pair_count, conserved_pairs, conserved_pct,
            monomer_pairs, conserved_monomer_pairs, monomer_conserved_pct,
            region_count, total_conserved_len);
    }

    free(regions);
}

/* Worker thread */
static void* worker_thread(void *arg) {
    WorkQueue *queue = (WorkQueue*)arg;

    while (1) {
        pthread_mutex_lock(&queue->mutex);
        if (queue->next_item >= queue->count) {
            pthread_mutex_unlock(&queue->mutex);
            break;
        }
        size_t idx = queue->next_item++;
        pthread_mutex_unlock(&queue->mutex);

        find_conserved_regions(&queue->arrays[idx],
                               queue->threshold,
                               queue->min_region_length);

        pthread_mutex_lock(&queue->mutex);
        queue->completed++;
        size_t done = queue->completed;
        size_t total = queue->count;
        pthread_mutex_unlock(&queue->mutex);

        if (done % 100 == 0 || done == total) {
            fprintf(stderr, "\rAnalyzing: %zu / %zu arrays (%.1f%%)   ",
                    done, total, 100.0 * done / total);
            fflush(stderr);
        }
    }

    return NULL;
}

/* Print usage */
static void print_usage(const char *prog) {
    fprintf(stderr, "Conserved Region Finder (Multi-threaded)\n\n");
    fprintf(stderr, "Usage: %s -i input.ed -o output_prefix [options]\n\n", prog);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -i, --input FILE     Input .ed file (from arraysplitter_edit_dist_c)\n");
    fprintf(stderr, "  -o, --output PREFIX  Output prefix for results\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -s, --similarity N   Minimum similarity threshold (default: 95.0)\n");
    fprintf(stderr, "  -m, --min-length N   Minimum region length in pairs (default: 3)\n");
    fprintf(stderr, "  -t, --threads N      Number of threads (default: 4)\n");
    fprintf(stderr, "  -h, --help           Show this help\n");
}

int main(int argc, char *argv[]) {
    const char *input_file = NULL;
    const char *output_prefix = NULL;
    double threshold = 95.0;
    int min_length = 3;
    int num_threads = 4;

    static struct option long_options[] = {
        {"input",      required_argument, 0, 'i'},
        {"output",     required_argument, 0, 'o'},
        {"similarity", required_argument, 0, 's'},
        {"min-length", required_argument, 0, 'm'},
        {"threads",    required_argument, 0, 't'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:s:m:t:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i': input_file = optarg; break;
            case 'o': output_prefix = optarg; break;
            case 's': threshold = atof(optarg); break;
            case 'm': min_length = atoi(optarg); break;
            case 't': num_threads = atoi(optarg); if (num_threads < 1) num_threads = 1; break;
            case 'h': print_usage(argv[0]); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    if (!input_file || !output_prefix) {
        fprintf(stderr, "Error: input file and output prefix are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    FILE *in = fopen(input_file, "r");
    if (!in) {
        fprintf(stderr, "Error: cannot open input file '%s'\n", input_file);
        return 1;
    }

    printf("Conserved Region Finder\n");
    printf("Input: %s\n", input_file);
    printf("Output prefix: %s\n", output_prefix);
    printf("Similarity threshold: %.1f%%\n", threshold);
    printf("Minimum region length: %d pairs\n", min_length);
    printf("Threads: %d\n", num_threads);
    printf("\n");

    /* Phase 1: Read data */
    printf("Phase 1: Reading input...\n");

    char *line = malloc(MAX_LINE_LENGTH);
    EdArray *all_arrays = NULL;
    size_t arrays_count = 0;
    size_t arrays_capacity = 0;
    EdArray current = {0};
    size_t line_num = 0;

    while (fgets(line, MAX_LINE_LENGTH, in)) {
        line_num++;

        /* Skip comments and empty lines */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }

        if (line[0] == '>') {
            /* Save previous array */
            if (current.sequence_id) {
                if (arrays_count >= arrays_capacity) {
                    size_t new_cap = arrays_capacity == 0 ? 1000 : arrays_capacity * 2;
                    EdArray *new_arrays = realloc(all_arrays, new_cap * sizeof(EdArray));
                    if (!new_arrays) break;
                    all_arrays = new_arrays;
                    arrays_capacity = new_cap;
                }
                all_arrays[arrays_count++] = current;
                memset(&current, 0, sizeof(current));
            }

            parse_header(line, &current);
        } else {
            /* Data line */
            EdPair pair = {0};
            if (parse_pair(line, &pair) == 0) {
                add_pair(&current, &pair);
            }
        }

        if (line_num % 100000 == 0) {
            fprintf(stderr, "\rPhase 1: %zu lines, %zu arrays...   ", line_num, arrays_count);
            fflush(stderr);
        }
    }

    /* Save last array */
    if (current.sequence_id) {
        if (arrays_count >= arrays_capacity) {
            size_t new_cap = arrays_capacity == 0 ? 1000 : arrays_capacity * 2;
            EdArray *new_arrays = realloc(all_arrays, new_cap * sizeof(EdArray));
            if (new_arrays) {
                all_arrays = new_arrays;
                arrays_capacity = new_cap;
            }
        }
        if (arrays_count < arrays_capacity) {
            all_arrays[arrays_count++] = current;
        }
    }

    free(line);
    fclose(in);

    fprintf(stderr, "\rPhase 1: %zu lines, %zu arrays.          \n", line_num, arrays_count);

    if (arrays_count == 0) {
        fprintf(stderr, "No arrays found\n");
        return 1;
    }

    /* Phase 2: Find conserved regions */
    printf("\nPhase 2: Finding conserved regions (%d threads)...\n", num_threads);

    WorkQueue queue = {
        .arrays = all_arrays,
        .count = arrays_count,
        .next_item = 0,
        .mutex = PTHREAD_MUTEX_INITIALIZER,
        .completed = 0,
        .threshold = threshold,
        .min_region_length = min_length
    };

    if (num_threads > (int)arrays_count) {
        num_threads = (int)arrays_count;
    }

    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker_thread, &queue);
    }
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    fprintf(stderr, "\rAnalyzing: %zu / %zu arrays (100.0%%)          \n", arrays_count, arrays_count);

    free(threads);
    pthread_mutex_destroy(&queue.mutex);

    /* Phase 3: Write output */
    printf("\nPhase 3: Writing results...\n");

    char regions_file[1024];
    char stats_file[1024];
    snprintf(regions_file, sizeof(regions_file), "%s.regions.tsv", output_prefix);
    snprintf(stats_file, sizeof(stats_file), "%s.stats.tsv", output_prefix);

    FILE *f_regions = fopen(regions_file, "w");
    FILE *f_stats = fopen(stats_file, "w");

    if (!f_regions || !f_stats) {
        fprintf(stderr, "Error: cannot open output files\n");
        return 1;
    }

    /* Write headers */
    fprintf(f_regions, "# Conserved regions (similarity >= %.1f%%, min_length >= %d)\n", threshold, min_length);
    fprintf(f_regions, "sequence_id\tstart_idx\tend_idx\tn_pairs\tn_monomers\tavg_similarity\tmin_similarity\ttotal_bp\tmin_mono_len\tmax_mono_len\tavg_mono_len\n");

    fprintf(f_stats, "# Statistics per array (threshold: %.1f%%)\n", threshold);
    fprintf(f_stats, "sequence_id\ttotal_pairs\tconserved_pairs\tconserved_pct\tmonomer_pairs\tconserved_monomer_pairs\tmonomer_conserved_pct\tnum_regions\ttotal_region_length\n");

    size_t total_regions = 0;
    for (size_t i = 0; i < arrays_count; i++) {
        if (all_arrays[i].result_regions && all_arrays[i].result_regions_len > 0) {
            fwrite(all_arrays[i].result_regions, 1, all_arrays[i].result_regions_len, f_regions);
        }
        if (all_arrays[i].result_stats && all_arrays[i].result_stats_len > 0) {
            fwrite(all_arrays[i].result_stats, 1, all_arrays[i].result_stats_len, f_stats);
        }

        /* Count regions */
        if (all_arrays[i].result_regions) {
            char *p = all_arrays[i].result_regions;
            while (*p) {
                if (*p == '\n') total_regions++;
                p++;
            }
        }
    }

    fclose(f_regions);
    fclose(f_stats);

    /* Cleanup */
    for (size_t i = 0; i < arrays_count; i++) {
        free_ed_array(&all_arrays[i]);
    }
    free(all_arrays);

    printf("\nResults:\n");
    printf("  Regions file: %s\n", regions_file);
    printf("  Stats file: %s\n", stats_file);
    printf("  Total conserved regions found: %zu\n", total_regions);
    printf("\nDone!\n");

    return 0;
}
