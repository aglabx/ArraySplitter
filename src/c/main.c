/*
 * ArraySplitter C - Main entry point with multithreading
 *
 * Usage: arraysplitter_c -i input.fa -o output_prefix [-t threads] [-d depth] [-v]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>
#include "arraysplitter.h"

/* Thread-safe progress counter */
typedef struct {
    size_t completed;
    size_t total;
    pthread_mutex_t mutex;
} ProgressCounter;

/* Thread work data */
typedef struct {
    FastaArray *arrays;
    DecompResult **results;
    size_t *work_queue;
    size_t queue_size;
    size_t next_work;
    pthread_mutex_t queue_mutex;
    int depth;
    bool verbose;
    ProgressCounter *progress;
} ThreadData;

/* Worker thread function */
static void* worker_thread(void *arg) {
    ThreadData *data = (ThreadData*)arg;

    while (1) {
        /* Get next work item */
        pthread_mutex_lock(&data->queue_mutex);
        if (data->next_work >= data->queue_size) {
            pthread_mutex_unlock(&data->queue_mutex);
            break;
        }
        size_t idx = data->work_queue[data->next_work++];
        pthread_mutex_unlock(&data->queue_mutex);

        /* Process array */
        FastaEntry *entry = &data->arrays->entries[idx];
        DecompResult *result = decompose_full(entry->header, entry->sequence,
                                               entry->length, data->depth,
                                               data->verbose);
        data->results[idx] = result;

        /* Update progress */
        pthread_mutex_lock(&data->progress->mutex);
        data->progress->completed++;
        size_t done = data->progress->completed;
        size_t total = data->progress->total;
        pthread_mutex_unlock(&data->progress->mutex);

        /* Print progress */
        if (!data->verbose && done % 100 == 0) {
            fprintf(stderr, "\rProcessed %zu / %zu arrays (%.1f%%)",
                    done, total, 100.0 * done / total);
        }
    }

    return NULL;
}

/* Process arrays in parallel */
DecompResult** process_parallel(FastaArray *arrays, int num_threads,
                                 int depth, bool verbose) {
    if (!arrays || arrays->count == 0) return NULL;

    /* Allocate results array */
    DecompResult **results = calloc(arrays->count, sizeof(DecompResult*));
    if (!results) return NULL;

    /* Create work queue (indices of arrays to process) */
    size_t *work_queue = malloc(arrays->count * sizeof(size_t));
    for (size_t i = 0; i < arrays->count; i++) {
        work_queue[i] = i;
    }

    /* Initialize progress counter */
    ProgressCounter progress = {0, arrays->count, PTHREAD_MUTEX_INITIALIZER};

    /* Initialize thread data */
    ThreadData data = {
        .arrays = arrays,
        .results = results,
        .work_queue = work_queue,
        .queue_size = arrays->count,
        .next_work = 0,
        .queue_mutex = PTHREAD_MUTEX_INITIALIZER,
        .depth = depth,
        .verbose = verbose,
        .progress = &progress
    };

    /* Limit threads to number of arrays */
    if (num_threads > (int)arrays->count) {
        num_threads = (int)arrays->count;
    }

    /* Create threads */
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));

    fprintf(stderr, "Starting %d threads for %zu arrays...\n",
            num_threads, arrays->count);

    for (int i = 0; i < num_threads; i++) {
        if (pthread_create(&threads[i], NULL, worker_thread, &data) != 0) {
            fprintf(stderr, "Error creating thread %d\n", i);
            num_threads = i;
            break;
        }
    }

    /* Wait for all threads */
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    fprintf(stderr, "\rProcessed %zu / %zu arrays (100.0%%)\n",
            arrays->count, arrays->count);

    /* Cleanup */
    free(threads);
    free(work_queue);
    pthread_mutex_destroy(&data.queue_mutex);
    pthread_mutex_destroy(&progress.mutex);

    return results;
}

/* Write results to output files */
int write_results(const char *output_prefix, FastaArray *arrays,
                  DecompResult **results, size_t count) {
    /* Create output filenames */
    size_t prefix_len = strlen(output_prefix);
    char *fasta_file = malloc(prefix_len + 20);
    char *tsv_file = malloc(prefix_len + 20);
    char *lengths_file = malloc(prefix_len + 20);

    sprintf(fasta_file, "%s.decomposed.fasta", output_prefix);
    sprintf(tsv_file, "%s.monomers.tsv", output_prefix);
    sprintf(lengths_file, "%s.lengths", output_prefix);

    /* Open files */
    FILE *f_fasta = fopen(fasta_file, "w");
    FILE *f_tsv = fopen(tsv_file, "w");
    FILE *f_lengths = fopen(lengths_file, "w");

    if (!f_fasta || !f_tsv || !f_lengths) {
        fprintf(stderr, "Error opening output files\n");
        free(fasta_file);
        free(tsv_file);
        free(lengths_file);
        return -1;
    }

    /* Write TSV header */
    fprintf(f_tsv, "sequence_id\torientation\tindex\ttype\tlength\tis_flank\tsequence\n");

    /* Write results */
    for (size_t i = 0; i < count; i++) {
        FastaEntry *entry = &arrays->entries[i];
        DecompResult *result = results[i];

        if (!result) continue;

        const char *orientation = result->was_reversed ? "rev" : "fwd";

        /* Calculate statistics for internal monomers */
        size_t internal_count = 0;
        size_t min_len = SIZE_MAX, max_len = 0;
        double sum_len = 0;

        for (size_t j = 0; j < result->monomer_count; j++) {
            const char *mono = result->monomers[j];
            size_t mono_len = result->monomer_lengths[j];

            /* Check if it starts with cut sequence */
            if (result->cut_seq && strlen(result->cut_seq) > 0 &&
                strncmp(mono, result->cut_seq, strlen(result->cut_seq)) == 0) {
                internal_count++;
                if (mono_len < min_len) min_len = mono_len;
                if (mono_len > max_len) max_len = mono_len;
                sum_len += mono_len;
            }
        }

        double avg_len = internal_count > 0 ? sum_len / internal_count : 0;

        /* Write FASTA header */
        fprintf(f_fasta, ">%s cut=%s orientation=%s n_monomers=%zu range=%zu-%zu avg=%.1f\n",
                entry->header, result->cut_seq, orientation,
                internal_count, min_len == SIZE_MAX ? 0 : min_len, max_len, avg_len);

        /* Write monomers space-separated */
        for (size_t j = 0; j < result->monomer_count; j++) {
            if (j > 0) fprintf(f_fasta, " ");
            fprintf(f_fasta, "%s", result->monomers[j]);
        }
        fprintf(f_fasta, "\n");

        /* Write lengths file */
        fprintf(f_lengths, ">%s cut=%s orientation=%s n_monomers=%zu range=%zu-%zu avg=%.1f\n",
                entry->header, result->cut_seq, orientation,
                internal_count, min_len == SIZE_MAX ? 0 : min_len, max_len, avg_len);

        for (size_t j = 0; j < result->monomer_count; j++) {
            if (j > 0) fprintf(f_lengths, " ");
            fprintf(f_lengths, "%zu", result->monomer_lengths[j]);
        }
        fprintf(f_lengths, "\n");

        /* Write TSV details */
        double flank_threshold = avg_len * 0.7;

        for (size_t j = 0; j < result->monomer_count; j++) {
            const char *mono = result->monomers[j];
            size_t mono_len = result->monomer_lengths[j];

            const char *piece_type;
            const char *is_flank;

            if (j == 0 && (strlen(result->cut_seq) == 0 ||
                           strncmp(mono, result->cut_seq, strlen(result->cut_seq)) != 0)) {
                piece_type = "LEFT_FLANK";
                is_flank = "TRUE";
            } else if (j == result->monomer_count - 1 && mono_len < flank_threshold) {
                piece_type = "RIGHT_FLANK";
                is_flank = "TRUE";
            } else {
                piece_type = "MONOMER";
                is_flank = "FALSE";
            }

            fprintf(f_tsv, "%s\t%s\t%zu\t%s\t%zu\t%s\t%s\n",
                    entry->header, orientation, j, piece_type,
                    mono_len, is_flank, mono);
        }
    }

    fclose(f_fasta);
    fclose(f_tsv);
    fclose(f_lengths);

    printf("Output file: %s\n", fasta_file);
    printf("Detail file: %s\n", tsv_file);
    printf("Lengths file: %s\n", lengths_file);

    free(fasta_file);
    free(tsv_file);
    free(lengths_file);

    return 0;
}

/* Print usage */
static void print_usage(const char *prog) {
    fprintf(stderr, "ArraySplitter C - De novo decomposition of satellite DNA arrays\n\n");
    fprintf(stderr, "Usage: %s -i input.fa -o output_prefix [options]\n\n", prog);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -i, --input FILE     Input FASTA file\n");
    fprintf(stderr, "  -o, --output PREFIX  Output prefix\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -t, --threads N      Number of threads (default: 4)\n");
    fprintf(stderr, "  -d, --depth N        Depth for hint discovery (default: 100)\n");
    fprintf(stderr, "  -v, --verbose        Verbose output\n");
    fprintf(stderr, "  -h, --help           Show this help\n");
}

int main(int argc, char *argv[]) {
    /* Default options */
    const char *input_file = NULL;
    const char *output_prefix = NULL;
    int num_threads = 4;
    int depth = 100;
    bool verbose = false;

    /* Parse arguments */
    static struct option long_options[] = {
        {"input",   required_argument, 0, 'i'},
        {"output",  required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"depth",   required_argument, 0, 'd'},
        {"verbose", no_argument,       0, 'v'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:t:d:vh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads < 1) num_threads = 1;
                break;
            case 'd':
                depth = atoi(optarg);
                if (depth < 1) depth = 100;
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    /* Check required arguments */
    if (!input_file || !output_prefix) {
        fprintf(stderr, "Error: input file and output prefix are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Check if input file exists */
    if (access(input_file, R_OK) != 0) {
        fprintf(stderr, "Error: Cannot read input file '%s'\n", input_file);
        return 1;
    }

    /* Remove .fasta/.fa suffix from output prefix if present */
    char *clean_prefix = strdup(output_prefix);
    size_t len = strlen(clean_prefix);
    if (len > 6 && strcmp(clean_prefix + len - 6, ".fasta") == 0) {
        clean_prefix[len - 6] = '\0';
        printf("Remove .fasta from output prefix\n");
    } else if (len > 3 && strcmp(clean_prefix + len - 3, ".fa") == 0) {
        clean_prefix[len - 3] = '\0';
        printf("Remove .fa from output prefix\n");
    }

    printf("ArraySplitter C v1.3.0\n");
    printf("Input: %s\n", input_file);
    printf("Output prefix: %s\n", clean_prefix);
    printf("Threads: %d\n", num_threads);
    printf("Depth: %d\n", depth);
    printf("\n");

    /* Read input file */
    printf("Reading input file...\n");
    FastaArray *arrays = fasta_read_file(input_file);
    if (!arrays) {
        fprintf(stderr, "Error: Failed to read input file\n");
        free(clean_prefix);
        return 1;
    }
    printf("Loaded %zu sequences\n\n", arrays->count);

    /* Process arrays */
    DecompResult **results = process_parallel(arrays, num_threads, depth, verbose);
    if (!results) {
        fprintf(stderr, "Error: Processing failed\n");
        fasta_free(arrays);
        free(clean_prefix);
        return 1;
    }

    /* Write output */
    printf("\nWriting output files...\n");
    int ret = write_results(clean_prefix, arrays, results, arrays->count);

    /* Cleanup */
    for (size_t i = 0; i < arrays->count; i++) {
        decomp_result_free(results[i]);
    }
    free(results);
    fasta_free(arrays);
    free(clean_prefix);

    printf("\nDone!\n");
    return ret;
}
