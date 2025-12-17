/*
 * Edit Distance Calculator for monomers TSV (Multi-threaded)
 *
 * Reads .monomers.tsv file and calculates edit distances between
 * consecutive monomers within each array.
 *
 * Usage: arraysplitter_edit_dist -i monomers.tsv -o output.ed [-t threads]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <pthread.h>

#define MAX_LINE_LENGTH 100000
#define MAX_SEQ_LENGTH 100000
#define MAX_ID_LENGTH 4096

/* External edit distance function */
extern size_t seq_edit_distance(const char *s1, size_t len1, const char *s2, size_t len2);

/* Monomer entry from TSV */
typedef struct {
    char *sequence_id;
    char *orientation;
    int index;
    char *type;
    size_t length;
    int is_flank;
    char *sequence;
} MonomerEntry;

/* Array to hold all monomers for one sequence */
typedef struct {
    MonomerEntry *entries;
    size_t count;
    size_t capacity;
    char *result_buffer;  /* Pre-allocated buffer for results */
    size_t result_len;
} MonomerArray;

/* Thread work item */
typedef struct {
    MonomerArray *array;
    int only_monomers;
} WorkItem;

/* Thread pool data */
typedef struct {
    WorkItem *items;
    size_t count;
    size_t next_item;
    pthread_mutex_t mutex;
    size_t completed;
} WorkQueue;

/* Free monomer entry */
static void free_entry(MonomerEntry *e) {
    if (e) {
        free(e->sequence_id);
        free(e->orientation);
        free(e->type);
        free(e->sequence);
    }
}

/* Free monomer array */
static void free_monomer_array(MonomerArray *arr) {
    if (arr) {
        for (size_t i = 0; i < arr->count; i++) {
            free_entry(&arr->entries[i]);
        }
        free(arr->entries);
        free(arr->result_buffer);
        arr->entries = NULL;
        arr->result_buffer = NULL;
        arr->count = 0;
        arr->capacity = 0;
    }
}

/* Parse TSV line */
static int parse_line(const char *line, MonomerEntry *entry) {
    char *line_copy = strdup(line);
    if (!line_copy) return -1;

    size_t len = strlen(line_copy);
    while (len > 0 && (line_copy[len-1] == '\n' || line_copy[len-1] == '\r')) {
        line_copy[--len] = '\0';
    }

    char *saveptr;
    char *token;

    token = strtok_r(line_copy, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->sequence_id = strdup(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->orientation = strdup(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->index = atoi(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->type = strdup(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->length = (size_t)atol(token);

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->is_flank = (strcmp(token, "TRUE") == 0) ? 1 : 0;

    token = strtok_r(NULL, "\t", &saveptr);
    if (!token) { free(line_copy); return -1; }
    entry->sequence = strdup(token);

    free(line_copy);
    return 0;
}

/* Add entry to array */
static int add_entry(MonomerArray *arr, MonomerEntry *entry) {
    if (arr->count >= arr->capacity) {
        size_t new_cap = arr->capacity == 0 ? 100 : arr->capacity * 2;
        MonomerEntry *new_entries = realloc(arr->entries, new_cap * sizeof(MonomerEntry));
        if (!new_entries) return -1;
        arr->entries = new_entries;
        arr->capacity = new_cap;
    }
    arr->entries[arr->count++] = *entry;
    return 0;
}

/* Process single array (thread-safe, writes to pre-allocated buffer) */
static void process_array(MonomerArray *arr, int only_monomers) {
    if (arr->count < 2) {
        arr->result_buffer = NULL;
        arr->result_len = 0;
        return;
    }

    /* Estimate buffer size: ~200 bytes per pair + header */
    size_t buf_size = 512 + arr->count * 200;
    arr->result_buffer = malloc(buf_size);
    if (!arr->result_buffer) {
        arr->result_len = 0;
        return;
    }

    char *ptr = arr->result_buffer;
    size_t remaining = buf_size;

    /* Write header */
    int written = snprintf(ptr, remaining, ">%s\torientation=%s\tn_monomers=%zu\n",
            arr->entries[0].sequence_id,
            arr->entries[0].orientation,
            arr->count);
    ptr += written;
    remaining -= written;

    /* Calculate edit distances */
    for (size_t i = 0; i < arr->count - 1; i++) {
        MonomerEntry *m1 = &arr->entries[i];
        MonomerEntry *m2 = &arr->entries[i + 1];

        if (only_monomers && (m1->is_flank || m2->is_flank)) {
            continue;
        }

        size_t edit_dist = seq_edit_distance(
            m1->sequence, m1->length,
            m2->sequence, m2->length
        );

        size_t max_len = m1->length > m2->length ? m1->length : m2->length;
        double similarity = max_len > 0 ? (1.0 - (double)edit_dist / max_len) * 100.0 : 100.0;

        written = snprintf(ptr, remaining, "%d->%d\t%zu\t%zu\t%zu\t%.2f%%\t%s\t%s\n",
                m1->index, m2->index,
                m1->length, m2->length,
                edit_dist,
                similarity,
                m1->type, m2->type);

        if (written >= (int)remaining) {
            /* Buffer overflow, reallocate */
            size_t current_len = ptr - arr->result_buffer;
            buf_size *= 2;
            char *new_buf = realloc(arr->result_buffer, buf_size);
            if (!new_buf) break;
            arr->result_buffer = new_buf;
            ptr = arr->result_buffer + current_len;
            remaining = buf_size - current_len;

            /* Retry write */
            written = snprintf(ptr, remaining, "%d->%d\t%zu\t%zu\t%zu\t%.2f%%\t%s\t%s\n",
                    m1->index, m2->index,
                    m1->length, m2->length,
                    edit_dist,
                    similarity,
                    m1->type, m2->type);
        }

        ptr += written;
        remaining -= written;
    }

    arr->result_len = ptr - arr->result_buffer;
}

/* Worker thread function */
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

        /* Process array */
        WorkItem *item = &queue->items[idx];
        process_array(item->array, item->only_monomers);

        /* Update progress */
        pthread_mutex_lock(&queue->mutex);
        queue->completed++;
        size_t done = queue->completed;
        size_t total = queue->count;
        pthread_mutex_unlock(&queue->mutex);

        if (done % 100 == 0 || done == total) {
            fprintf(stderr, "\rPhase 2: %zu / %zu arrays (%.1f%%)   ",
                    done, total, 100.0 * done / total);
            fflush(stderr);
        }
    }

    return NULL;
}

/* Print usage */
static void print_usage(const char *prog) {
    fprintf(stderr, "Edit Distance Calculator for monomers (Multi-threaded)\n\n");
    fprintf(stderr, "Usage: %s -i input.monomers.tsv -o output.ed [options]\n\n", prog);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -i, --input FILE     Input .monomers.tsv file\n");
    fprintf(stderr, "  -o, --output FILE    Output file for edit distances\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -t, --threads N      Number of threads (default: 4)\n");
    fprintf(stderr, "  -m, --monomers-only  Only compare monomers (skip flanks)\n");
    fprintf(stderr, "  -h, --help           Show this help\n");
}

int main(int argc, char *argv[]) {
    const char *input_file = NULL;
    const char *output_file = NULL;
    int num_threads = 4;
    int only_monomers = 0;

    static struct option long_options[] = {
        {"input",         required_argument, 0, 'i'},
        {"output",        required_argument, 0, 'o'},
        {"threads",       required_argument, 0, 't'},
        {"monomers-only", no_argument,       0, 'm'},
        {"help",          no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:t:mh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 't':
                num_threads = atoi(optarg);
                if (num_threads < 1) num_threads = 1;
                break;
            case 'm':
                only_monomers = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (!input_file || !output_file) {
        fprintf(stderr, "Error: input and output files are required\n\n");
        print_usage(argv[0]);
        return 1;
    }

    FILE *in = fopen(input_file, "r");
    if (!in) {
        fprintf(stderr, "Error: cannot open input file '%s'\n", input_file);
        return 1;
    }

    printf("Edit Distance Calculator (Multi-threaded)\n");
    printf("Input: %s\n", input_file);
    printf("Output: %s\n", output_file);
    printf("Threads: %d\n", num_threads);
    if (only_monomers) {
        printf("Mode: monomers only (skipping flanks)\n");
    }
    printf("\n");

    /* Phase 1: Read all data into memory */
    printf("Phase 1: Reading input file...\n");

    char *line = malloc(MAX_LINE_LENGTH);
    if (!line) {
        fprintf(stderr, "Error: out of memory\n");
        fclose(in);
        return 1;
    }

    /* Array of all sequence arrays */
    MonomerArray *all_arrays = NULL;
    size_t arrays_count = 0;
    size_t arrays_capacity = 0;

    MonomerArray current_array = {0};
    char *current_seq_id = NULL;
    size_t line_num = 0;

    while (fgets(line, MAX_LINE_LENGTH, in)) {
        line_num++;

        if (line_num == 1 && strncmp(line, "sequence_id", 11) == 0) {
            continue;
        }

        MonomerEntry entry = {0};
        if (parse_line(line, &entry) != 0) {
            continue;
        }

        if (current_seq_id == NULL || strcmp(current_seq_id, entry.sequence_id) != 0) {
            if (current_array.count > 0) {
                /* Save current array */
                if (arrays_count >= arrays_capacity) {
                    size_t new_cap = arrays_capacity == 0 ? 10000 : arrays_capacity * 2;
                    MonomerArray *new_arrays = realloc(all_arrays, new_cap * sizeof(MonomerArray));
                    if (!new_arrays) {
                        fprintf(stderr, "Error: out of memory\n");
                        break;
                    }
                    all_arrays = new_arrays;
                    arrays_capacity = new_cap;
                }
                all_arrays[arrays_count++] = current_array;
                memset(&current_array, 0, sizeof(current_array));
            }

            free(current_seq_id);
            current_seq_id = strdup(entry.sequence_id);
        }

        add_entry(&current_array, &entry);

        if (line_num % 50000 == 0) {
            fprintf(stderr, "\rPhase 1: %zu lines, %zu arrays...   ", line_num, arrays_count);
            fflush(stderr);
        }
    }

    /* Save last array */
    if (current_array.count > 0) {
        if (arrays_count >= arrays_capacity) {
            size_t new_cap = arrays_capacity == 0 ? 1000 : arrays_capacity * 2;
            MonomerArray *new_arrays = realloc(all_arrays, new_cap * sizeof(MonomerArray));
            if (new_arrays) {
                all_arrays = new_arrays;
                arrays_capacity = new_cap;
            }
        }
        if (arrays_count < arrays_capacity) {
            all_arrays[arrays_count++] = current_array;
        }
    }

    free(current_seq_id);
    free(line);
    fclose(in);

    fprintf(stderr, "\rRead %zu lines, %zu arrays.          \n", line_num, arrays_count);

    if (arrays_count == 0) {
        fprintf(stderr, "No arrays found in input file\n");
        return 1;
    }

    /* Phase 2: Process arrays in parallel */
    printf("\nPhase 2: Computing edit distances (%d threads)...\n", num_threads);

    /* Create work items */
    WorkItem *work_items = malloc(arrays_count * sizeof(WorkItem));
    for (size_t i = 0; i < arrays_count; i++) {
        work_items[i].array = &all_arrays[i];
        work_items[i].only_monomers = only_monomers;
    }

    /* Create work queue */
    WorkQueue queue = {
        .items = work_items,
        .count = arrays_count,
        .next_item = 0,
        .mutex = PTHREAD_MUTEX_INITIALIZER,
        .completed = 0
    };

    /* Limit threads */
    if (num_threads > (int)arrays_count) {
        num_threads = (int)arrays_count;
    }

    /* Create threads */
    pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
    for (int i = 0; i < num_threads; i++) {
        pthread_create(&threads[i], NULL, worker_thread, &queue);
    }

    /* Wait for completion */
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    fprintf(stderr, "\rProcessed %zu / %zu arrays (100.0%%)          \n", arrays_count, arrays_count);

    free(threads);
    free(work_items);
    pthread_mutex_destroy(&queue.mutex);

    /* Phase 3: Write results */
    printf("\nPhase 3: Writing output...\n");

    FILE *out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Error: cannot open output file '%s'\n", output_file);
        return 1;
    }

    fprintf(out, "# Edit distances between consecutive monomers\n");
    fprintf(out, "# Format: idx1->idx2\\tlength1\\tlength2\\tedit_dist\\tsimilarity%%\\ttype1\\ttype2\n\n");

    for (size_t i = 0; i < arrays_count; i++) {
        if (all_arrays[i].result_buffer && all_arrays[i].result_len > 0) {
            fwrite(all_arrays[i].result_buffer, 1, all_arrays[i].result_len, out);
        }
    }

    fclose(out);

    /* Cleanup */
    for (size_t i = 0; i < arrays_count; i++) {
        free_monomer_array(&all_arrays[i]);
    }
    free(all_arrays);

    printf("\nDone! Output written to: %s\n", output_file);

    return 0;
}
