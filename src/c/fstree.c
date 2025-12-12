/*
 * Frequency Suffix Tree implementation
 *
 * Builds a suffix tree-like structure where each node tracks
 * positions of sequences that share a common prefix.
 * Used to find frequent subsequences (anchors/hints).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arraysplitter.h"

/* Priority queue node for BFS traversal */
typedef struct HeapNode {
    size_t level;
    size_t *names;       /* Starting positions */
    size_t *positions;   /* Current positions */
    size_t count;
} HeapNode;

/* Min-heap for BFS by level */
typedef struct {
    HeapNode *nodes;
    size_t count;
    size_t capacity;
} MinHeap;

static void heap_init(MinHeap *h, size_t capacity) {
    h->nodes = malloc(capacity * sizeof(HeapNode));
    h->count = 0;
    h->capacity = capacity;
}

static void heap_free(MinHeap *h) {
    for (size_t i = 0; i < h->count; i++) {
        free(h->nodes[i].names);
        free(h->nodes[i].positions);
    }
    free(h->nodes);
    h->nodes = NULL;
    h->count = 0;
}

static void heap_push(MinHeap *h, HeapNode node) {
    if (h->count >= h->capacity) {
        h->capacity *= 2;
        h->nodes = realloc(h->nodes, h->capacity * sizeof(HeapNode));
    }

    /* Insert at end */
    size_t i = h->count++;
    h->nodes[i] = node;

    /* Bubble up */
    while (i > 0) {
        size_t parent = (i - 1) / 2;
        if (h->nodes[parent].level <= h->nodes[i].level) break;

        HeapNode tmp = h->nodes[parent];
        h->nodes[parent] = h->nodes[i];
        h->nodes[i] = tmp;
        i = parent;
    }
}

static HeapNode heap_pop(MinHeap *h) {
    HeapNode result = h->nodes[0];

    /* Move last to root */
    h->nodes[0] = h->nodes[--h->count];

    /* Bubble down */
    size_t i = 0;
    while (1) {
        size_t left = 2 * i + 1;
        size_t right = 2 * i + 2;
        size_t smallest = i;

        if (left < h->count && h->nodes[left].level < h->nodes[smallest].level)
            smallest = left;
        if (right < h->count && h->nodes[right].level < h->nodes[smallest].level)
            smallest = right;

        if (smallest == i) break;

        HeapNode tmp = h->nodes[i];
        h->nodes[i] = h->nodes[smallest];
        h->nodes[smallest] = tmp;
        i = smallest;
    }

    return result;
}

/* Check if pattern is composed of repeated smaller units */
static char* is_self_repeating(const char *pattern, size_t len) {
    for (size_t sub_len = 1; sub_len <= len / 2; sub_len++) {
        if (len % sub_len != 0) continue;

        bool is_repeat = true;
        for (size_t i = sub_len; i < len && is_repeat; i++) {
            if (pattern[i] != pattern[i % sub_len]) {
                is_repeat = false;
            }
        }

        if (is_repeat) {
            char *sub = malloc(sub_len + 1);
            memcpy(sub, pattern, sub_len);
            sub[sub_len] = '\0';
            return sub;
        }
    }
    return NULL;
}

/* Add hint to array if not duplicate */
static void hints_add(HintArray *hints, const char *seq, size_t len, size_t freq) {
    /* Check for duplicates */
    for (size_t i = 0; i < hints->count; i++) {
        if (hints->items[i].length == len &&
            memcmp(hints->items[i].sequence, seq, len) == 0) {
            /* Update frequency if higher */
            if (freq > hints->items[i].frequency) {
                hints->items[i].frequency = freq;
            }
            return;
        }
    }

    /* Expand if needed */
    if (hints->count >= hints->capacity) {
        hints->capacity *= 2;
        hints->items = realloc(hints->items, hints->capacity * sizeof(Hint));
    }

    /* Add new hint */
    hints->items[hints->count].sequence = malloc(len + 1);
    memcpy(hints->items[hints->count].sequence, seq, len);
    hints->items[hints->count].sequence[len] = '\0';
    hints->items[hints->count].length = len;
    hints->items[hints->count].frequency = freq;
    hints->count++;
}

/* Build and iterate FS-tree to get hints for one starting nucleotide */
static void fstree_iterate_nucleotide(const char *array, size_t len,
                                       char start_nuc, int depth, size_t cutoff,
                                       HintArray *hints) {
    /* Find all positions of starting nucleotide */
    size_t *start_positions = malloc(len * sizeof(size_t));
    size_t start_count = 0;

    for (size_t i = 0; i < len; i++) {
        if (array[i] == start_nuc) {
            start_positions[start_count++] = i;
        }
    }

    if (start_count <= cutoff) {
        free(start_positions);
        return;
    }

    /* Initialize heap with starting node */
    MinHeap heap;
    heap_init(&heap, 1000);

    HeapNode initial;
    initial.level = 1;
    initial.names = malloc(start_count * sizeof(size_t));
    initial.positions = malloc(start_count * sizeof(size_t));
    memcpy(initial.names, start_positions, start_count * sizeof(size_t));
    memcpy(initial.positions, start_positions, start_count * sizeof(size_t));
    initial.count = start_count;

    heap_push(&heap, initial);
    free(start_positions);

    /* Track found patterns to avoid duplicates from self-repeating */
    char **found_patterns = malloc(MAX_HINTS * sizeof(char*));
    size_t found_count = 0;

    /* BFS traversal */
    size_t current_level = 0;
    size_t *level_buffer = malloc(1000 * sizeof(size_t));  /* Store (start, end, count) */
    size_t buffer_count = 0;

    while (heap.count > 0) {
        HeapNode node = heap_pop(&heap);

        /* Process level change - emit best hint from previous level */
        if (node.level != current_level) {
            if (buffer_count > 0 && current_level > 0) {
                /* Find best in buffer (highest count) */
                size_t best_idx = 0;
                size_t best_count = 0;
                for (size_t i = 0; i < buffer_count; i += 3) {
                    if (level_buffer[i + 2] > best_count) {
                        best_count = level_buffer[i + 2];
                        best_idx = i;
                    }
                }

                size_t start = level_buffer[best_idx];
                size_t seq_len = current_level;

                /* Extract sequence */
                char *found_seq = malloc(seq_len + 1);
                memcpy(found_seq, array + start, seq_len);
                found_seq[seq_len] = '\0';

                /* Check if self-repeating */
                char *minimal = is_self_repeating(found_seq, seq_len);
                if (minimal) {
                    /* Check if minimal unit already found */
                    bool already_found = false;
                    for (size_t i = 0; i < found_count; i++) {
                        if (strcmp(found_patterns[i], minimal) == 0) {
                            already_found = true;
                            break;
                        }
                    }

                    if (!already_found) {
                        size_t min_len = strlen(minimal);
                        size_t min_count = 0;
                        for (size_t i = 0; i <= len - min_len; i++) {
                            if (memcmp(array + i, minimal, min_len) == 0) {
                                min_count++;
                            }
                        }
                        hints_add(hints, minimal, min_len, min_count);

                        if (found_count < MAX_HINTS) {
                            found_patterns[found_count++] = strdup(minimal);
                        }
                    }
                    free(minimal);
                } else {
                    hints_add(hints, found_seq, seq_len, best_count);
                    if (found_count < MAX_HINTS) {
                        found_patterns[found_count++] = strdup(found_seq);
                    }
                }

                free(found_seq);
            }

            buffer_count = 0;
            current_level = node.level;

            if ((int)current_level > depth) {
                free(node.names);
                free(node.positions);
                break;
            }
        }

        /* Add to buffer */
        if (buffer_count + 3 <= 1000) {
            level_buffer[buffer_count++] = node.names[0];
            level_buffer[buffer_count++] = node.positions[0];
            level_buffer[buffer_count++] = node.count;
        }

        /* Expand node - group by next nucleotide */
        size_t next_a = 0, next_c = 0, next_g = 0, next_t = 0;
        size_t *names_a = malloc(node.count * sizeof(size_t));
        size_t *names_c = malloc(node.count * sizeof(size_t));
        size_t *names_g = malloc(node.count * sizeof(size_t));
        size_t *names_t = malloc(node.count * sizeof(size_t));
        size_t *pos_a = malloc(node.count * sizeof(size_t));
        size_t *pos_c = malloc(node.count * sizeof(size_t));
        size_t *pos_g = malloc(node.count * sizeof(size_t));
        size_t *pos_t = malloc(node.count * sizeof(size_t));

        for (size_t i = 0; i < node.count; i++) {
            size_t pos = node.positions[i];
            if (pos + 1 >= len) continue;

            char nuc = array[pos + 1];
            switch (nuc) {
                case 'A':
                    names_a[next_a] = node.names[i];
                    pos_a[next_a++] = pos + 1;
                    break;
                case 'C':
                    names_c[next_c] = node.names[i];
                    pos_c[next_c++] = pos + 1;
                    break;
                case 'G':
                    names_g[next_g] = node.names[i];
                    pos_g[next_g++] = pos + 1;
                    break;
                case 'T':
                    names_t[next_t] = node.names[i];
                    pos_t[next_t++] = pos + 1;
                    break;
            }
        }

        /* Push children that pass cutoff */
        if (next_a > cutoff) {
            HeapNode child;
            child.level = node.level + 1;
            child.names = realloc(names_a, next_a * sizeof(size_t));
            child.positions = realloc(pos_a, next_a * sizeof(size_t));
            child.count = next_a;
            heap_push(&heap, child);
        } else {
            free(names_a);
            free(pos_a);
        }

        if (next_c > cutoff) {
            HeapNode child;
            child.level = node.level + 1;
            child.names = realloc(names_c, next_c * sizeof(size_t));
            child.positions = realloc(pos_c, next_c * sizeof(size_t));
            child.count = next_c;
            heap_push(&heap, child);
        } else {
            free(names_c);
            free(pos_c);
        }

        if (next_g > cutoff) {
            HeapNode child;
            child.level = node.level + 1;
            child.names = realloc(names_g, next_g * sizeof(size_t));
            child.positions = realloc(pos_g, next_g * sizeof(size_t));
            child.count = next_g;
            heap_push(&heap, child);
        } else {
            free(names_g);
            free(pos_g);
        }

        if (next_t > cutoff) {
            HeapNode child;
            child.level = node.level + 1;
            child.names = realloc(names_t, next_t * sizeof(size_t));
            child.positions = realloc(pos_t, next_t * sizeof(size_t));
            child.count = next_t;
            heap_push(&heap, child);
        } else {
            free(names_t);
            free(pos_t);
        }

        free(node.names);
        free(node.positions);
    }

    /* Cleanup */
    heap_free(&heap);
    for (size_t i = 0; i < found_count; i++) {
        free(found_patterns[i]);
    }
    free(found_patterns);
    free(level_buffer);
}

/* Build and iterate FS-tree to get hints from all nucleotides */
HintArray* fstree_get_hints(const char *array, size_t len, int depth, size_t cutoff) {
    HintArray *hints = malloc(sizeof(HintArray));
    hints->capacity = 100;
    hints->count = 0;
    hints->items = malloc(hints->capacity * sizeof(Hint));

    /* Iterate over all starting nucleotides */
    const char nucleotides[] = "ACGT";
    for (int i = 0; i < 4; i++) {
        fstree_iterate_nucleotide(array, len, nucleotides[i], depth, cutoff, hints);
    }

    return hints;
}

/* Free hint array */
void hints_free(HintArray *hints) {
    if (!hints) return;

    for (size_t i = 0; i < hints->count; i++) {
        free(hints->items[i].sequence);
    }
    free(hints->items);
    free(hints);
}
