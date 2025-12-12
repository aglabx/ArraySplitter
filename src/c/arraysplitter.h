/*
 * ArraySplitter - C implementation with multithreading
 *
 * De novo decomposition of satellite DNA arrays into monomers.
 *
 * Author: Aleksey Komissarov
 */

#ifndef ARRAYSPLITTER_H
#define ARRAYSPLITTER_H

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <pthread.h>

/* Maximum sequence length and other limits */
#define MAX_SEQ_LENGTH    100000000  /* 100 Mb max array */
#define MAX_HEADER_LENGTH 4096
#define MAX_CUT_LENGTH    100
#define MAX_HINTS         10000
#define MAX_MONOMERS      1000000

/* Nucleotide to index conversion */
#define NUC_A 0
#define NUC_C 1
#define NUC_G 2
#define NUC_T 3
#define NUC_N 4

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/* Single sequence from FASTA file */
typedef struct {
    char *header;
    char *sequence;
    size_t length;
} FastaEntry;

/* Array of FASTA entries */
typedef struct {
    FastaEntry *entries;
    size_t count;
    size_t capacity;
} FastaArray;

/* Hint from FS-tree (candidate cut sequence) */
typedef struct {
    char *sequence;      /* The cut sequence */
    size_t length;       /* Length of cut sequence */
    size_t frequency;    /* How many times it appears */
} Hint;

/* Array of hints */
typedef struct {
    Hint *items;
    size_t count;
    size_t capacity;
} HintArray;

/* Candidate cut with scoring */
typedef struct {
    char *cut_seq;
    size_t cut_len;
    size_t mode_period;
    double base_score;
    double adjusted_score;
    size_t num_segments;
    double fragmentation;
} CutCandidate;

/* Decomposition result for single array */
typedef struct {
    char **monomers;         /* Array of monomer strings */
    size_t *monomer_lengths; /* Length of each monomer */
    size_t monomer_count;
    char *cut_seq;
    size_t period;
    double score;
    bool was_reversed;       /* True if array was reverse-complemented */
} DecompResult;

/* Thread work item */
typedef struct {
    FastaEntry *entry;
    DecompResult *result;
    int depth;
    bool verbose;
} ThreadWorkItem;

/* Thread pool work queue */
typedef struct WorkQueue {
    ThreadWorkItem *items;
    size_t count;
    size_t capacity;
    size_t next_item;      /* Next item to process */
    pthread_mutex_t mutex;
} WorkQueue;

/* ============================================================================
 * FASTA I/O Functions
 * ============================================================================ */

/* Read all entries from FASTA file */
FastaArray* fasta_read_file(const char *filename);

/* Free FASTA array */
void fasta_free(FastaArray *fa);

/* ============================================================================
 * Sequence Utilities
 * ============================================================================ */

/* Get reverse complement of sequence */
char* seq_revcomp(const char *seq, size_t len);

/* Check if sequence is in canonical orientation (A > T, C > G) */
bool seq_is_canonical(const char *seq, size_t len);

/* Convert to uppercase */
void seq_to_upper(char *seq, size_t len);

/* Count nucleotide in sequence */
size_t seq_count_char(const char *seq, size_t len, char c);

/* Get most frequent nucleotide */
char seq_most_frequent_nucleotide(const char *seq, size_t len);

/* ============================================================================
 * FS-Tree Functions
 * ============================================================================ */

/* Position list for FS-tree nodes */
typedef struct {
    size_t *positions;
    size_t count;
    size_t capacity;
} PositionList;

/* FS-tree node */
typedef struct FSNode {
    PositionList names;      /* Starting positions (names) */
    PositionList positions;  /* Current positions */
    size_t level;            /* Current depth in tree */
} FSNode;

/* Build and iterate FS-tree to get hints */
HintArray* fstree_get_hints(const char *array, size_t len, int depth, size_t cutoff);

/* Free hint array */
void hints_free(HintArray *hints);

/* ============================================================================
 * Decomposition Functions
 * ============================================================================ */

/* Find best cut sequence from hints */
CutCandidate find_best_cut(const char *array, size_t len, HintArray *hints);

/* Decompose array using cut sequence */
DecompResult* decompose_array(const char *array, size_t len,
                              const char *cut_seq, size_t cut_len);

/* Full decomposition pipeline for single array */
DecompResult* decompose_full(const char *header, const char *array, size_t len,
                             int depth, bool verbose);

/* Free decomposition result */
void decomp_result_free(DecompResult *result);

/* ============================================================================
 * Multithreading
 * ============================================================================ */

/* Process arrays in parallel */
DecompResult** process_parallel(FastaArray *arrays, int num_threads,
                                 int depth, bool verbose);

/* ============================================================================
 * Output Functions
 * ============================================================================ */

/* Write results to output files */
int write_results(const char *output_prefix, FastaArray *arrays,
                  DecompResult **results, size_t count);

#endif /* ARRAYSPLITTER_H */
