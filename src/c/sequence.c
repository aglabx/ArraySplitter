/*
 * Sequence utility functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "arraysplitter.h"

/* Complement table */
static const char COMPLEMENT[256] = {
    ['A'] = 'T', ['T'] = 'A', ['C'] = 'G', ['G'] = 'C',
    ['a'] = 't', ['t'] = 'a', ['c'] = 'g', ['g'] = 'c',
    ['N'] = 'N', ['n'] = 'n',
    [' '] = ' ', ['~'] = '~', ['['] = ']', [']'] = '['
};

/* Get reverse complement of sequence */
char* seq_revcomp(const char *seq, size_t len) {
    char *result = malloc(len + 1);
    if (!result) return NULL;

    for (size_t i = 0; i < len; i++) {
        char c = COMPLEMENT[(unsigned char)seq[len - 1 - i]];
        result[i] = c ? c : seq[len - 1 - i];
    }
    result[len] = '\0';

    return result;
}

/* Check if sequence is in canonical orientation (A > T, C > G) */
bool seq_is_canonical(const char *seq, size_t len) {
    size_t a_count = 0, t_count = 0, c_count = 0, g_count = 0;

    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
            case 'A': case 'a': a_count++; break;
            case 'T': case 't': t_count++; break;
            case 'C': case 'c': c_count++; break;
            case 'G': case 'g': g_count++; break;
        }
    }

    /* Primary criterion: A > T */
    if (a_count != t_count) {
        return a_count > t_count;
    }

    /* Secondary criterion: C > G (when A == T) */
    return c_count >= g_count;
}

/* Convert to uppercase in place */
void seq_to_upper(char *seq, size_t len) {
    for (size_t i = 0; i < len; i++) {
        seq[i] = toupper((unsigned char)seq[i]);
    }
}

/* Count nucleotide in sequence */
size_t seq_count_char(const char *seq, size_t len, char c) {
    size_t count = 0;
    char upper_c = toupper((unsigned char)c);

    for (size_t i = 0; i < len; i++) {
        if (toupper((unsigned char)seq[i]) == upper_c) {
            count++;
        }
    }

    return count;
}

/* Get most frequent nucleotide */
char seq_most_frequent_nucleotide(const char *seq, size_t len) {
    size_t counts[4] = {0};

    for (size_t i = 0; i < len; i++) {
        switch (seq[i]) {
            case 'A': case 'a': counts[0]++; break;
            case 'C': case 'c': counts[1]++; break;
            case 'G': case 'g': counts[2]++; break;
            case 'T': case 't': counts[3]++; break;
        }
    }

    size_t max_idx = 0;
    for (int i = 1; i < 4; i++) {
        if (counts[i] > counts[max_idx]) {
            max_idx = i;
        }
    }

    const char nucleotides[] = "ACGT";
    return nucleotides[max_idx];
}

/* Find all positions of a substring using explicit length */
size_t* seq_find_all(const char *seq, size_t seq_len,
                     const char *pattern, size_t pat_len,
                     size_t *out_count) {
    if (!seq || !pattern || pat_len == 0 || seq_len < pat_len) {
        *out_count = 0;
        return NULL;
    }

    /* Initial allocation */
    size_t capacity = 1000;
    size_t *positions = malloc(capacity * sizeof(size_t));
    if (!positions) {
        *out_count = 0;
        return NULL;
    }

    size_t count = 0;

    /* Search using memcmp with explicit pattern length */
    for (size_t i = 0; i <= seq_len - pat_len; i++) {
        if (memcmp(seq + i, pattern, pat_len) == 0) {
            /* Expand if needed */
            if (count >= capacity) {
                capacity *= 2;
                positions = realloc(positions, capacity * sizeof(size_t));
                if (!positions) {
                    *out_count = 0;
                    return NULL;
                }
            }

            positions[count++] = i;
        }
    }

    *out_count = count;
    return positions;
}

/* Split sequence by delimiter, returning array of parts */
char** seq_split(const char *seq, size_t seq_len,
                 const char *delim, size_t delim_len,
                 size_t *out_count, size_t **out_lengths) {
    if (!seq || !delim || delim_len == 0) {
        *out_count = 0;
        return NULL;
    }

    /* Find all occurrences of delimiter */
    size_t pos_count;
    size_t *positions = seq_find_all(seq, seq_len, delim, delim_len, &pos_count);

    /* Calculate number of parts */
    size_t num_parts = pos_count + 1;

    /* Allocate arrays */
    char **parts = malloc(num_parts * sizeof(char*));
    size_t *lengths = malloc(num_parts * sizeof(size_t));
    if (!parts || !lengths) {
        free(positions);
        free(parts);
        free(lengths);
        *out_count = 0;
        return NULL;
    }

    /* Filter out overlapping positions */
    size_t *filtered_pos = malloc(pos_count * sizeof(size_t));
    size_t filtered_count = 0;
    size_t last_end = 0;

    for (size_t i = 0; i < pos_count; i++) {
        if (positions[i] >= last_end) {
            filtered_pos[filtered_count++] = positions[i];
            last_end = positions[i] + delim_len;
        }
    }
    free(positions);
    positions = filtered_pos;
    pos_count = filtered_count;

    /* Recalculate number of parts */
    num_parts = pos_count + 1;

    /* Reallocate arrays for filtered size */
    free(parts);
    free(lengths);
    parts = malloc(num_parts * sizeof(char*));
    lengths = malloc(num_parts * sizeof(size_t));
    if (!parts || !lengths) {
        free(positions);
        free(parts);
        free(lengths);
        *out_count = 0;
        return NULL;
    }

    /* Extract parts */
    size_t prev_end = 0;
    for (size_t i = 0; i < pos_count; i++) {
        size_t start = prev_end;
        size_t end = positions[i];
        size_t len = end - start;

        parts[i] = malloc(len + 1);
        if (parts[i]) {
            memcpy(parts[i], seq + start, len);
            parts[i][len] = '\0';
        }
        lengths[i] = len;

        prev_end = positions[i] + delim_len;
    }

    /* Last part */
    size_t last_len = seq_len - prev_end;
    parts[pos_count] = malloc(last_len + 1);
    if (parts[pos_count]) {
        memcpy(parts[pos_count], seq + prev_end, last_len);
        parts[pos_count][last_len] = '\0';
    }
    lengths[pos_count] = last_len;

    free(positions);
    *out_count = num_parts;
    *out_lengths = lengths;
    return parts;
}

/* Free split result */
void seq_split_free(char **parts, size_t count, size_t *lengths) {
    if (parts) {
        for (size_t i = 0; i < count; i++) {
            free(parts[i]);
        }
        free(parts);
    }
    free(lengths);
}
