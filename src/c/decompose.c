/*
 * Decomposition algorithm
 *
 * Finds optimal cut sequence and splits array into monomers.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "arraysplitter.h"

/* External sequence functions */
extern size_t* seq_find_all(const char *seq, size_t seq_len,
                            const char *pattern, size_t pat_len,
                            size_t *out_count);
extern char** seq_split(const char *seq, size_t seq_len,
                        const char *delim, size_t delim_len,
                        size_t *out_count, size_t **out_lengths);
extern void seq_split_free(char **parts, size_t count, size_t *lengths);

/* Period counter for finding mode */
typedef struct {
    size_t period;
    size_t count;
} PeriodCount;

/* Find mode (most common value) in array */
static size_t find_mode(size_t *values, size_t count, size_t *mode_count) {
    if (count == 0) {
        *mode_count = 0;
        return 0;
    }

    /* Count frequencies (simple approach for reasonable period ranges) */
    size_t max_period = 0;
    for (size_t i = 0; i < count; i++) {
        if (values[i] > max_period) max_period = values[i];
    }

    if (max_period == 0 || max_period > 100000) {
        *mode_count = count;
        return values[0];
    }

    /* Count occurrences */
    size_t *freq = calloc(max_period + 1, sizeof(size_t));
    for (size_t i = 0; i < count; i++) {
        if (values[i] <= max_period) {
            freq[values[i]]++;
        }
    }

    /* Find mode */
    size_t mode = 0;
    size_t max_freq = 0;
    for (size_t i = 1; i <= max_period; i++) {
        if (freq[i] > max_freq) {
            max_freq = freq[i];
            mode = i;
        }
    }

    free(freq);
    *mode_count = max_freq;
    return mode;
}

/* Evaluate a cut candidate */
static CutCandidate evaluate_cut(const char *array, size_t array_len,
                                  const char *cut_seq, size_t cut_len,
                                  size_t frequency) {
    (void)frequency;  /* Currently unused, kept for future scoring */
    CutCandidate cand = {0};

    cand.cut_seq = strdup(cut_seq);
    cand.cut_len = cut_len;

    /* Split array by cut sequence */
    size_t num_parts;
    size_t *part_lengths;
    char **parts = seq_split(array, array_len, cut_seq, cut_len, &num_parts, &part_lengths);

    if (!parts || num_parts == 0) {
        cand.base_score = 0;
        return cand;
    }

    /* Calculate periods (length + cut_len for each part) */
    size_t *periods = malloc(num_parts * sizeof(size_t));
    size_t *non_empty_periods = malloc(num_parts * sizeof(size_t));
    size_t period_count = 0;
    size_t non_empty_count = 0;
    size_t empty_count = 0;

    for (size_t i = 0; i < num_parts; i++) {
        /* All except possibly empty last */
        if (i < num_parts - 1 || part_lengths[i] > 0) {
            size_t period = part_lengths[i] + cut_len;
            periods[period_count++] = period;

            if (part_lengths[i] > 0) {
                non_empty_periods[non_empty_count++] = period;
            } else {
                empty_count++;
            }
        }
    }

    if (period_count == 0) {
        seq_split_free(parts, num_parts, part_lengths);
        free(periods);
        free(non_empty_periods);
        cand.base_score = 0;
        return cand;
    }

    /* Calculate empty ratio */
    double empty_ratio = (double)empty_count / period_count;

    size_t mode_period, mode_count;
    size_t total_segments;

    if (empty_ratio >= 0.8) {
        /* Perfect or near-perfect repeat */
        mode_period = cut_len;
        mode_count = period_count;
        total_segments = period_count;
    } else if (non_empty_count > 0) {
        /* Use non-empty periods */
        mode_period = find_mode(non_empty_periods, non_empty_count, &mode_count);
        total_segments = non_empty_count;
    } else {
        /* Fallback */
        mode_period = find_mode(periods, period_count, &mode_count);
        total_segments = period_count;
    }

    /* Calculate scores */
    double base_score = (double)mode_count / total_segments;

    /* Fragmentation penalty */
    double short_threshold = mode_period * 0.5;
    size_t short_fragments = 0;
    for (size_t i = 0; i < period_count; i++) {
        if (periods[i] < short_threshold) {
            short_fragments++;
        }
    }
    double fragmentation = (double)short_fragments / total_segments;

    cand.mode_period = mode_period;
    cand.base_score = base_score;
    cand.fragmentation = fragmentation;
    cand.adjusted_score = base_score * (1.0 - fragmentation * 0.5);
    cand.num_segments = total_segments;

    seq_split_free(parts, num_parts, part_lengths);
    free(periods);
    free(non_empty_periods);

    return cand;
}

/* Find best cut sequence from hints */
CutCandidate find_best_cut(const char *array, size_t len, HintArray *hints) {
    CutCandidate best = {0};
    best.base_score = -1;

    if (!hints || hints->count == 0) {
        /* No hints, return whole array */
        best.cut_seq = strdup("");
        best.mode_period = len;
        return best;
    }

    double score_threshold = 0.05;

    /* Evaluate all hints */
    CutCandidate *candidates = malloc(hints->count * sizeof(CutCandidate));
    size_t cand_count = 0;

    for (size_t i = 0; i < hints->count; i++) {
        Hint *h = &hints->items[i];

        /* Skip if cut not in array - use explicit length search */
        size_t pos_count;
        size_t *positions = seq_find_all(array, len, h->sequence, h->length, &pos_count);
        if (pos_count == 0) {
            free(positions);
            continue;
        }
        free(positions);

        CutCandidate cand = evaluate_cut(array, len, h->sequence, h->length, h->frequency);
        if (cand.base_score >= 0) {
            candidates[cand_count++] = cand;
        }
    }

    if (cand_count == 0) {
        free(candidates);
        best.cut_seq = strdup("");
        best.mode_period = len;
        return best;
    }

    /* Find best base score */
    double best_base_score = 0;
    for (size_t i = 0; i < cand_count; i++) {
        if (candidates[i].base_score > best_base_score) {
            best_base_score = candidates[i].base_score;
        }
    }

    /* Get candidates within threshold */
    CutCandidate *similar = malloc(cand_count * sizeof(CutCandidate));
    size_t similar_count = 0;

    for (size_t i = 0; i < cand_count; i++) {
        if (candidates[i].base_score >= best_base_score - score_threshold) {
            similar[similar_count++] = candidates[i];
        } else {
            free(candidates[i].cut_seq);
        }
    }

    /* Sort by adjusted score (descending), then by period (ascending) */
    for (size_t i = 0; i < similar_count - 1; i++) {
        for (size_t j = i + 1; j < similar_count; j++) {
            bool swap = false;
            if (similar[j].adjusted_score > similar[i].adjusted_score) {
                swap = true;
            } else if (similar[j].adjusted_score == similar[i].adjusted_score &&
                       similar[j].mode_period < similar[i].mode_period) {
                swap = true;
            }

            if (swap) {
                CutCandidate tmp = similar[i];
                similar[i] = similar[j];
                similar[j] = tmp;
            }
        }
    }

    /* Take best */
    best = similar[0];

    /* Free others */
    for (size_t i = 1; i < similar_count; i++) {
        free(similar[i].cut_seq);
    }

    free(candidates);
    free(similar);

    return best;
}

/* Decompose array using cut sequence */
DecompResult* decompose_array(const char *array, size_t len,
                              const char *cut_seq, size_t cut_len) {
    DecompResult *result = calloc(1, sizeof(DecompResult));
    if (!result) return NULL;

    /* Copy only cut_len characters */
    result->cut_seq = malloc(cut_len + 1);
    memcpy(result->cut_seq, cut_seq, cut_len);
    result->cut_seq[cut_len] = '\0';
    result->period = cut_len;

    /* Check if cut exists using explicit length */
    size_t pos_count;
    size_t *positions = seq_find_all(array, len, cut_seq, cut_len, &pos_count);
    if (cut_len == 0 || pos_count == 0) {
        /* No cut, return whole array */
        free(positions);
        result->monomers = malloc(sizeof(char*));
        result->monomer_lengths = malloc(sizeof(size_t));
        result->monomers[0] = strdup(array);
        result->monomer_lengths[0] = len;
        result->monomer_count = 1;
        return result;
    }
    free(positions);

    /* Split by cut sequence */
    size_t num_parts;
    size_t *part_lengths;
    char **parts = seq_split(array, len, cut_seq, cut_len, &num_parts, &part_lengths);

    if (!parts) {
        result->monomers = malloc(sizeof(char*));
        result->monomer_lengths = malloc(sizeof(size_t));
        result->monomers[0] = strdup(array);
        result->monomer_lengths[0] = len;
        result->monomer_count = 1;
        return result;
    }

    /* Build monomers: cut + following part (except first which is flank) */
    result->monomers = malloc(num_parts * sizeof(char*));
    result->monomer_lengths = malloc(num_parts * sizeof(size_t));
    result->monomer_count = 0;

    /* First part (flank) - if not empty */
    if (part_lengths[0] > 0) {
        result->monomers[result->monomer_count] = strdup(parts[0]);
        result->monomer_lengths[result->monomer_count] = part_lengths[0];
        result->monomer_count++;
    }

    /* Other parts: prepend cut sequence */
    for (size_t i = 1; i < num_parts; i++) {
        size_t mono_len = cut_len + part_lengths[i];
        char *monomer = malloc(mono_len + 1);

        memcpy(monomer, cut_seq, cut_len);
        memcpy(monomer + cut_len, parts[i], part_lengths[i]);
        monomer[mono_len] = '\0';

        result->monomers[result->monomer_count] = monomer;
        result->monomer_lengths[result->monomer_count] = mono_len;
        result->monomer_count++;
    }

    seq_split_free(parts, num_parts, part_lengths);

    return result;
}

/* Calculate cutoff based on array size */
static size_t calculate_cutoff(size_t array_len) {
    if (array_len > 1000000) return 1000;
    if (array_len > 100000) return 250;
    if (array_len > 10000) return 10;
    return 3;
}

/* Full decomposition pipeline for single array */
DecompResult* decompose_full(const char *header, const char *array, size_t len,
                             int depth, bool verbose) {
    /* Check canonical orientation */
    bool is_canonical = seq_is_canonical(array, len);

    char *work_array;
    bool was_reversed = false;

    if (!is_canonical) {
        work_array = seq_revcomp(array, len);
        was_reversed = true;
    } else {
        work_array = strdup(array);
    }

    /* Calculate cutoff */
    size_t cutoff = calculate_cutoff(len);

    if (verbose) {
        printf("Processing: %s (len=%zu, cutoff=%zu)\n", header, len, cutoff);
    }

    /* Get hints from FS-tree */
    HintArray *hints = fstree_get_hints(work_array, len, depth, cutoff);

    if (verbose) {
        printf("  Found %zu hints\n", hints->count);
    }

    /* Find best cut */
    CutCandidate best = find_best_cut(work_array, len, hints);

    if (verbose) {
        printf("  Best cut: '%s' (score=%.3f, period=%zu)\n",
               best.cut_seq, best.base_score, best.mode_period);
    }

    /* Decompose */
    DecompResult *result = decompose_array(work_array, len,
                                            best.cut_seq, strlen(best.cut_seq));

    result->was_reversed = was_reversed;
    result->score = best.base_score;
    result->period = best.mode_period;

    /* Update cut_seq in result */
    free(result->cut_seq);
    result->cut_seq = best.cut_seq;

    /* Verify reconstruction */
    size_t total_len = 0;
    for (size_t i = 0; i < result->monomer_count; i++) {
        total_len += result->monomer_lengths[i];
    }

    if (total_len != len) {
        fprintf(stderr, "ERROR: Reconstruction failed! %zu != %zu\n", total_len, len);
    }

    /* Cleanup */
    hints_free(hints);
    free(work_array);

    return result;
}

/* Free decomposition result */
void decomp_result_free(DecompResult *result) {
    if (!result) return;

    for (size_t i = 0; i < result->monomer_count; i++) {
        free(result->monomers[i]);
    }
    free(result->monomers);
    free(result->monomer_lengths);
    free(result->cut_seq);
    free(result);
}
