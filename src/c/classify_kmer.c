/*
 * classify_kmer.c - 8-mer bitvector building and index construction
 *
 * Each monomer is represented as a 65536-bit vector (bag of k-mers).
 */

#include "classify.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Decode 16-bit integer to 8-mer string */
void decode_8mer(uint16_t val, char *out) {
    static const char bases[] = "ACGT";
    for (int i = 7; i >= 0; i--) {
        out[i] = bases[val & 0x3];
        val >>= 2;
    }
    out[8] = '\0';
}

/* Build bitvector for a monomer (all 8-mers, sliding window step=1) */
void build_bitvector(Monomer *mono) {
    bitvec_clear(&mono->bitvec);

    if (mono->length < KMER_SIZE) {
        return;
    }

    /* Slide through sequence with step=1 */
    for (uint32_t pos = 0; pos + KMER_SIZE <= mono->length; pos++) {
        /* Check for N in this 8-mer */
        bool has_n = false;
        for (int i = 0; i < KMER_SIZE; i++) {
            char c = mono->sequence[pos + i];
            if (c == 'N' || c == 'n') {
                has_n = true;
                break;
            }
        }

        if (!has_n) {
            uint16_t kmer = encode_8mer(mono->sequence + pos);
            bitvec_set(&mono->bitvec, kmer);
        }
    }

    mono->bitvec.popcount = bitvec_popcount(&mono->bitvec);
}

/* Add hit to bucket (caller ensures no duplicates via bitvec_test) */
static void bucket_add(KmerBucket *bucket, uint32_t monomer_id) {
    if (bucket->count >= bucket->capacity) {
        uint32_t new_cap = bucket->capacity == 0 ? 16 : bucket->capacity * 2;
        bucket->hits = realloc(bucket->hits, new_cap * sizeof(KmerHit));
        bucket->capacity = new_cap;
    }

    bucket->hits[bucket->count].monomer_id = monomer_id;
    bucket->count++;
}

/* Build bitvector AND add to index during sequence scan */
static void build_bitvector_and_index(ClassifyContext *ctx, uint32_t mono_idx) {
    Monomer *mono = &ctx->monomers[mono_idx];
    bitvec_clear(&mono->bitvec);

    if (mono->length < KMER_SIZE) {
        return;
    }

    /* Slide through sequence with step=1 */
    for (uint32_t pos = 0; pos + KMER_SIZE <= mono->length; pos++) {
        /* Check for N in this 8-mer */
        bool has_n = false;
        for (int i = 0; i < KMER_SIZE; i++) {
            char c = mono->sequence[pos + i];
            if (c == 'N' || c == 'n') {
                has_n = true;
                break;
            }
        }

        if (!has_n) {
            uint16_t kmer = encode_8mer(mono->sequence + pos);

            /* Only add to index if this is a new k-mer for this monomer */
            if (!bitvec_test(&mono->bitvec, kmer)) {
                bitvec_set(&mono->bitvec, kmer);
                bucket_add(&ctx->buckets[kmer], mono_idx);
            }
        }
    }

    mono->bitvec.popcount = bitvec_popcount(&mono->bitvec);
}

/* Build inverted index from all monomers */
void build_kmer_index(ClassifyContext *ctx) {
    /* Initialize buckets */
    for (int i = 0; i < KMER_BUCKETS; i++) {
        ctx->buckets[i].hits = NULL;
        ctx->buckets[i].count = 0;
        ctx->buckets[i].capacity = 0;
    }

    /* Count monomers to process for progress */
    uint32_t total = 0;
    for (uint32_t m = 0; m < ctx->monomer_count; m++) {
        if (!ctx->skip_flanks || ctx->monomers[m].type == TYPE_MONOMER) {
            total++;
        }
    }

    uint32_t progress_step = total / 20;
    if (progress_step == 0) progress_step = 1;
    uint32_t processed = 0;

    /* Process each monomer - build bitvector and index simultaneously */
    for (uint32_t m = 0; m < ctx->monomer_count; m++) {
        Monomer *mono = &ctx->monomers[m];

        /* Skip flanks if requested */
        if (ctx->skip_flanks && mono->type != TYPE_MONOMER) {
            continue;
        }

        processed++;

        /* Progress bar */
        if (ctx->verbose && processed % progress_step == 0) {
            uint32_t pct = (processed * 100) / total;
            fprintf(stderr, "\r  Building index: %u%% (%u/%u monomers)", pct, processed, total);
        }

        if (mono->length >= KMER_SIZE) {
            build_bitvector_and_index(ctx, m);
        }
    }

    if (ctx->verbose) {
        /* Count total indexed k-mers */
        uint64_t total_hits = 0;
        uint32_t non_empty_buckets = 0;
        for (int i = 0; i < KMER_BUCKETS; i++) {
            if (ctx->buckets[i].count > 0) {
                total_hits += ctx->buckets[i].count;
                non_empty_buckets++;
            }
        }
        fprintf(stderr, "\r  Index built: %llu monomer-kmer pairs in %u non-empty buckets\n",
                (unsigned long long)total_hits, non_empty_buckets);
    }
}
