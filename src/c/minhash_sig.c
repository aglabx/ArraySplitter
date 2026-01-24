/*
 * minhash_sig.c - MinHash signature computation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <pthread.h>
#include "minhash.h"

/* Pre-computed random seeds for hash functions */
static uint32_t hash_seeds[MINHASH_SIZE];
static bool seeds_initialized = false;

/* Initialize hash seeds (deterministic for reproducibility) */
static void init_hash_seeds(void) {
    if (seeds_initialized) return;

    /* Use a simple LCG to generate deterministic seeds */
    uint32_t seed = 12345;
    for (int i = 0; i < MINHASH_SIZE; i++) {
        seed = seed * 1103515245 + 12345;
        hash_seeds[i] = seed;
    }
    seeds_initialized = true;
}

/*
 * Check if k-mer has low complexity (should be filtered)
 * Returns true if k-mer should be SKIPPED
 *
 * Criteria:
 * 1. Less than 3 distinct nucleotides
 * 2. Any single nucleotide > 75% of k-mer
 */
static inline bool is_low_complexity_kmer(uint64_t kmer, int k) {
    int counts[4] = {0, 0, 0, 0};  /* A, C, G, T */

    /* Count each nucleotide */
    uint64_t tmp = kmer;
    for (int i = 0; i < k; i++) {
        counts[tmp & 3]++;
        tmp >>= 2;
    }

    /* Check distinct count */
    int distinct = (counts[0] > 0) + (counts[1] > 0) +
                   (counts[2] > 0) + (counts[3] > 0);
    if (distinct < 3) {
        return true;  /* Filter: only 1-2 nucleotide types */
    }

    /* Check for dominant nucleotide (>75%) */
    int threshold = (k * 3) / 4;  /* 75% of k */
    for (int i = 0; i < 4; i++) {
        if (counts[i] > threshold) {
            return true;  /* Filter: one nucleotide dominates */
        }
    }

    return false;  /* Keep this k-mer */
}

/* Global statistics for low-complexity filter */
static uint64_t total_kmers_processed = 0;
static uint64_t total_kmers_filtered = 0;
static pthread_mutex_t stats_mutex = PTHREAD_MUTEX_INITIALIZER;

/*
 * Compute MinHash signature for a monomer
 *
 * For each of 128 hash functions:
 *   - Hash all canonical k-mers in sequence
 *   - Keep minimum hash value
 *   - Skip low-complexity k-mers (AT-rich, homopolymers)
 */
void compute_minhash(MHMonomer *mono) {
    init_hash_seeds();

    /* Initialize all hashes to maximum */
    for (int i = 0; i < MINHASH_SIZE; i++) {
        mono->signature.hashes[i] = UINT64_MAX;
    }

    /* Need at least k nucleotides */
    if (mono->length < MINHASH_K) {
        return;
    }

    /* Sliding window over sequence */
    const char *seq = mono->sequence;
    int num_kmers = mono->length - MINHASH_K + 1;

    uint32_t local_total = 0;
    uint32_t local_filtered = 0;

    for (int pos = 0; pos < num_kmers; pos++) {
        /* Check for N in k-mer window */
        bool has_n = false;
        for (int j = 0; j < MINHASH_K; j++) {
            char c = seq[pos + j];
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                c != 'a' && c != 'c' && c != 'g' && c != 't') {
                has_n = true;
                break;
            }
        }
        if (has_n) continue;

        local_total++;

        /* Encode k-mer */
        uint64_t kmer = encode_kmer(seq + pos, MINHASH_K);

        /* Get canonical form */
        uint64_t canon = canonical_kmer(kmer, MINHASH_K);

        /* Filter low-complexity k-mers */
        if (is_low_complexity_kmer(canon, MINHASH_K)) {
            local_filtered++;
            continue;  /* Skip AT-rich, homopolymer-like k-mers */
        }

        /* Update minimum for each hash function */
        for (int i = 0; i < MINHASH_SIZE; i++) {
            uint64_t h = hash_kmer(canon, hash_seeds[i]);
            if (h < mono->signature.hashes[i]) {
                mono->signature.hashes[i] = h;
            }
        }
    }

    /* Update global statistics (thread-safe) */
    pthread_mutex_lock(&stats_mutex);
    total_kmers_processed += local_total;
    total_kmers_filtered += local_filtered;
    pthread_mutex_unlock(&stats_mutex);
}

/* Get filter statistics */
void get_kmer_filter_stats(uint64_t *processed, uint64_t *filtered) {
    *processed = total_kmers_processed;
    *filtered = total_kmers_filtered;
}

/* Reset filter statistics */
void reset_kmer_filter_stats(void) {
    pthread_mutex_lock(&stats_mutex);
    total_kmers_processed = 0;
    total_kmers_filtered = 0;
    pthread_mutex_unlock(&stats_mutex);
}

/* Thread work for parallel signature computation */
typedef struct {
    MHContext *ctx;
    uint32_t start;
    uint32_t end;
    uint32_t thread_id;
} SigThreadWork;

static void* compute_signatures_thread(void *arg) {
    SigThreadWork *work = (SigThreadWork*)arg;

    for (uint32_t i = work->start; i < work->end; i++) {
        compute_minhash(&work->ctx->monomers[i]);
    }

    return NULL;
}

/*
 * Compute signatures for all monomers with progress bar (parallel)
 */
void compute_all_signatures(MHContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "Computing MinHash signatures for %u monomers (%d threads)...\n",
                ctx->monomer_count, ctx->num_threads);
    }

    int num_threads = ctx->num_threads;
    if (num_threads <= 1) {
        /* Single-threaded fallback */
        uint32_t progress_step = ctx->monomer_count / 100;
        if (progress_step == 0) progress_step = 1;

        for (uint32_t i = 0; i < ctx->monomer_count; i++) {
            compute_minhash(&ctx->monomers[i]);
            if (ctx->verbose && (i % progress_step == 0 || i == ctx->monomer_count - 1)) {
                int pct = (i + 1) * 100 / ctx->monomer_count;
                fprintf(stderr, "\r  Signatures: %3d%% (%u/%u)",
                        pct, i + 1, ctx->monomer_count);
            }
        }
    } else {
        /* Multi-threaded */
        pthread_t *threads = malloc(num_threads * sizeof(pthread_t));
        SigThreadWork *works = malloc(num_threads * sizeof(SigThreadWork));

        uint32_t chunk = ctx->monomer_count / num_threads;
        for (int t = 0; t < num_threads; t++) {
            works[t].ctx = ctx;
            works[t].start = t * chunk;
            works[t].end = (t == num_threads - 1) ? ctx->monomer_count : (t + 1) * chunk;
            works[t].thread_id = t;
            pthread_create(&threads[t], NULL, compute_signatures_thread, &works[t]);
        }

        for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }

        free(threads);
        free(works);
    }

    if (ctx->verbose) {
        fprintf(stderr, "\r  Signatures: 100%% (%u/%u)\n", ctx->monomer_count, ctx->monomer_count);

        /* Print low-complexity filter statistics */
        uint64_t processed, filtered;
        get_kmer_filter_stats(&processed, &filtered);
        if (processed > 0) {
            fprintf(stderr, "  Low-complexity k-mers filtered: %llu / %llu (%.1f%%)\n",
                    (unsigned long long)filtered, (unsigned long long)processed,
                    100.0 * filtered / processed);
        }
    }
}
