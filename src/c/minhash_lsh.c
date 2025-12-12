/*
 * minhash_lsh.c - Locality Sensitive Hashing for MinHash
 *
 * LSH with b bands of r rows:
 * - Two signatures match in a band if all r hash values match
 * - Probability of being candidates = 1 - (1 - s^r)^b
 * where s is the Jaccard similarity
 *
 * With b=16 bands, r=8 rows:
 * - s=0.8 -> P(candidate) ≈ 0.996 (almost certain)
 * - s=0.5 -> P(candidate) ≈ 0.470
 * - s=0.2 -> P(candidate) ≈ 0.006 (very unlikely)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minhash.h"

/* Hash a band (r consecutive hash values) to a bucket index */
static uint32_t hash_band(const uint64_t *hashes, int start, int r) {
    uint64_t h = 0;
    for (int i = 0; i < r; i++) {
        h ^= hashes[start + i];
        h = (h * 0x9e3779b97f4a7c15ULL) >> 32;
    }
    return h % LSH_BUCKETS;
}

/* Initialize LSH index */
void lsh_init(LSHIndex *lsh) {
    for (int b = 0; b < LSH_BANDS; b++) {
        for (int i = 0; i < LSH_BUCKETS; i++) {
            lsh->bands[b].buckets[i] = NULL;
        }
    }
}

/* Free LSH index */
void lsh_free(LSHIndex *lsh) {
    for (int b = 0; b < LSH_BANDS; b++) {
        for (int i = 0; i < LSH_BUCKETS; i++) {
            LSHEntry *e = lsh->bands[b].buckets[i];
            while (e) {
                LSHEntry *next = e->next;
                free(e);
                e = next;
            }
        }
    }
}

/* Clear LSH index (free entries, keep structure) */
void lsh_clear(LSHIndex *lsh) {
    for (int b = 0; b < LSH_BANDS; b++) {
        for (int i = 0; i < LSH_BUCKETS; i++) {
            LSHEntry *e = lsh->bands[b].buckets[i];
            while (e) {
                LSHEntry *next = e->next;
                free(e);
                e = next;
            }
            lsh->bands[b].buckets[i] = NULL;
        }
    }
}

/* Insert a monomer into LSH index */
void lsh_insert(LSHIndex *lsh, uint32_t monomer_id, const MinHashSig *sig) {
    for (int b = 0; b < LSH_BANDS; b++) {
        uint32_t bucket = hash_band(sig->hashes, b * LSH_ROWS, LSH_ROWS);

        LSHEntry *entry = malloc(sizeof(LSHEntry));
        entry->monomer_id = monomer_id;
        entry->next = lsh->bands[b].buckets[bucket];
        lsh->bands[b].buckets[bucket] = entry;
    }
}

/*
 * Query for candidate matches
 *
 * Returns array of unique monomer IDs that share at least one band.
 * Caller must free the returned array.
 */
uint32_t* lsh_query_candidates(LSHIndex *lsh, const MinHashSig *sig, uint32_t *count) {
    /* Use a simple seen array to deduplicate
     * For large datasets, could use a hash set */
    static uint32_t *seen = NULL;
    static uint32_t seen_capacity = 0;
    static uint32_t *candidates = NULL;
    static uint32_t candidates_capacity = 0;

    /* Lazy initialization */
    if (seen_capacity < LSH_BUCKETS * 10) {
        seen_capacity = LSH_BUCKETS * 10;
        seen = realloc(seen, seen_capacity * sizeof(uint32_t));
        candidates_capacity = 10000;
        candidates = realloc(candidates, candidates_capacity * sizeof(uint32_t));
    }

    uint32_t seen_count = 0;
    *count = 0;

    /* Check each band */
    for (int b = 0; b < LSH_BANDS; b++) {
        uint32_t bucket = hash_band(sig->hashes, b * LSH_ROWS, LSH_ROWS);

        LSHEntry *e = lsh->bands[b].buckets[bucket];
        while (e) {
            uint32_t id = e->monomer_id;

            /* Check if already seen (linear search, ok for typical candidate counts) */
            bool found = false;
            for (uint32_t i = 0; i < seen_count; i++) {
                if (seen[i] == id) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                /* Add to seen and candidates */
                if (seen_count >= seen_capacity) {
                    seen_capacity *= 2;
                    seen = realloc(seen, seen_capacity * sizeof(uint32_t));
                }
                seen[seen_count++] = id;

                if (*count >= candidates_capacity) {
                    candidates_capacity *= 2;
                    candidates = realloc(candidates, candidates_capacity * sizeof(uint32_t));
                }
                candidates[(*count)++] = id;
            }

            e = e->next;
        }
    }

    /* Return a copy of candidates array */
    uint32_t *result = malloc((*count) * sizeof(uint32_t));
    memcpy(result, candidates, (*count) * sizeof(uint32_t));
    return result;
}

/*
 * Build LSH index for a subset of monomers
 */
void lsh_build(MHContext *ctx, uint32_t *ids, uint32_t count) {
    lsh_clear(&ctx->lsh);

    for (uint32_t i = 0; i < count; i++) {
        uint32_t id = ids[i];
        lsh_insert(&ctx->lsh, id, &ctx->monomers[id].signature);
    }
}
