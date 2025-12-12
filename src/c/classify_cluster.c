/*
 * classify_cluster.c - Component and subcomponent discovery
 *
 * Algorithm:
 * 1. Find connected components where all members share >= MIN_SHARED_KMERS
 * 2. Within each component, find subcomponents with Jaccard >= threshold
 */

#include "classify.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Union-Find data structure */
typedef struct {
    uint32_t *parent;
    uint32_t *rank;
    uint32_t size;
} UnionFind;

static UnionFind* uf_create(uint32_t size) {
    UnionFind *uf = malloc(sizeof(UnionFind));
    uf->parent = malloc(size * sizeof(uint32_t));
    uf->rank = calloc(size, sizeof(uint32_t));
    uf->size = size;

    for (uint32_t i = 0; i < size; i++) {
        uf->parent[i] = i;
    }
    return uf;
}

static void uf_free(UnionFind *uf) {
    free(uf->parent);
    free(uf->rank);
    free(uf);
}

static uint32_t uf_find(UnionFind *uf, uint32_t x) {
    if (uf->parent[x] != x) {
        uf->parent[x] = uf_find(uf, uf->parent[x]);
    }
    return uf->parent[x];
}

static void uf_union(UnionFind *uf, uint32_t x, uint32_t y) {
    uint32_t rx = uf_find(uf, x);
    uint32_t ry = uf_find(uf, y);
    if (rx == ry) return;

    if (uf->rank[rx] < uf->rank[ry]) {
        uf->parent[rx] = ry;
    } else if (uf->rank[rx] > uf->rank[ry]) {
        uf->parent[ry] = rx;
    } else {
        uf->parent[ry] = rx;
        uf->rank[rx]++;
    }
}

/* Count shared k-mers between two monomers using bitvectors */
static uint32_t count_shared_kmers(const Monomer *a, const Monomer *b) {
    return bitvec_intersection_count(&a->bitvec, &b->bitvec);
}

/* Jaccard threshold for component membership */
#define COMPONENT_JACCARD_THRESHOLD 0.8f

/* Greedy clustering driven by large buckets
 *
 * Algorithm:
 * 1. Sort buckets by size (largest first)
 * 2. For each bucket:
 *    - Take first unassigned monomer as seed
 *    - Linear scan: find all with Jaccard >= 80%
 *    - If matched monomer has cluster -> merge via Union-Find
 *    - Continue until bucket exhausted
 */
void find_components(ClassifyContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "Finding components (greedy, Jaccard >= %.0f%%)...\n",
                COMPONENT_JACCARD_THRESHOLD * 100);
    }

    /* Build mapping: monomer_id -> local index for Union-Find */
    uint32_t n_monomers = 0;
    uint32_t *mono_indices = malloc(ctx->monomer_count * sizeof(uint32_t));
    uint32_t *mono_to_idx = malloc(ctx->monomer_count * sizeof(uint32_t));

    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        mono_to_idx[i] = UINT32_MAX;
        if (ctx->skip_flanks && ctx->monomers[i].type != TYPE_MONOMER) {
            continue;
        }
        mono_to_idx[i] = n_monomers;
        mono_indices[n_monomers++] = i;
    }

    if (n_monomers == 0) {
        free(mono_indices);
        free(mono_to_idx);
        ctx->components = NULL;
        ctx->component_count = 0;
        return;
    }

    UnionFind *uf = uf_create(n_monomers);

    /* Sort buckets by size (largest first) */
    typedef struct { uint16_t kmer; uint32_t count; } BucketInfo;
    BucketInfo *sorted_buckets = malloc(KMER_BUCKETS * sizeof(BucketInfo));
    uint32_t n_useful_buckets = 0;

    for (int k = 0; k < KMER_BUCKETS; k++) {
        if (ctx->buckets[k].count >= 2) {
            sorted_buckets[n_useful_buckets].kmer = (uint16_t)k;
            sorted_buckets[n_useful_buckets].count = ctx->buckets[k].count;
            n_useful_buckets++;
        }
    }

    /* Sort by count ASCENDING - rare k-mers first (more informative) */
    for (uint32_t i = 0; i < n_useful_buckets; i++) {
        for (uint32_t j = i + 1; j < n_useful_buckets; j++) {
            if (sorted_buckets[j].count < sorted_buckets[i].count) {
                BucketInfo tmp = sorted_buckets[i];
                sorted_buckets[i] = sorted_buckets[j];
                sorted_buckets[j] = tmp;
            }
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "  Processing %u buckets (SMALLEST first - rare k-mers), %u monomers...\n",
                n_useful_buckets, n_monomers);
        fprintf(stderr, "  First 5 buckets: %u, %u, %u, %u, %u monomers\n",
                sorted_buckets[0].count, sorted_buckets[1].count,
                sorted_buckets[2].count, sorted_buckets[3].count,
                sorted_buckets[4].count);
        fprintf(stderr, "  Last 5 buckets: %u, %u, %u, %u, %u monomers\n",
                sorted_buckets[n_useful_buckets-5].count,
                sorted_buckets[n_useful_buckets-4].count,
                sorted_buckets[n_useful_buckets-3].count,
                sorted_buckets[n_useful_buckets-2].count,
                sorted_buckets[n_useful_buckets-1].count);
    }

    uint64_t pairs_checked = 0;
    uint32_t merges = 0;
    uint32_t progress_step = n_useful_buckets / 20;
    if (progress_step == 0) progress_step = 1;

    /* Global "processed" flag - once a monomer is seed in any bucket, don't use as seed again */
    bool *processed = calloc(n_monomers, sizeof(bool));

    /* Process buckets from largest to smallest */
    for (uint32_t b = 0; b < n_useful_buckets; b++) {
        uint16_t kmer = sorted_buckets[b].kmer;
        KmerBucket *bucket = &ctx->buckets[kmer];

        /* Progress: show bucket info for large buckets */
        if (ctx->verbose) {
            if (bucket->count >= 1000) {
                fprintf(stderr, "\r  Bucket %u/%u: %u monomers (kmer %u)...                    \n",
                        b + 1, n_useful_buckets, bucket->count, kmer);
            } else if (b % progress_step == 0) {
                uint32_t pct = (b * 100) / n_useful_buckets;
                fprintf(stderr, "\r  Progress: %u%% (%u/%u buckets, %llu pairs, %u merges)       ",
                        pct, b, n_useful_buckets, (unsigned long long)pairs_checked, merges);
            }
        }

        uint32_t seeds_in_bucket = 0;
        uint32_t bucket_pairs = 0;

        /* For each monomer in bucket, compare with all others */
        for (uint32_t i = 0; i < bucket->count; i++) {
            uint32_t mono_i = bucket->hits[i].monomer_id;
            uint32_t idx_i = mono_to_idx[mono_i];
            if (idx_i == UINT32_MAX) continue;

            /* Skip if already processed as seed - but still compare others with existing clusters */
            if (processed[idx_i]) continue;

            seeds_in_bucket++;

            /* Progress within large bucket */
            if (ctx->verbose && bucket->count >= 1000 && seeds_in_bucket % 100 == 0) {
                fprintf(stderr, "\r    Seed %u, %u pairs checked, %u merges       ",
                        seeds_in_bucket, bucket_pairs, merges);
            }

            /* Compare with ALL other monomers in bucket (including processed ones) */
            for (uint32_t j = 0; j < bucket->count; j++) {
                if (j == i) continue;

                uint32_t mono_j = bucket->hits[j].monomer_id;
                uint32_t idx_j = mono_to_idx[mono_j];
                if (idx_j == UINT32_MAX) continue;

                /* Skip if already in same component */
                if (uf_find(uf, idx_i) == uf_find(uf, idx_j)) continue;

                pairs_checked++;
                bucket_pairs++;

                /* Check Jaccard */
                float jaccard = bitvec_jaccard(&ctx->monomers[mono_i].bitvec,
                                               &ctx->monomers[mono_j].bitvec);

                if (jaccard >= COMPONENT_JACCARD_THRESHOLD) {
                    uf_union(uf, idx_i, idx_j);
                    merges++;
                }
            }

            /* Mark i as processed */
            processed[idx_i] = true;
        }

        if (ctx->verbose && bucket->count >= 1000) {
            fprintf(stderr, "\r    Done: %u seeds, %u pairs, %u total merges\n",
                    seeds_in_bucket, bucket_pairs, merges);
        }
    }

    free(processed);

    free(sorted_buckets);

    if (ctx->verbose) {
        fprintf(stderr, "\r  Progress: 100%% (%u/%u buckets, %llu pairs, %u merges)\n",
                n_useful_buckets, n_useful_buckets, (unsigned long long)pairs_checked, merges);
    }

    /* Extract components from Union-Find */
    uint32_t *root_to_comp = calloc(n_monomers, sizeof(uint32_t));
    bool *root_seen = calloc(n_monomers, sizeof(bool));
    uint32_t n_components = 0;

    /* First pass: count components and assign IDs */
    for (uint32_t i = 0; i < n_monomers; i++) {
        uint32_t root = uf_find(uf, i);
        if (!root_seen[root]) {
            root_to_comp[root] = n_components++;
            root_seen[root] = true;
        }
    }

    /* Allocate components */
    ctx->components = malloc(n_components * sizeof(Component));
    ctx->component_count = n_components;
    ctx->component_capacity = n_components;

    for (uint32_t c = 0; c < n_components; c++) {
        ctx->components[c].member_ids = NULL;
        ctx->components[c].member_count = 0;
        ctx->components[c].member_capacity = 0;
    }

    /* Second pass: count members per component */
    uint32_t *comp_sizes = calloc(n_components, sizeof(uint32_t));
    for (uint32_t i = 0; i < n_monomers; i++) {
        uint32_t root = uf_find(uf, i);
        uint32_t comp_id = root_to_comp[root];
        comp_sizes[comp_id]++;
    }

    /* Allocate member arrays */
    for (uint32_t c = 0; c < n_components; c++) {
        ctx->components[c].member_ids = malloc(comp_sizes[c] * sizeof(uint32_t));
        ctx->components[c].member_capacity = comp_sizes[c];
    }

    /* Third pass: fill member arrays */
    memset(comp_sizes, 0, n_components * sizeof(uint32_t));
    for (uint32_t i = 0; i < n_monomers; i++) {
        uint32_t root = uf_find(uf, i);
        uint32_t comp_id = root_to_comp[root];
        uint32_t mono_id = mono_indices[i];

        ctx->components[comp_id].member_ids[comp_sizes[comp_id]++] = mono_id;
        ctx->monomers[mono_id].component_id = comp_id;
    }

    for (uint32_t c = 0; c < n_components; c++) {
        ctx->components[c].member_count = comp_sizes[c];
        ctx->components[c].representative_id = ctx->components[c].member_ids[0];
    }

    free(comp_sizes);
    free(root_to_comp);
    free(root_seen);
    free(mono_to_idx);
    free(mono_indices);
    uf_free(uf);

    if (ctx->verbose) {
        fprintf(stderr, "  Found %u components\n", n_components);

        /* Size distribution */
        uint32_t singles = 0, small = 0, medium = 0, large = 0;
        for (uint32_t c = 0; c < n_components; c++) {
            uint32_t size = ctx->components[c].member_count;
            if (size == 1) singles++;
            else if (size <= 10) small++;
            else if (size <= 100) medium++;
            else large++;
        }
        fprintf(stderr, "  Distribution: %u singletons, %u small (2-10), "
                "%u medium (11-100), %u large (>100)\n",
                singles, small, medium, large);
    }
}

/* Find subcomponents within a single component using Jaccard threshold */
static void cluster_component(ClassifyContext *ctx, uint32_t comp_id) {
    Component *comp = &ctx->components[comp_id];

    if (comp->member_count <= 1) {
        /* Single member = single subcomponent */
        uint32_t mono_id = comp->member_ids[0];

        if (ctx->subcomponent_count >= ctx->subcomponent_capacity) {
            ctx->subcomponent_capacity = ctx->subcomponent_capacity == 0 ? 64 : ctx->subcomponent_capacity * 2;
            ctx->subcomponents = realloc(ctx->subcomponents,
                                          ctx->subcomponent_capacity * sizeof(Subcomponent));
        }

        Subcomponent *sub = &ctx->subcomponents[ctx->subcomponent_count];
        sub->component_id = comp_id;
        sub->member_ids = malloc(sizeof(uint32_t));
        sub->member_ids[0] = mono_id;
        sub->member_count = 1;
        sub->member_capacity = 1;
        sub->representative_id = mono_id;
        sub->avg_jaccard = 1.0f;

        ctx->monomers[mono_id].subcomponent_id = ctx->subcomponent_count;
        ctx->monomers[mono_id].jaccard_to_rep = 1.0f;

        ctx->subcomponent_count++;
        return;
    }

    /* Use Union-Find to cluster by Jaccard */
    UnionFind *uf = uf_create(comp->member_count);

    /* Compare all pairs (within component, should be manageable) */
    for (uint32_t i = 0; i < comp->member_count; i++) {
        uint32_t mono_i = comp->member_ids[i];

        for (uint32_t j = i + 1; j < comp->member_count; j++) {
            uint32_t mono_j = comp->member_ids[j];

            /* Skip if already merged */
            if (uf_find(uf, i) == uf_find(uf, j)) continue;

            float jaccard = bitvec_jaccard(&ctx->monomers[mono_i].bitvec,
                                           &ctx->monomers[mono_j].bitvec);

            if (jaccard >= ctx->jaccard_threshold) {
                uf_union(uf, i, j);
            }
        }
    }

    /* Extract subcomponents */
    uint32_t *root_to_sub = calloc(comp->member_count, sizeof(uint32_t));
    bool *root_seen = calloc(comp->member_count, sizeof(bool));
    uint32_t n_subs = 0;

    for (uint32_t i = 0; i < comp->member_count; i++) {
        uint32_t root = uf_find(uf, i);
        if (!root_seen[root]) {
            root_to_sub[root] = n_subs++;
            root_seen[root] = true;
        }
    }

    /* Create subcomponents */
    uint32_t *sub_sizes = calloc(n_subs, sizeof(uint32_t));
    for (uint32_t i = 0; i < comp->member_count; i++) {
        uint32_t root = uf_find(uf, i);
        sub_sizes[root_to_sub[root]]++;
    }

    uint32_t first_sub_id = ctx->subcomponent_count;

    for (uint32_t s = 0; s < n_subs; s++) {
        if (ctx->subcomponent_count >= ctx->subcomponent_capacity) {
            ctx->subcomponent_capacity = ctx->subcomponent_capacity == 0 ? 64 : ctx->subcomponent_capacity * 2;
            ctx->subcomponents = realloc(ctx->subcomponents,
                                          ctx->subcomponent_capacity * sizeof(Subcomponent));
        }

        Subcomponent *sub = &ctx->subcomponents[ctx->subcomponent_count];
        sub->component_id = comp_id;
        sub->member_ids = malloc(sub_sizes[s] * sizeof(uint32_t));
        sub->member_count = 0;
        sub->member_capacity = sub_sizes[s];
        sub->avg_jaccard = 0.0f;

        ctx->subcomponent_count++;
    }

    /* Fill member arrays and find representatives */
    for (uint32_t i = 0; i < comp->member_count; i++) {
        uint32_t root = uf_find(uf, i);
        uint32_t sub_local = root_to_sub[root];
        uint32_t sub_id = first_sub_id + sub_local;
        uint32_t mono_id = comp->member_ids[i];

        Subcomponent *sub = &ctx->subcomponents[sub_id];
        sub->member_ids[sub->member_count++] = mono_id;
        ctx->monomers[mono_id].subcomponent_id = sub_id;
    }

    /* Set representatives (first member) and compute jaccard to rep */
    for (uint32_t s = 0; s < n_subs; s++) {
        uint32_t sub_id = first_sub_id + s;
        Subcomponent *sub = &ctx->subcomponents[sub_id];
        sub->representative_id = sub->member_ids[0];

        float total_jaccard = 0.0f;
        for (uint32_t i = 0; i < sub->member_count; i++) {
            uint32_t mono_id = sub->member_ids[i];
            float jac = bitvec_jaccard(&ctx->monomers[mono_id].bitvec,
                                       &ctx->monomers[sub->representative_id].bitvec);
            ctx->monomers[mono_id].jaccard_to_rep = jac;
            total_jaccard += jac;
        }
        sub->avg_jaccard = total_jaccard / sub->member_count;
    }

    free(sub_sizes);
    free(root_to_sub);
    free(root_seen);
    uf_free(uf);
}

/* Find subcomponents within each component */
void find_subcomponents(ClassifyContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "Finding subcomponents (Jaccard >= %.0f%%)...\n",
                ctx->jaccard_threshold * 100);
    }

    ctx->subcomponents = NULL;
    ctx->subcomponent_count = 0;
    ctx->subcomponent_capacity = 0;

    uint32_t progress_step = ctx->component_count / 20;
    if (progress_step == 0) progress_step = 1;

    for (uint32_t c = 0; c < ctx->component_count; c++) {
        if (ctx->verbose && c % progress_step == 0) {
            uint32_t pct = (c * 100) / ctx->component_count;
            fprintf(stderr, "\r  Progress: %u%% (%u/%u components, %u subcomponents)",
                    pct, c, ctx->component_count, ctx->subcomponent_count);
        }

        /* Skip very large components - treat as single subcomponent */
        if (ctx->components[c].member_count > 5000) {
            Component *comp = &ctx->components[c];

            if (ctx->subcomponent_count >= ctx->subcomponent_capacity) {
                ctx->subcomponent_capacity = ctx->subcomponent_capacity == 0 ? 64 : ctx->subcomponent_capacity * 2;
                ctx->subcomponents = realloc(ctx->subcomponents,
                                              ctx->subcomponent_capacity * sizeof(Subcomponent));
            }

            Subcomponent *sub = &ctx->subcomponents[ctx->subcomponent_count];
            sub->component_id = c;
            sub->member_ids = malloc(comp->member_count * sizeof(uint32_t));
            memcpy(sub->member_ids, comp->member_ids, comp->member_count * sizeof(uint32_t));
            sub->member_count = comp->member_count;
            sub->member_capacity = comp->member_count;
            sub->representative_id = comp->member_ids[0];
            sub->avg_jaccard = 0.0f;

            for (uint32_t i = 0; i < comp->member_count; i++) {
                ctx->monomers[comp->member_ids[i]].subcomponent_id = ctx->subcomponent_count;
                ctx->monomers[comp->member_ids[i]].jaccard_to_rep = 0.0f;
            }

            ctx->subcomponent_count++;
            continue;
        }

        cluster_component(ctx, c);
    }

    if (ctx->verbose) {
        fprintf(stderr, "\r  Progress: 100%% (%u/%u components, %u subcomponents)\n",
                ctx->component_count, ctx->component_count, ctx->subcomponent_count);
    }
}

/* Main clustering entry point */
void cluster_monomers(ClassifyContext *ctx) {
    find_components(ctx);
    find_subcomponents(ctx);
}
