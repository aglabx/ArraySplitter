/*
 * minhash_cluster.c - Cascade clustering with MinHash/LSH
 *
 * Three-level hierarchical clustering:
 * 1. Subfamilies (threshold ~0.8) - nearly identical sequences
 * 2. Families (threshold ~0.5) - cluster representatives
 * 3. Superfamilies (threshold ~0.2) - distant homology
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include "minhash.h"

/* Forward declaration */
void lsh_build(MHContext *ctx, uint32_t *ids, uint32_t count);

/*
 * Select centroid for a cluster
 *
 * Centroid = member with highest average similarity to all others
 */
uint32_t select_centroid(MHContext *ctx, uint32_t *members, uint32_t count) {
    if (count == 1) {
        return members[0];
    }

    uint32_t best_id = members[0];
    float best_score = -1.0f;

    for (uint32_t i = 0; i < count; i++) {
        uint32_t id_i = members[i];
        float sum_sim = 0.0f;

        for (uint32_t j = 0; j < count; j++) {
            if (i == j) continue;
            uint32_t id_j = members[j];
            sum_sim += minhash_similarity(&ctx->monomers[id_i].signature,
                                          &ctx->monomers[id_j].signature);
        }

        float avg_sim = sum_sim / (count - 1);
        if (avg_sim > best_score) {
            best_score = avg_sim;
            best_id = id_i;
        }
    }

    return best_id;
}

/*
 * Greedy Cover clustering (MMseqs2/Linclust style)
 *
 * Instead of all-pairs comparison within LSH buckets (O(S²)),
 * we use star topology: pick center, compare others only to center (O(S))
 */
void cluster_level_greedy(MHContext *ctx, float threshold, uint32_t *ids, uint32_t count,
                          MHCluster **out_clusters, uint32_t *out_count) {
    if (count == 0) {
        *out_clusters = NULL;
        *out_count = 0;
        return;
    }

    /* Build LSH index for this subset */
    lsh_build(ctx, ids, count);

    if (ctx->verbose) {
        fprintf(stderr, "  Greedy clustering %u items at threshold %.2f...\n", count, threshold);
    }

    /* Covered flag: true if monomer already assigned to a cluster */
    bool *is_covered = calloc(ctx->monomer_count, sizeof(bool));

    /* Sort by length (longer = better center candidates) */
    /* Create sorted index array */
    uint32_t *sorted_ids = malloc(count * sizeof(uint32_t));
    memcpy(sorted_ids, ids, count * sizeof(uint32_t));

    /* Simple insertion sort by length descending (stable for equal lengths) */
    for (uint32_t i = 1; i < count; i++) {
        uint32_t key = sorted_ids[i];
        uint32_t key_len = ctx->monomers[key].length;
        int j = i - 1;
        while (j >= 0 && ctx->monomers[sorted_ids[j]].length < key_len) {
            sorted_ids[j + 1] = sorted_ids[j];
            j--;
        }
        sorted_ids[j + 1] = key;
    }

    /* Allocate clusters array */
    uint32_t cluster_capacity = 1024;
    MHCluster *clusters = malloc(cluster_capacity * sizeof(MHCluster));
    uint32_t num_clusters = 0;

    uint32_t total_compared = 0;
    uint32_t total_assigned = 0;
    time_t last_report = time(NULL);

    /* Greedy cover: iterate by quality (length) */
    for (uint32_t i = 0; i < count; i++) {
        uint32_t center_id = sorted_ids[i];

        /* Skip if already covered */
        if (is_covered[center_id]) continue;

        /* Create new cluster with this center */
        if (num_clusters >= cluster_capacity) {
            cluster_capacity *= 2;
            clusters = realloc(clusters, cluster_capacity * sizeof(MHCluster));
        }

        MHCluster *c = &clusters[num_clusters];
        c->id = num_clusters;
        c->level = 1;
        c->member_capacity = 64;
        c->members = malloc(c->member_capacity * sizeof(uint32_t));
        c->member_count = 1;
        c->members[0] = center_id;
        c->representative_id = center_id;

        is_covered[center_id] = true;
        total_assigned++;

        /* Query LSH for candidates */
        MinHashSig *center_sig = &ctx->monomers[center_id].signature;
        uint32_t cand_count = 0;
        uint32_t *candidates = lsh_query_candidates(&ctx->lsh, center_sig, &cand_count);

        /* Check each candidate against center only (star topology) */
        for (uint32_t j = 0; j < cand_count; j++) {
            uint32_t cand_id = candidates[j];

            /* Skip self and already covered */
            if (cand_id == center_id || is_covered[cand_id]) continue;

            total_compared++;

            /* Quick pre-filter: check first 16 hashes */
            uint32_t quick_matches = 0;
            for (int h = 0; h < 16; h++) {
                if (center_sig->hashes[h] == ctx->monomers[cand_id].signature.hashes[h]) {
                    quick_matches++;
                }
            }
            /* If <10 of first 16 match, full Jaccard likely < 0.6 */
            if (quick_matches < (uint32_t)(threshold * 16 * 0.8)) continue;

            /* Full similarity check */
            float sim = minhash_similarity(center_sig, &ctx->monomers[cand_id].signature);

            if (sim >= threshold) {
                /* Add to this cluster */
                if (c->member_count >= c->member_capacity) {
                    c->member_capacity *= 2;
                    c->members = realloc(c->members, c->member_capacity * sizeof(uint32_t));
                }
                c->members[c->member_count++] = cand_id;
                is_covered[cand_id] = true;
                total_assigned++;
            }
        }

        free(candidates);
        num_clusters++;

        /* Progress: report every second */
        if (ctx->verbose) {
            time_t now = time(NULL);
            if (now > last_report || total_assigned == count) {
                last_report = now;
                fprintf(stderr, "\r    Assigned: %u/%u (%.1f%%) | Clusters: %u | Comparisons: %u    ",
                        total_assigned, count, 100.0 * total_assigned / count,
                        num_clusters, total_compared);
                fflush(stderr);
            }
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "\r    Comparisons: %u | Clusters: %u                    \n",
                total_compared, num_clusters);
    }

    free(is_covered);
    free(sorted_ids);

    /* === Post-processing: Merge similar clusters (transitive closure) === */
    if (num_clusters > 1 && ctx->verbose) {
        fprintf(stderr, "  Merging similar clusters...\n");
    }

    /* Use Union-Find on clusters */
    uint32_t *cluster_parent = malloc(num_clusters * sizeof(uint32_t));
    for (uint32_t i = 0; i < num_clusters; i++) {
        cluster_parent[i] = i;
    }

    /* Find with path compression */
    #define CLUSTER_FIND(x) ({ \
        uint32_t _x = (x); \
        while (cluster_parent[_x] != _x) { \
            cluster_parent[_x] = cluster_parent[cluster_parent[_x]]; \
            _x = cluster_parent[_x]; \
        } \
        _x; \
    })

    /* Build LSH index of cluster representatives */
    uint32_t *rep_ids = malloc(num_clusters * sizeof(uint32_t));
    for (uint32_t i = 0; i < num_clusters; i++) {
        rep_ids[i] = clusters[i].representative_id;
    }
    lsh_build(ctx, rep_ids, num_clusters);

    uint32_t merges = 0;
    for (uint32_t i = 0; i < num_clusters; i++) {
        uint32_t rep_id = clusters[i].representative_id;
        MinHashSig *sig = &ctx->monomers[rep_id].signature;

        uint32_t cand_count = 0;
        uint32_t *candidates = lsh_query_candidates(&ctx->lsh, sig, &cand_count);

        for (uint32_t j = 0; j < cand_count; j++) {
            uint32_t cand_rep = candidates[j];
            if (cand_rep == rep_id) continue;

            /* Find which cluster this candidate belongs to */
            uint32_t cand_cluster = UINT32_MAX;
            for (uint32_t k = 0; k < num_clusters; k++) {
                if (clusters[k].representative_id == cand_rep) {
                    cand_cluster = k;
                    break;
                }
            }
            if (cand_cluster == UINT32_MAX) continue;

            uint32_t root_i = CLUSTER_FIND(i);
            uint32_t root_j = CLUSTER_FIND(cand_cluster);
            if (root_i == root_j) continue;

            float sim = minhash_similarity(sig, &ctx->monomers[cand_rep].signature);
            if (sim >= threshold) {
                cluster_parent[root_j] = root_i;
                merges++;
            }
        }
        free(candidates);
    }

    if (ctx->verbose && merges > 0) {
        fprintf(stderr, "    Merged %u cluster pairs\n", merges);
    }

    /* Rebuild clusters after merging */
    if (merges > 0) {
        /* Count final clusters */
        uint32_t final_count = 0;
        uint32_t *new_id = malloc(num_clusters * sizeof(uint32_t));
        for (uint32_t i = 0; i < num_clusters; i++) {
            if (CLUSTER_FIND(i) == i) {
                new_id[i] = final_count++;
            }
        }
        for (uint32_t i = 0; i < num_clusters; i++) {
            new_id[i] = new_id[CLUSTER_FIND(i)];
        }

        /* Create merged clusters */
        MHCluster *merged = calloc(final_count, sizeof(MHCluster));
        for (uint32_t i = 0; i < final_count; i++) {
            merged[i].id = i;
            merged[i].level = 1;
            merged[i].member_capacity = 64;
            merged[i].members = malloc(64 * sizeof(uint32_t));
            merged[i].member_count = 0;
        }

        /* Assign members to merged clusters */
        for (uint32_t i = 0; i < num_clusters; i++) {
            uint32_t dest = new_id[i];
            MHCluster *src = &clusters[i];
            MHCluster *dst = &merged[dest];

            for (uint32_t m = 0; m < src->member_count; m++) {
                if (dst->member_count >= dst->member_capacity) {
                    dst->member_capacity *= 2;
                    dst->members = realloc(dst->members, dst->member_capacity * sizeof(uint32_t));
                }
                dst->members[dst->member_count++] = src->members[m];
            }

            /* Keep representative from largest original cluster */
            if (CLUSTER_FIND(i) == i) {
                dst->representative_id = src->representative_id;
            }

            free(src->members);
        }
        free(clusters);
        free(new_id);

        clusters = merged;
        num_clusters = final_count;

        if (ctx->verbose) {
            fprintf(stderr, "    After merge: %u clusters\n", num_clusters);
        }
    }

    free(cluster_parent);
    free(rep_ids);
    #undef CLUSTER_FIND

    *out_clusters = clusters;
    *out_count = num_clusters;
}

/* Legacy function for compatibility */
void cluster_level(MHContext *ctx, float threshold, uint32_t *ids, uint32_t count) {
    (void)ctx; (void)threshold; (void)ids; (void)count;
    /* Deprecated - use cluster_level_greedy instead */
}

/*
 * Extract clusters from Union-Find
 * Returns array of cluster structs
 */
static MHCluster* extract_clusters(MHContext *ctx, uint32_t *ids, uint32_t count,
                                   uint32_t *cluster_count, uint32_t level) {
    /* Count unique roots */
    uint32_t *roots = malloc(count * sizeof(uint32_t));
    uint32_t *root_idx = malloc(ctx->monomer_count * sizeof(uint32_t));
    memset(root_idx, 0xFF, ctx->monomer_count * sizeof(uint32_t));

    uint32_t num_clusters = 0;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t root = uf_find(&ctx->uf, ids[i]);
        if (root_idx[root] == 0xFFFFFFFF) {
            root_idx[root] = num_clusters;
            roots[num_clusters++] = root;
        }
    }

    /* Allocate clusters */
    MHCluster *clusters = malloc(num_clusters * sizeof(MHCluster));
    for (uint32_t i = 0; i < num_clusters; i++) {
        clusters[i].id = i;
        clusters[i].level = level;
        clusters[i].member_count = 0;
        clusters[i].member_capacity = 16;
        clusters[i].members = malloc(16 * sizeof(uint32_t));
        clusters[i].representative_id = 0;
        clusters[i].avg_similarity = 0.0f;
    }

    /* Assign members to clusters */
    for (uint32_t i = 0; i < count; i++) {
        uint32_t id = ids[i];
        uint32_t root = uf_find(&ctx->uf, id);
        uint32_t ci = root_idx[root];

        MHCluster *c = &clusters[ci];
        if (c->member_count >= c->member_capacity) {
            c->member_capacity *= 2;
            c->members = realloc(c->members, c->member_capacity * sizeof(uint32_t));
        }
        c->members[c->member_count++] = id;
    }

    /* Select centroids */
    for (uint32_t i = 0; i < num_clusters; i++) {
        clusters[i].representative_id = select_centroid(ctx,
            clusters[i].members, clusters[i].member_count);
    }

    free(roots);
    free(root_idx);

    *cluster_count = num_clusters;
    return clusters;
}

/*
 * Level 1: Find subfamilies (high stringency) - Greedy Cover
 */
void find_subfamilies(MHContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "\n=== Level 1: Subfamilies (threshold %.2f) ===\n",
                ctx->level1_threshold);
    }

    /* All monomers participate */
    uint32_t *all_ids = malloc(ctx->monomer_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        all_ids[i] = i;
    }

    /* Greedy clustering */
    cluster_level_greedy(ctx, ctx->level1_threshold, all_ids, ctx->monomer_count,
                         &ctx->subfamilies, &ctx->subfamily_count);

    /* Assign subfamily IDs to monomers */
    for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
        MHCluster *c = &ctx->subfamilies[i];
        for (uint32_t j = 0; j < c->member_count; j++) {
            ctx->monomers[c->members[j]].subfamily_id = i;
        }
    }

    if (ctx->verbose) {
        uint32_t singletons = 0;
        for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
            if (ctx->subfamilies[i].member_count == 1) singletons++;
        }
        fprintf(stderr, "  Result: %u subfamilies (%u singletons)\n",
                ctx->subfamily_count, singletons);
    }

    free(all_ids);
}

/*
 * Level 2: Find families (medium stringency) - Greedy Cover
 * Only cluster representatives from Level 1
 */
void find_families(MHContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "\n=== Level 2: Families (threshold %.2f) ===\n",
                ctx->level2_threshold);
    }

    /* Collect representatives */
    uint32_t *rep_ids = malloc(ctx->subfamily_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
        rep_ids[i] = ctx->subfamilies[i].representative_id;
    }

    /* Greedy clustering of representatives */
    cluster_level_greedy(ctx, ctx->level2_threshold, rep_ids, ctx->subfamily_count,
                         &ctx->families, &ctx->family_count);

    /* Build reverse mapping: rep_id -> subfamily index */
    uint32_t *rep_to_subfamily = calloc(ctx->monomer_count, sizeof(uint32_t));
    for (uint32_t k = 0; k < ctx->subfamily_count; k++) {
        rep_to_subfamily[ctx->subfamilies[k].representative_id] = k;
    }

    /* Map subfamilies to families, assign family_id to monomers */
    for (uint32_t i = 0; i < ctx->family_count; i++) {
        MHCluster *fam = &ctx->families[i];
        for (uint32_t j = 0; j < fam->member_count; j++) {
            uint32_t rep_id = fam->members[j];
            uint32_t sf_idx = rep_to_subfamily[rep_id];
            MHCluster *sf = &ctx->subfamilies[sf_idx];
            for (uint32_t m = 0; m < sf->member_count; m++) {
                ctx->monomers[sf->members[m]].family_id = i;
            }
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "  Result: %u families\n", ctx->family_count);
    }

    free(rep_to_subfamily);
    free(rep_ids);
}

/*
 * Level 3: Find superfamilies (low stringency)
 * Only cluster representatives from Level 2
 */
void find_superfamilies(MHContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "\n=== Level 3: Superfamilies (threshold %.2f) ===\n",
                ctx->level3_threshold);
    }

    /* Collect family representatives */
    uint32_t *fam_reps = malloc(ctx->family_count * sizeof(uint32_t));
    for (uint32_t i = 0; i < ctx->family_count; i++) {
        fam_reps[i] = ctx->families[i].representative_id;
    }

    /* For small count, do all-pairs; otherwise use greedy */
    if (ctx->family_count < 5000) {
        if (ctx->verbose) {
            fprintf(stderr, "  Using all-pairs comparison for %u families\n",
                    ctx->family_count);
        }

        /* Initialize Union-Find */
        uf_init(&ctx->uf, ctx->monomer_count);

        /* All-pairs comparison */
        uint32_t pairs_merged = 0;
        for (uint32_t i = 0; i < ctx->family_count; i++) {
            for (uint32_t j = i + 1; j < ctx->family_count; j++) {
                float sim = minhash_similarity(
                    &ctx->monomers[fam_reps[i]].signature,
                    &ctx->monomers[fam_reps[j]].signature);
                if (sim >= ctx->level3_threshold) {
                    uf_union(&ctx->uf, fam_reps[i], fam_reps[j]);
                    pairs_merged++;
                }
            }
        }

        /* Extract superfamilies using Union-Find */
        ctx->superfamilies = extract_clusters(ctx, fam_reps, ctx->family_count,
                                              &ctx->superfamily_count, 3);
        uf_free(&ctx->uf);

        if (ctx->verbose) {
            fprintf(stderr, "  Merged: %u pairs\n", pairs_merged);
        }
    } else {
        /* Use greedy for large counts */
        cluster_level_greedy(ctx, ctx->level3_threshold, fam_reps, ctx->family_count,
                             &ctx->superfamilies, &ctx->superfamily_count);
    }

    /* Build reverse mapping: family rep -> family index */
    uint32_t *rep_to_family = calloc(ctx->monomer_count, sizeof(uint32_t));
    for (uint32_t k = 0; k < ctx->family_count; k++) {
        rep_to_family[ctx->families[k].representative_id] = k;
    }

    /* Map superfamilies to monomers via family -> subfamily chain */
    for (uint32_t i = 0; i < ctx->superfamily_count; i++) {
        MHCluster *sf = &ctx->superfamilies[i];
        for (uint32_t j = 0; j < sf->member_count; j++) {
            uint32_t fam_rep_id = sf->members[j];
            uint32_t fam_idx = rep_to_family[fam_rep_id];

            /* All monomers with this family_id get this superfamily_id */
            for (uint32_t m = 0; m < ctx->monomer_count; m++) {
                if (ctx->monomers[m].family_id == fam_idx) {
                    ctx->monomers[m].superfamily_id = i;
                }
            }
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "  Result: %u superfamilies\n", ctx->superfamily_count);
    }

    free(rep_to_family);
    free(fam_reps);
}
