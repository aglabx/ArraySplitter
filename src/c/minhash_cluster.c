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

/* ==================== Array Grouping ==================== */

/* Group of monomers belonging to one array */
typedef struct {
    char *array_id;           /* sequence_id string */
    uint32_t *monomer_ids;    /* monomer IDs in this array */
    uint32_t count;
    uint32_t capacity;
} ArrayGroup;

/* Collection of all array groups */
typedef struct {
    ArrayGroup *groups;
    uint32_t count;
    uint32_t capacity;
} ArrayGroups;

/* Initialize array groups */
static void array_groups_init(ArrayGroups *ag) {
    ag->capacity = 256;
    ag->groups = malloc(ag->capacity * sizeof(ArrayGroup));
    ag->count = 0;
}

/* Free array groups */
static void array_groups_free(ArrayGroups *ag) {
    for (uint32_t i = 0; i < ag->count; i++) {
        free(ag->groups[i].array_id);
        free(ag->groups[i].monomer_ids);
    }
    free(ag->groups);
}

/* Find or create group for array_id */
static ArrayGroup* array_groups_get(ArrayGroups *ag, const char *array_id) {
    /* Linear search (could use hash table for many arrays) */
    for (uint32_t i = 0; i < ag->count; i++) {
        if (strcmp(ag->groups[i].array_id, array_id) == 0) {
            return &ag->groups[i];
        }
    }

    /* Create new group */
    if (ag->count >= ag->capacity) {
        ag->capacity *= 2;
        ag->groups = realloc(ag->groups, ag->capacity * sizeof(ArrayGroup));
    }

    ArrayGroup *g = &ag->groups[ag->count++];
    g->array_id = strdup(array_id);
    g->capacity = 64;
    g->monomer_ids = malloc(g->capacity * sizeof(uint32_t));
    g->count = 0;

    return g;
}

/* Add monomer to its array group */
static void array_groups_add(ArrayGroups *ag, const char *array_id, uint32_t monomer_id) {
    ArrayGroup *g = array_groups_get(ag, array_id);

    if (g->count >= g->capacity) {
        g->capacity *= 2;
        g->monomer_ids = realloc(g->monomer_ids, g->capacity * sizeof(uint32_t));
    }
    g->monomer_ids[g->count++] = monomer_id;
}

/* Group monomers by their sequence_id (array) */
static void group_monomers_by_array(MHContext *ctx, uint32_t *ids, uint32_t count, ArrayGroups *ag) {
    array_groups_init(ag);

    for (uint32_t i = 0; i < count; i++) {
        uint32_t id = ids[i];
        const char *array_id = ctx->monomers[id].sequence_id;
        array_groups_add(ag, array_id, id);
    }
}

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
 * Array-local clustering: cluster each array independently
 *
 * This respects the biological reality that monomers within an array
 * are more related (concerted evolution) and should cluster together.
 * Cross-array clustering happens at higher levels (families, superfamilies).
 */
void cluster_array_local(MHContext *ctx, float threshold, uint32_t *ids, uint32_t count,
                         MHCluster **out_clusters, uint32_t *out_count) {
    if (count == 0) {
        *out_clusters = NULL;
        *out_count = 0;
        return;
    }

    /* Group monomers by array */
    ArrayGroups ag;
    group_monomers_by_array(ctx, ids, count, &ag);

    if (ctx->verbose) {
        fprintf(stderr, "  Array-local clustering: %u arrays, %u monomers\n", ag.count, count);
    }

    /* Allocate result arrays */
    uint32_t total_clusters = 0;
    uint32_t clusters_capacity = 1024;
    MHCluster *all_clusters = malloc(clusters_capacity * sizeof(MHCluster));

    /* Process each array */
    for (uint32_t a = 0; a < ag.count; a++) {
        ArrayGroup *g = &ag.groups[a];

        MHCluster *array_clusters = NULL;
        uint32_t array_cluster_count = 0;

        if (g->count < 3) {
            /* Small array: put all in one cluster */
            array_cluster_count = 1;
            array_clusters = malloc(sizeof(MHCluster));
            array_clusters[0].id = 0;
            array_clusters[0].level = 1;
            array_clusters[0].member_count = g->count;
            array_clusters[0].member_capacity = g->count;
            array_clusters[0].members = malloc(g->count * sizeof(uint32_t));
            memcpy(array_clusters[0].members, g->monomer_ids, g->count * sizeof(uint32_t));
            array_clusters[0].representative_id = g->monomer_ids[0];
        } else {
            /* Build LSH index for this array only */
            lsh_build(ctx, g->monomer_ids, g->count);

            /* Run greedy clustering within array */
            cluster_level_greedy(ctx, threshold, g->monomer_ids, g->count,
                                 &array_clusters, &array_cluster_count);
        }

        /* Copy array clusters to global list with renumbered IDs */
        for (uint32_t c = 0; c < array_cluster_count; c++) {
            if (total_clusters >= clusters_capacity) {
                clusters_capacity *= 2;
                all_clusters = realloc(all_clusters, clusters_capacity * sizeof(MHCluster));
            }

            MHCluster *dst = &all_clusters[total_clusters];
            MHCluster *src = &array_clusters[c];

            dst->id = total_clusters;
            dst->level = 1;
            dst->member_count = src->member_count;
            dst->member_capacity = src->member_capacity;
            dst->members = src->members;  /* Transfer ownership */
            dst->representative_id = src->representative_id;
            dst->avg_similarity = 0.0f;

            total_clusters++;
        }

        /* Free array_clusters array (members transferred) */
        free(array_clusters);

        /* Progress */
        if (ctx->verbose && (a + 1) % 1000 == 0) {
            fprintf(stderr, "\r    Processed %u/%u arrays, %u clusters    ",
                    a + 1, ag.count, total_clusters);
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "\r    Processed %u arrays -> %u clusters        \n",
                ag.count, total_clusters);
    }

    array_groups_free(&ag);

    *out_clusters = all_clusters;
    *out_count = total_clusters;
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

/* Comparison function for sorting monomers by sequence_id + index */
static MHContext *sort_ctx;  /* Global for qsort comparison */

static int cmp_monomers_by_position(const void *a, const void *b) {
    uint32_t id_a = *(const uint32_t*)a;
    uint32_t id_b = *(const uint32_t*)b;
    MHMonomer *ma = &sort_ctx->monomers[id_a];
    MHMonomer *mb = &sort_ctx->monomers[id_b];

    int cmp = strcmp(ma->sequence_id, mb->sequence_id);
    if (cmp != 0) return cmp;
    return (int)ma->index - (int)mb->index;
}

/*
 * Level 1: Find subfamilies (high stringency) - Greedy Cover with Syntenic Scaffolding
 */
void find_subfamilies(MHContext *ctx) {
    if (ctx->verbose) {
        fprintf(stderr, "\n=== Level 1: Subfamilies (threshold %.2f, synteny %.2f) ===\n",
                ctx->level1_threshold, ctx->synteny_threshold);
    }

    /* Only non-flank monomers with length >= MINHASH_K participate */
    uint32_t *all_ids = malloc(ctx->monomer_count * sizeof(uint32_t));
    uint32_t non_flank_count = 0;
    uint32_t flank_count = 0;
    uint32_t short_count = 0;
    for (uint32_t i = 0; i < ctx->monomer_count; i++) {
        if (ctx->monomers[i].type != 1) {
            /* Mark flanks as unassigned (-1) */
            ctx->monomers[i].subfamily_id = UINT32_MAX;
            ctx->monomers[i].family_id = UINT32_MAX;
            ctx->monomers[i].superfamily_id = UINT32_MAX;
            flank_count++;
        } else if (ctx->monomers[i].length < MINHASH_K) {
            /* Mark short monomers as unassigned (< k, can't produce valid MinHash) */
            ctx->monomers[i].subfamily_id = UINT32_MAX;
            ctx->monomers[i].family_id = UINT32_MAX;
            ctx->monomers[i].superfamily_id = UINT32_MAX;
            short_count++;
        } else {
            /* Valid monomer for clustering */
            all_ids[non_flank_count++] = i;
        }
    }

    if (ctx->verbose) {
        fprintf(stderr, "  Valid monomers: %u (flanks: %u, short < %d bp: %u)\n",
                non_flank_count, flank_count, MINHASH_K, short_count);
    }

    /* === Syntenic Scaffolding: Pre-union adjacent monomers === */
    if (ctx->synteny_threshold > 0 && non_flank_count > 1) {
        if (ctx->verbose) {
            fprintf(stderr, "  Syntenic scaffolding: sorting by position...\n");
        }

        /* Sort by sequence_id + index */
        sort_ctx = ctx;
        qsort(all_ids, non_flank_count, sizeof(uint32_t), cmp_monomers_by_position);

        /* Initialize Union-Find for pre-union */
        uf_init(&ctx->uf, ctx->monomer_count);

        uint32_t synteny_unions = 0;
        uint32_t synteny_checked = 0;

        if (ctx->verbose) {
            fprintf(stderr, "  Pre-unioning adjacent monomers...\n");
        }

        /* Pre-union adjacent monomers in same array */
        for (uint32_t i = 0; i < non_flank_count - 1; i++) {
            uint32_t id_a = all_ids[i];
            uint32_t id_b = all_ids[i + 1];
            MHMonomer *ma = &ctx->monomers[id_a];
            MHMonomer *mb = &ctx->monomers[id_b];

            /* Check if same array and adjacent (flanks already filtered) */
            if (strcmp(ma->sequence_id, mb->sequence_id) == 0 &&
                mb->index == ma->index + 1) {

                synteny_checked++;

                /* Check similarity with low threshold */
                float sim = minhash_similarity(&ma->signature, &mb->signature);
                if (sim >= ctx->synteny_threshold) {
                    uf_union(&ctx->uf, id_a, id_b);
                    synteny_unions++;
                }
            }
        }

        if (ctx->verbose) {
            fprintf(stderr, "    Checked: %u adjacent pairs, United: %u (%.1f%%)\n",
                    synteny_checked, synteny_unions,
                    synteny_checked > 0 ? 100.0 * synteny_unions / synteny_checked : 0);
        }

        /* Now run greedy clustering, but respect pre-unions */
        /* We modify greedy to check uf_find before creating new clusters */
    }

    /* Array-local clustering: cluster each array independently */
    cluster_array_local(ctx, ctx->level1_threshold, all_ids, non_flank_count,
                        &ctx->subfamilies, &ctx->subfamily_count);

    /* === Merge clusters based on synteny pre-unions === */
    if (ctx->synteny_threshold > 0 && ctx->subfamily_count > 1) {
        if (ctx->verbose) {
            fprintf(stderr, "  Merging clusters by synteny...\n");
        }

        /* Build monomer_id -> cluster_id mapping */
        uint32_t *mono_to_cluster = malloc(ctx->monomer_count * sizeof(uint32_t));
        for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
            for (uint32_t j = 0; j < ctx->subfamilies[i].member_count; j++) {
                mono_to_cluster[ctx->subfamilies[i].members[j]] = i;
            }
        }

        /* Use Union-Find on clusters */
        uint32_t *cluster_parent = malloc(ctx->subfamily_count * sizeof(uint32_t));
        for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
            cluster_parent[i] = i;
        }

        /* Find with path compression (inline) */
        #define CFIND(x) ({ \
            uint32_t _r = (x); \
            while (cluster_parent[_r] != _r) { \
                cluster_parent[_r] = cluster_parent[cluster_parent[_r]]; \
                _r = cluster_parent[_r]; \
            } \
            _r; \
        })

        /* For each monomer, check if its synteny-root is in a different cluster */
        uint32_t synteny_merges = 0;
        for (uint32_t i = 0; i < ctx->monomer_count; i++) {
            uint32_t root = uf_find(&ctx->uf, i);
            if (root != i) {
                uint32_t cluster_i = mono_to_cluster[i];
                uint32_t cluster_root = mono_to_cluster[root];

                uint32_t ci = CFIND(cluster_i);
                uint32_t cr = CFIND(cluster_root);

                if (ci != cr) {
                    cluster_parent[cr] = ci;
                    synteny_merges++;
                }
            }
        }

        if (ctx->verbose) {
            fprintf(stderr, "    Synteny-based cluster merges: %u\n", synteny_merges);
        }

        /* Rebuild clusters if merges happened */
        if (synteny_merges > 0) {
            /* Count final clusters */
            uint32_t final_count = 0;
            uint32_t *new_id = malloc(ctx->subfamily_count * sizeof(uint32_t));
            for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
                if (CFIND(i) == i) {
                    new_id[i] = final_count++;
                }
            }
            for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
                new_id[i] = new_id[CFIND(i)];
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

            /* Move members to merged clusters */
            for (uint32_t i = 0; i < ctx->subfamily_count; i++) {
                uint32_t dest = new_id[i];
                MHCluster *src = &ctx->subfamilies[i];
                MHCluster *dst = &merged[dest];

                for (uint32_t m = 0; m < src->member_count; m++) {
                    if (dst->member_count >= dst->member_capacity) {
                        dst->member_capacity *= 2;
                        dst->members = realloc(dst->members, dst->member_capacity * sizeof(uint32_t));
                    }
                    dst->members[dst->member_count++] = src->members[m];
                }

                if (CFIND(i) == i) {
                    dst->representative_id = src->representative_id;
                }

                free(src->members);
            }
            free(ctx->subfamilies);
            free(new_id);

            ctx->subfamilies = merged;
            ctx->subfamily_count = final_count;

            if (ctx->verbose) {
                fprintf(stderr, "    After synteny merge: %u subfamilies\n", final_count);
            }
        }

        #undef CFIND
        free(cluster_parent);
        free(mono_to_cluster);
        uf_free(&ctx->uf);
    }

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
