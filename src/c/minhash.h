/*
 * minhash.h - MinHash-based hierarchical classification
 *
 * Cascade Clustering approach:
 * Level 1 (Subfamilies): threshold 0.75-0.80, high stringency
 * Level 2 (Families): threshold 0.4-0.5, medium stringency
 * Level 3 (Superfamilies): threshold 0.2, low stringency
 */

#ifndef MINHASH_H
#define MINHASH_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <pthread.h>

/* ==================== Constants ==================== */

#define MINHASH_K 12                /* K-mer size for MinHash */
#define MINHASH_SIZE 128            /* Number of hash functions */
#define LSH_BANDS 16                /* Number of LSH bands */
#define LSH_ROWS 8                  /* Rows per band (BANDS * ROWS = MINHASH_SIZE) */

/* Thresholds for cascade levels */
#define LEVEL1_THRESHOLD 0.80       /* Subfamilies */
#define LEVEL2_THRESHOLD 0.50       /* Families */
#define LEVEL3_THRESHOLD 0.20       /* Superfamilies */

/* Hash table sizes (prime numbers for better distribution) */
#define LSH_BUCKETS 100003

/* ==================== Data Structures ==================== */

/* MinHash signature for a monomer */
typedef struct {
    uint64_t hashes[MINHASH_SIZE];
} MinHashSig;

/* Monomer with MinHash signature */
typedef struct {
    uint32_t id;                    /* Original monomer ID */
    char *sequence;                 /* DNA sequence */
    uint32_t length;                /* Sequence length */
    MinHashSig signature;           /* MinHash signature */

    /* Hierarchy assignments */
    uint32_t subfamily_id;          /* Level 1 cluster */
    uint32_t family_id;             /* Level 2 cluster */
    uint32_t superfamily_id;        /* Level 3 cluster */

    /* Metadata from input */
    char *sequence_id;              /* Original array ID */
    char *orientation;              /* fwd/rev */
    uint32_t index;                 /* Position in array */
    int type;                       /* 0=left_flank, 1=monomer, 2=right_flank */
    bool is_flank;
} MHMonomer;

/* Cluster at any level */
typedef struct {
    uint32_t id;                    /* Cluster ID at this level */
    uint32_t *members;              /* Array of member IDs */
    uint32_t member_count;
    uint32_t member_capacity;
    uint32_t representative_id;     /* Centroid monomer ID */
    uint32_t level;                 /* 1=subfamily, 2=family, 3=superfamily */
    float avg_similarity;           /* Average internal similarity */
} MHCluster;

/* LSH bucket entry */
typedef struct LSHEntry {
    uint32_t monomer_id;
    struct LSHEntry *next;
} LSHEntry;

/* LSH index for one band */
typedef struct {
    LSHEntry *buckets[LSH_BUCKETS];
} LSHBand;

/* Full LSH index (all bands) */
typedef struct {
    LSHBand bands[LSH_BANDS];
} LSHIndex;

/* Union-Find for clustering */
typedef struct {
    uint32_t *parent;
    uint32_t *rank;
    uint32_t size;
} UnionFind;

/* Main context for MinHash classification */
typedef struct {
    /* Input data */
    MHMonomer *monomers;
    uint32_t monomer_count;
    uint32_t monomer_capacity;

    /* LSH index */
    LSHIndex lsh;

    /* Union-Find for current clustering level */
    UnionFind uf;
    pthread_mutex_t uf_mutex;       /* Mutex for Union-Find operations */

    /* Clusters at each level */
    MHCluster *subfamilies;         /* Level 1 */
    uint32_t subfamily_count;
    uint32_t subfamily_capacity;

    MHCluster *families;            /* Level 2 */
    uint32_t family_count;
    uint32_t family_capacity;

    MHCluster *superfamilies;       /* Level 3 */
    uint32_t superfamily_count;
    uint32_t superfamily_capacity;

    /* Parameters */
    float level1_threshold;
    float level2_threshold;
    float level3_threshold;
    int num_threads;
    bool verbose;
} MHContext;

/* ==================== Hash Functions ==================== */

/* MurmurHash3 finalizer for 64-bit */
static inline uint64_t murmur_mix(uint64_t h) {
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return h;
}

/* Hash a k-mer with seed */
static inline uint64_t hash_kmer(uint64_t kmer, uint32_t seed) {
    return murmur_mix(kmer ^ ((uint64_t)seed * 0x9e3779b97f4a7c15ULL));
}

/* ==================== K-mer Encoding ==================== */

/* Nucleotide to 2-bit encoding */
static inline uint8_t nuc_to_2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;  /* N treated as A */
    }
}

/* Complement of 2-bit nucleotide */
static inline uint8_t complement_2bit(uint8_t b) {
    return b ^ 3;  /* A(0)↔T(3), C(1)↔G(2) */
}

/* Encode k-mer (k=12) to 24-bit value stored in uint64_t */
static inline uint64_t encode_kmer(const char *seq, int k) {
    uint64_t val = 0;
    for (int i = 0; i < k; i++) {
        val = (val << 2) | nuc_to_2bit(seq[i]);
    }
    return val;
}

/* Reverse complement of encoded k-mer */
static inline uint64_t revcomp_kmer(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        rc = (rc << 2) | complement_2bit(kmer & 3);
        kmer >>= 2;
    }
    return rc;
}

/* Canonical k-mer (min of forward and revcomp) */
static inline uint64_t canonical_kmer(uint64_t kmer, int k) {
    uint64_t rc = revcomp_kmer(kmer, k);
    return kmer < rc ? kmer : rc;
}

/* ==================== MinHash Similarity ==================== */

/* Estimate Jaccard similarity from MinHash signatures */
static inline float minhash_similarity(const MinHashSig *a, const MinHashSig *b) {
    uint32_t matches = 0;
    for (int i = 0; i < MINHASH_SIZE; i++) {
        if (a->hashes[i] == b->hashes[i]) {
            matches++;
        }
    }
    return (float)matches / MINHASH_SIZE;
}

/* ==================== Union-Find Operations ==================== */

static inline void uf_init(UnionFind *uf, uint32_t size) {
    uf->parent = malloc(size * sizeof(uint32_t));
    uf->rank = calloc(size, sizeof(uint32_t));
    uf->size = size;
    for (uint32_t i = 0; i < size; i++) {
        uf->parent[i] = i;
    }
}

static inline uint32_t uf_find(UnionFind *uf, uint32_t x) {
    if (uf->parent[x] != x) {
        uf->parent[x] = uf_find(uf, uf->parent[x]);  /* Path compression */
    }
    return uf->parent[x];
}

static inline void uf_union(UnionFind *uf, uint32_t x, uint32_t y) {
    uint32_t rx = uf_find(uf, x);
    uint32_t ry = uf_find(uf, y);
    if (rx == ry) return;

    /* Union by rank */
    if (uf->rank[rx] < uf->rank[ry]) {
        uf->parent[rx] = ry;
    } else if (uf->rank[rx] > uf->rank[ry]) {
        uf->parent[ry] = rx;
    } else {
        uf->parent[ry] = rx;
        uf->rank[rx]++;
    }
}

static inline void uf_free(UnionFind *uf) {
    free(uf->parent);
    free(uf->rank);
}

/* ==================== Function Declarations ==================== */

/* minhash_sig.c - Signature computation */
void compute_minhash(MHMonomer *mono);
void compute_all_signatures(MHContext *ctx);

/* minhash_lsh.c - LSH indexing */
void lsh_init(LSHIndex *lsh);
void lsh_free(LSHIndex *lsh);
void lsh_insert(LSHIndex *lsh, uint32_t monomer_id, const MinHashSig *sig);
void lsh_clear(LSHIndex *lsh);
uint32_t* lsh_query_candidates(LSHIndex *lsh, const MinHashSig *sig, uint32_t *count);

/* minhash_cluster.c - Cascade clustering */
void cluster_level(MHContext *ctx, float threshold, uint32_t *ids, uint32_t count);
void find_subfamilies(MHContext *ctx);
void find_families(MHContext *ctx);
void find_superfamilies(MHContext *ctx);
uint32_t select_centroid(MHContext *ctx, uint32_t *members, uint32_t count);

/* minhash_main.c - I/O and main pipeline */
MHContext* mh_init(void);
void mh_free(MHContext *ctx);
int mh_load_monomers(MHContext *ctx, const char *filename);
int mh_classify(MHContext *ctx);
int mh_write_results(MHContext *ctx, const char *prefix);

#endif /* MINHASH_H */
