/*
 * classify.h - Monomer classification using 8-mer bags
 *
 * Classifies monomers into hierarchical components based on
 * shared 8-mer content (bag of k-mers), then refines with Jaccard similarity.
 */

#ifndef CLASSIFY_H
#define CLASSIFY_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <pthread.h>

/* Constants */
#define KMER_SIZE 8
#define KMER_BUCKETS 65536      /* 4^8 possible 8-mers */
#define BITVEC_WORDS (KMER_BUCKETS / 64)  /* 1024 uint64_t words for bitvector */
#define DEFAULT_MIN_SHARED_KMERS 3        /* Min shared k-mers for component */
#define DEFAULT_JACCARD_THRESHOLD 0.8     /* Jaccard for subcomponents */
#define SHORT_MONOMER_THRESHOLD 32

/* Monomer types */
typedef enum {
    TYPE_LEFT_FLANK = 0,
    TYPE_MONOMER = 1,
    TYPE_RIGHT_FLANK = 2
} MonomerType;

/* 8-mer occurrence in inverted index */
typedef struct {
    uint32_t monomer_id;
} KmerHit;

/* Bucket in inverted index */
typedef struct {
    KmerHit *hits;
    uint32_t count;
    uint32_t capacity;
} KmerBucket;

/* Bitvector for bag of k-mers (65536 bits = 1024 uint64_t) */
typedef struct {
    uint64_t bits[BITVEC_WORDS];
    uint32_t popcount;  /* Number of set bits (distinct k-mers) */
} KmerBitvec;

/* Monomer record */
typedef struct {
    char *sequence_id;      /* Original array ID */
    char *orientation;      /* fwd or rev */
    char *sequence;         /* DNA sequence */
    uint32_t length;        /* Sequence length */
    uint32_t index;         /* Position in original array */
    MonomerType type;       /* LEFT_FLANK, MONOMER, RIGHT_FLANK */
    bool is_flank;

    /* Classification results */
    uint32_t component_id;      /* Connected component */
    uint32_t subcomponent_id;   /* Subcomponent within component */
    float jaccard_to_rep;       /* Jaccard to subcomponent representative */

    /* Bag of k-mers representation */
    KmerBitvec bitvec;
} Monomer;

/* Component (connected by shared k-mers) */
typedef struct {
    uint32_t *member_ids;
    uint32_t member_count;
    uint32_t member_capacity;
    uint32_t representative_id;
} Component;

/* Subcomponent (high Jaccard within component) */
typedef struct {
    uint32_t component_id;
    uint32_t *member_ids;
    uint32_t member_count;
    uint32_t member_capacity;
    uint32_t representative_id;
    float avg_jaccard;
} Subcomponent;

/* Global classification context */
typedef struct {
    /* Input data */
    Monomer *monomers;
    uint32_t monomer_count;
    uint32_t monomer_capacity;

    /* 8-mer inverted index */
    KmerBucket buckets[KMER_BUCKETS];

    /* Components */
    Component *components;
    uint32_t component_count;
    uint32_t component_capacity;

    /* Subcomponents */
    Subcomponent *subcomponents;
    uint32_t subcomponent_count;
    uint32_t subcomponent_capacity;

    /* Parameters */
    int min_shared_kmers;       /* Min shared k-mers for connectivity */
    float jaccard_threshold;    /* Jaccard threshold for subcomponents */
    bool skip_flanks;
    int num_threads;
    bool verbose;
} ClassifyContext;

/* ==================== Bitvector operations ==================== */

/* Set bit for k-mer */
static inline void bitvec_set(KmerBitvec *bv, uint16_t kmer) {
    bv->bits[kmer / 64] |= (1ULL << (kmer % 64));
}

/* Test bit for k-mer */
static inline bool bitvec_test(const KmerBitvec *bv, uint16_t kmer) {
    return (bv->bits[kmer / 64] & (1ULL << (kmer % 64))) != 0;
}

/* Count set bits (popcount) */
static inline uint32_t bitvec_popcount(const KmerBitvec *bv) {
    uint32_t count = 0;
    for (int i = 0; i < BITVEC_WORDS; i++) {
        count += __builtin_popcountll(bv->bits[i]);
    }
    return count;
}

/* Count intersection (AND) */
static inline uint32_t bitvec_intersection_count(const KmerBitvec *a, const KmerBitvec *b) {
    uint32_t count = 0;
    for (int i = 0; i < BITVEC_WORDS; i++) {
        count += __builtin_popcountll(a->bits[i] & b->bits[i]);
    }
    return count;
}

/* Count union (OR) */
static inline uint32_t bitvec_union_count(const KmerBitvec *a, const KmerBitvec *b) {
    uint32_t count = 0;
    for (int i = 0; i < BITVEC_WORDS; i++) {
        count += __builtin_popcountll(a->bits[i] | b->bits[i]);
    }
    return count;
}

/* Jaccard similarity */
static inline float bitvec_jaccard(const KmerBitvec *a, const KmerBitvec *b) {
    uint32_t intersection = bitvec_intersection_count(a, b);
    uint32_t union_size = bitvec_union_count(a, b);
    if (union_size == 0) return 0.0f;
    return (float)intersection / (float)union_size;
}

/* Clear bitvector */
static inline void bitvec_clear(KmerBitvec *bv) {
    for (int i = 0; i < BITVEC_WORDS; i++) {
        bv->bits[i] = 0;
    }
    bv->popcount = 0;
}

/* ==================== K-mer encoding ==================== */

/* Encode nucleotide to 2-bit value */
static inline uint8_t nuc_to_bits(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;  /* N treated as A */
    }
}

/* Compute reverse complement of encoded 8-mer
 * Complement: A(00)↔T(11), C(01)↔G(10) => XOR with 0xFFFF
 * Then reverse order of 8 2-bit pairs */
static inline uint16_t revcomp_8mer(uint16_t kmer) {
    uint16_t comp = kmer ^ 0xFFFF;

    /* Reverse 8 2-bit pairs */
    uint16_t rev = 0;
    rev |= (comp & 0x0003) << 14;  /* pos 0 -> pos 7 */
    rev |= (comp & 0x000C) << 10;  /* pos 1 -> pos 6 */
    rev |= (comp & 0x0030) << 6;   /* pos 2 -> pos 5 */
    rev |= (comp & 0x00C0) << 2;   /* pos 3 -> pos 4 */
    rev |= (comp & 0x0300) >> 2;   /* pos 4 -> pos 3 */
    rev |= (comp & 0x0C00) >> 6;   /* pos 5 -> pos 2 */
    rev |= (comp & 0x3000) >> 10;  /* pos 6 -> pos 1 */
    rev |= (comp & 0xC000) >> 14;  /* pos 7 -> pos 0 */

    return rev;
}

/* Encode 8-mer to canonical 16-bit integer (min of forward and revcomp) */
static inline uint16_t encode_8mer(const char *seq) {
    uint16_t val = 0;
    for (int i = 0; i < 8; i++) {
        val = (val << 2) | nuc_to_bits(seq[i]);
    }
    uint16_t rc = revcomp_8mer(val);
    return val < rc ? val : rc;  /* Canonical = lexicographically smaller */
}

/* ==================== classify_kmer.c ==================== */

/* Decode 16-bit integer to 8-mer string */
void decode_8mer(uint16_t val, char *out);

/* Build bitvector for a monomer (all 8-mers, step=1) */
void build_bitvector(Monomer *mono);

/* Build inverted index from all monomers */
void build_kmer_index(ClassifyContext *ctx);

/* ==================== classify_cluster.c ==================== */

/* Find connected components (monomers sharing >= min_shared_kmers) */
void find_components(ClassifyContext *ctx);

/* Find subcomponents within each component (Jaccard >= threshold) */
void find_subcomponents(ClassifyContext *ctx);

/* Main clustering entry point */
void cluster_monomers(ClassifyContext *ctx);

/* ==================== classify_main.c ==================== */

/* Initialize classification context */
ClassifyContext* classify_init(void);

/* Free classification context */
void classify_free(ClassifyContext *ctx);

/* Load monomers from TSV file */
int load_monomers_tsv(ClassifyContext *ctx, const char *filename);

/* Main classification pipeline */
int classify_monomers(ClassifyContext *ctx);

/* Output functions */
int write_components_tsv(ClassifyContext *ctx, const char *filename);
int write_component_stats(ClassifyContext *ctx, const char *filename);

#endif /* CLASSIFY_H */
