/*
 * K-mer Enrichment Analysis for Centromeric Regions (v2 with TF-DF)
 *
 * Compares k-mer frequencies between CDR (deep centromere) and flanking regions.
 * Now tracks both TF (total count) and DF (number of unique regions).
 * Multi-threaded implementation for efficiency.
 *
 * Usage: kmer_enrichment -g genome.fa -c cdr.bed -f centromere.bed -o output.tsv [-k max_k] [-m min_count] [-t threads]
 *
 * Author: ArraySplitter project
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <pthread.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>

#define MAX_LINE 65536
#define MAX_KMER 20
#define MAX_REGIONS 100000
#define INITIAL_SEQ_SIZE 1000000000  // 1GB initial buffer

// Region structure
typedef struct {
    char chrom[64];
    uint64_t start;
    uint64_t end;
    int region_id;  // unique ID for DF tracking
} Region;

// Chromosome sequence
typedef struct {
    char name[64];
    char *seq;
    uint64_t length;
    uint64_t offset;  // offset in file for fai
} ChromSeq;

// K-mer count entry with TF and DF
typedef struct {
    uint64_t kmer;           // encoded k-mer
    uint64_t cdr_count;      // TF in CDR
    uint64_t flank_count;    // TF in flanking
    uint32_t cdr_regions;    // DF - number of unique CDR regions
    uint32_t flank_regions;  // DF - number of unique flanking regions
    int k;
    // Bitmap for tracking which regions contain this k-mer (for small region counts)
    uint64_t *cdr_bitmap;    // bitmap of CDR regions (allocated if needed)
    uint64_t *flank_bitmap;  // bitmap of flanking regions
} KmerCount;

// Global data
ChromSeq *chromosomes = NULL;
int num_chromosomes = 0;

Region *cdr_regions = NULL;
int num_cdr_regions = 0;

Region *centromere_regions = NULL;
int num_centromere_regions = 0;

// Flanking regions (computed from centromere - CDR)
Region *flank_regions = NULL;
int num_flank_regions = 0;

// K-mer hash tables (one per k value)
typedef struct {
    KmerCount *entries;
    uint64_t capacity;
    uint64_t count;
} KmerTable;

KmerTable kmer_tables[MAX_KMER + 1];
pthread_mutex_t table_mutexes[MAX_KMER + 1];

// Bitmap size (in 64-bit words) for region tracking
int cdr_bitmap_words = 0;
int flank_bitmap_words = 0;

// Parameters
int max_k = 20;
int min_count = 40;
int min_df = 1;  // minimum DF to report
int num_threads = 8;
char *genome_file = NULL;
char *cdr_bed = NULL;
char *centromere_bed = NULL;
char *output_file = NULL;

// Convert nucleotide to 2-bit encoding
static inline int base_to_code(char c) {
    switch (toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;  // N or other
    }
}

// Convert 2-bit code to nucleotide
static inline char code_to_base(int code) {
    const char bases[] = "ACGT";
    return bases[code & 3];
}

// Decode k-mer to string
void decode_kmer(uint64_t kmer, int k, char *out) {
    for (int i = k - 1; i >= 0; i--) {
        out[i] = code_to_base(kmer & 3);
        kmer >>= 2;
    }
    out[k] = '\0';
}

// Reverse complement of encoded k-mer
static inline uint64_t reverse_complement(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        rc = (rc << 2) | (3 - (kmer & 3));
        kmer >>= 2;
    }
    return rc;
}

// Canonical k-mer (smaller of kmer and its reverse complement)
static inline uint64_t canonical_kmer(uint64_t kmer, int k) {
    uint64_t rc = reverse_complement(kmer, k);
    return (kmer < rc) ? kmer : rc;
}

// Hash function for k-mer
static inline uint64_t hash_kmer(uint64_t kmer, uint64_t capacity) {
    kmer ^= kmer >> 33;
    kmer *= 0xff51afd7ed558ccdULL;
    kmer ^= kmer >> 33;
    kmer *= 0xc4ceb9fe1a85ec53ULL;
    kmer ^= kmer >> 33;
    return kmer % capacity;
}

// Initialize k-mer table
void init_kmer_table(int k) {
    uint64_t capacity;
    if (k <= 14) {
        capacity = 1ULL << (2 * k);
    } else {
        capacity = 1ULL << 28;  // ~256M entries for large k
    }

    kmer_tables[k].entries = calloc(capacity, sizeof(KmerCount));
    kmer_tables[k].capacity = capacity;
    kmer_tables[k].count = 0;
    pthread_mutex_init(&table_mutexes[k], NULL);
}

// Set bit in bitmap
static inline void set_bit(uint64_t *bitmap, int bit) {
    if (bitmap) {
        bitmap[bit / 64] |= (1ULL << (bit % 64));
    }
}

// Check bit in bitmap
static inline int get_bit(uint64_t *bitmap, int bit) {
    if (!bitmap) return 0;
    return (bitmap[bit / 64] >> (bit % 64)) & 1;
}

// Add or update k-mer count with region tracking
void add_kmer(int k, uint64_t kmer, int is_cdr, int region_id) {
    KmerTable *table = &kmer_tables[k];
    uint64_t idx;

    if (k <= 14) {
        idx = kmer;
    } else {
        idx = hash_kmer(kmer, table->capacity);
        while (table->entries[idx].k != 0 && table->entries[idx].kmer != kmer) {
            idx = (idx + 1) % table->capacity;
        }
    }

    pthread_mutex_lock(&table_mutexes[k]);

    KmerCount *entry = &table->entries[idx];

    if (entry->k == 0) {
        entry->kmer = kmer;
        entry->k = k;
        entry->cdr_count = 0;
        entry->flank_count = 0;
        entry->cdr_regions = 0;
        entry->flank_regions = 0;

        // Allocate bitmaps for region tracking
        if (cdr_bitmap_words > 0) {
            entry->cdr_bitmap = calloc(cdr_bitmap_words, sizeof(uint64_t));
        }
        if (flank_bitmap_words > 0) {
            entry->flank_bitmap = calloc(flank_bitmap_words, sizeof(uint64_t));
        }
        table->count++;
    }

    if (is_cdr) {
        entry->cdr_count++;
        // Update DF if this is a new region for this k-mer
        if (entry->cdr_bitmap && !get_bit(entry->cdr_bitmap, region_id)) {
            set_bit(entry->cdr_bitmap, region_id);
            entry->cdr_regions++;
        }
    } else {
        entry->flank_count++;
        if (entry->flank_bitmap && !get_bit(entry->flank_bitmap, region_id)) {
            set_bit(entry->flank_bitmap, region_id);
            entry->flank_regions++;
        }
    }

    pthread_mutex_unlock(&table_mutexes[k]);
}

// Read full genome FASTA
int read_genome(const char *fasta_path) {
    FILE *fp = fopen(fasta_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open FASTA file: %s\n", fasta_path);
        return -1;
    }

    chromosomes = malloc(sizeof(ChromSeq) * 1000);
    char line[MAX_LINE];
    int current_chr = -1;
    size_t seq_pos = 0;
    size_t seq_capacity = 0;

    fprintf(stderr, "Reading genome...\n");

    while (fgets(line, MAX_LINE, fp)) {
        if (line[0] == '>') {
            if (current_chr >= 0) {
                chromosomes[current_chr].length = seq_pos;
            }

            current_chr = num_chromosomes++;
            ChromSeq *chr = &chromosomes[current_chr];

            char *name_end = strchr(line + 1, ' ');
            if (!name_end) name_end = strchr(line + 1, '\n');
            if (!name_end) name_end = strchr(line + 1, '\t');

            int name_len = name_end ? (int)(name_end - line - 1) : (int)strlen(line + 1);
            if (name_len > 63) name_len = 63;
            strncpy(chr->name, line + 1, name_len);
            chr->name[name_len] = '\0';

            while (name_len > 0 && (chr->name[name_len-1] == '\n' || chr->name[name_len-1] == '\r')) {
                chr->name[--name_len] = '\0';
            }

            seq_capacity = 300000000;
            chr->seq = malloc(seq_capacity);
            seq_pos = 0;

            fprintf(stderr, "  Reading %s...\n", chr->name);
        } else if (current_chr >= 0) {
            ChromSeq *chr = &chromosomes[current_chr];
            size_t len = strlen(line);

            while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
                line[--len] = '\0';
            }

            if (seq_pos + len >= seq_capacity) {
                seq_capacity *= 2;
                chr->seq = realloc(chr->seq, seq_capacity);
            }

            for (size_t i = 0; i < len; i++) {
                chr->seq[seq_pos++] = toupper(line[i]);
            }
        }
    }

    if (current_chr >= 0) {
        chromosomes[current_chr].length = seq_pos;
    }

    fclose(fp);
    fprintf(stderr, "Loaded %d chromosomes\n", num_chromosomes);
    return 0;
}

// Find chromosome by name
ChromSeq* find_chromosome(const char *name) {
    for (int i = 0; i < num_chromosomes; i++) {
        if (strcmp(chromosomes[i].name, name) == 0) {
            return &chromosomes[i];
        }
    }
    return NULL;
}

// Read BED file
int read_bed(const char *bed_path, Region **regions, int *num_regions) {
    FILE *fp = fopen(bed_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open BED file: %s\n", bed_path);
        return -1;
    }

    *regions = malloc(sizeof(Region) * MAX_REGIONS);
    *num_regions = 0;
    char line[MAX_LINE];

    while (fgets(line, MAX_LINE, fp) && *num_regions < MAX_REGIONS) {
        if (line[0] == '#' || line[0] == '\n') continue;

        Region *r = &(*regions)[*num_regions];
        if (sscanf(line, "%63s\t%llu\t%llu", r->chrom,
                   (unsigned long long*)&r->start,
                   (unsigned long long*)&r->end) >= 3) {
            r->region_id = *num_regions;
            (*num_regions)++;
        }
    }

    fclose(fp);
    fprintf(stderr, "Loaded %d regions from %s\n", *num_regions, bed_path);
    return 0;
}

// Compute flanking regions (centromere - CDR)
void compute_flanking_regions() {
    flank_regions = malloc(sizeof(Region) * MAX_REGIONS * 2);
    num_flank_regions = 0;

    for (int i = 0; i < num_centromere_regions; i++) {
        Region *centro = &centromere_regions[i];
        uint64_t current_pos = centro->start;

        // Find overlapping CDR regions
        for (int j = 0; j < num_cdr_regions; j++) {
            Region *cdr = &cdr_regions[j];

            if (strcmp(cdr->chrom, centro->chrom) != 0) continue;
            if (cdr->end <= centro->start || cdr->start >= centro->end) continue;

            uint64_t cdr_start = cdr->start < centro->start ? centro->start : cdr->start;
            uint64_t cdr_end = cdr->end > centro->end ? centro->end : cdr->end;

            // Add flanking region before CDR
            if (current_pos < cdr_start) {
                Region *flank = &flank_regions[num_flank_regions];
                strcpy(flank->chrom, centro->chrom);
                flank->start = current_pos;
                flank->end = cdr_start;
                flank->region_id = num_flank_regions;
                num_flank_regions++;
            }
            current_pos = cdr_end;
        }

        // Add remaining flanking after last CDR
        if (current_pos < centro->end) {
            Region *flank = &flank_regions[num_flank_regions];
            strcpy(flank->chrom, centro->chrom);
            flank->start = current_pos;
            flank->end = centro->end;
            flank->region_id = num_flank_regions;
            num_flank_regions++;
        }
    }

    fprintf(stderr, "Computed %d flanking regions\n", num_flank_regions);
}

// Process sequence for k-mers with region tracking
void process_sequence_with_region(const char *seq, uint64_t len, int is_cdr, int region_id, int min_k, int max_k_local) {
    if (len < (uint64_t)min_k) return;

    for (int k = min_k; k <= max_k_local && k <= (int)len; k++) {
        uint64_t valid_pos = 0;
        uint64_t current_kmer = 0;
        uint64_t mask = (1ULL << (2*k)) - 1;

        for (uint64_t i = 0; i < len; i++) {
            int code = base_to_code(seq[i]);

            if (code < 0) {
                valid_pos = 0;
                current_kmer = 0;
                continue;
            }

            current_kmer = ((current_kmer << 2) | code) & mask;
            valid_pos++;

            if (valid_pos >= (uint64_t)k) {
                uint64_t canon = canonical_kmer(current_kmer, k);
                add_kmer(k, canon, is_cdr, region_id);
            }
        }
    }
}

// Thread argument
typedef struct {
    int thread_id;
    int num_threads;
    int is_cdr;
    int min_k;
    int max_k;
} ThreadArg;

// Thread worker for CDR regions
void* count_kmers_cdr_thread(void *arg) {
    ThreadArg *targ = (ThreadArg*)arg;

    for (int i = targ->thread_id; i < num_cdr_regions; i += targ->num_threads) {
        Region *r = &cdr_regions[i];
        ChromSeq *chr = find_chromosome(r->chrom);

        if (!chr || !chr->seq) continue;

        uint64_t start = r->start;
        uint64_t end = r->end;
        if (end > chr->length) end = chr->length;
        if (start >= end) continue;

        process_sequence_with_region(chr->seq + start, end - start, 1, r->region_id,
                                      targ->min_k, targ->max_k);
    }

    return NULL;
}

// Thread worker for flanking regions
void* count_kmers_flank_thread(void *arg) {
    ThreadArg *targ = (ThreadArg*)arg;

    for (int i = targ->thread_id; i < num_flank_regions; i += targ->num_threads) {
        Region *r = &flank_regions[i];
        ChromSeq *chr = find_chromosome(r->chrom);

        if (!chr || !chr->seq) continue;

        uint64_t start = r->start;
        uint64_t end = r->end;
        if (end > chr->length) end = chr->length;
        if (start >= end) continue;

        process_sequence_with_region(chr->seq + start, end - start, 0, r->region_id,
                                      targ->min_k, targ->max_k);
    }

    return NULL;
}

// Comparison function for sorting by CDR DF (most regions first), then by enrichment
int compare_cdr_df(const void *a, const void *b) {
    const KmerCount *ka = (const KmerCount*)a;
    const KmerCount *kb = (const KmerCount*)b;

    // First sort by CDR DF (descending)
    if (kb->cdr_regions != ka->cdr_regions) {
        return (int)kb->cdr_regions - (int)ka->cdr_regions;
    }

    // Then by log2 fold change
    double fc_a = (ka->cdr_count + 1.0) / (ka->flank_count + 1.0);
    double fc_b = (kb->cdr_count + 1.0) / (kb->flank_count + 1.0);

    if (fc_b > fc_a) return 1;
    if (fc_b < fc_a) return -1;
    return 0;
}

// Comparison for flanking DF
int compare_flank_df(const void *a, const void *b) {
    const KmerCount *ka = (const KmerCount*)a;
    const KmerCount *kb = (const KmerCount*)b;

    if (kb->flank_regions != ka->flank_regions) {
        return (int)kb->flank_regions - (int)ka->flank_regions;
    }

    double fc_a = (ka->flank_count + 1.0) / (ka->cdr_count + 1.0);
    double fc_b = (kb->flank_count + 1.0) / (kb->cdr_count + 1.0);

    if (fc_b > fc_a) return 1;
    if (fc_b < fc_a) return -1;
    return 0;
}

// Write results with TF and DF
void write_results(const char *output_path) {
    FILE *fp = fopen(output_path, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_path);
        return;
    }

    // Header with TF and DF columns
    fprintf(fp, "k\trank\tkmer\tcdr_tf\tcdr_df\tflank_tf\tflank_df\ttotal_tf\tlog2_fc\tcdr_df_pct\tflank_df_pct\tenrichment\n");

    for (int k = 1; k <= max_k; k++) {
        KmerTable *table = &kmer_tables[k];

        // Collect significant k-mers (by total count)
        KmerCount *significant = malloc(sizeof(KmerCount) * table->count);
        int num_significant = 0;

        for (uint64_t i = 0; i < table->capacity; i++) {
            if (table->entries[i].k == k) {
                uint64_t total = table->entries[i].cdr_count + table->entries[i].flank_count;
                if (total >= (uint64_t)min_count) {
                    significant[num_significant++] = table->entries[i];
                }
            }
        }

        // Sort by CDR DF and output top CDR-enriched
        qsort(significant, num_significant, sizeof(KmerCount), compare_cdr_df);

        char kmer_str[MAX_KMER + 1];
        int top_n = num_significant < 20 ? num_significant : 20;

        for (int i = 0; i < top_n; i++) {
            KmerCount *kc = &significant[i];

            // Skip if CDR DF is 0
            if (kc->cdr_regions == 0) continue;

            decode_kmer(kc->kmer, k, kmer_str);

            uint64_t total = kc->cdr_count + kc->flank_count;
            double fc = (kc->cdr_count + 1.0) / (kc->flank_count + 1.0);
            double log2_fc = log2(fc);
            double cdr_df_pct = 100.0 * kc->cdr_regions / num_cdr_regions;
            double flank_df_pct = num_flank_regions > 0 ? 100.0 * kc->flank_regions / num_flank_regions : 0;
            const char *enrichment = (log2_fc > 0.5) ? "CDR" : (log2_fc < -0.5) ? "Flanking" : "Neutral";

            fprintf(fp, "%d\t%d\t%s\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%.3f\t%.1f\t%.1f\t%s\n",
                    k, i + 1, kmer_str,
                    kc->cdr_count, kc->cdr_regions,
                    kc->flank_count, kc->flank_regions,
                    total, log2_fc, cdr_df_pct, flank_df_pct, enrichment);
        }

        // Sort by Flanking DF and output top flanking-enriched
        qsort(significant, num_significant, sizeof(KmerCount), compare_flank_df);

        for (int i = 0; i < top_n; i++) {
            KmerCount *kc = &significant[i];

            if (kc->flank_regions == 0) continue;

            // Skip if already CDR-enriched
            double fc = (kc->cdr_count + 1.0) / (kc->flank_count + 1.0);
            if (fc > 1.0) continue;

            decode_kmer(kc->kmer, k, kmer_str);

            uint64_t total = kc->cdr_count + kc->flank_count;
            double log2_fc = log2(fc);
            double cdr_df_pct = 100.0 * kc->cdr_regions / num_cdr_regions;
            double flank_df_pct = num_flank_regions > 0 ? 100.0 * kc->flank_regions / num_flank_regions : 0;

            fprintf(fp, "%d\t%d\t%s\t%" PRIu64 "\t%u\t%" PRIu64 "\t%u\t%" PRIu64 "\t%.3f\t%.1f\t%.1f\tFlanking\n",
                    k, -(i + 1), kmer_str,
                    kc->cdr_count, kc->cdr_regions,
                    kc->flank_count, kc->flank_regions,
                    total, log2_fc, cdr_df_pct, flank_df_pct);
        }

        free(significant);
        fprintf(stderr, "k=%d: %d significant k-mers (CDR regions: %d, Flank regions: %d)\n",
                k, num_significant, num_cdr_regions, num_flank_regions);
    }

    fclose(fp);
    fprintf(stderr, "Results written to %s\n", output_path);
}

void print_usage(const char *prog) {
    fprintf(stderr, "K-mer Enrichment Analysis with TF-DF (v2)\n\n");
    fprintf(stderr, "Usage: %s [options]\n\n", prog);
    fprintf(stderr, "Required:\n");
    fprintf(stderr, "  -g, --genome FILE     Genome FASTA file\n");
    fprintf(stderr, "  -c, --cdr FILE        CDR regions BED file\n");
    fprintf(stderr, "  -f, --flanking FILE   Centromere regions BED file\n");
    fprintf(stderr, "  -o, --output FILE     Output TSV file\n\n");
    fprintf(stderr, "Optional:\n");
    fprintf(stderr, "  -k, --max-k INT       Maximum k-mer size [default: 20]\n");
    fprintf(stderr, "  -m, --min-count INT   Minimum total count threshold [default: 40]\n");
    fprintf(stderr, "  -t, --threads INT     Number of threads [default: 8]\n");
    fprintf(stderr, "  -h, --help            Show this help\n\n");
    fprintf(stderr, "Output columns:\n");
    fprintf(stderr, "  k          - k-mer size\n");
    fprintf(stderr, "  rank       - rank (positive=CDR-enriched, negative=Flanking-enriched)\n");
    fprintf(stderr, "  kmer       - k-mer sequence\n");
    fprintf(stderr, "  cdr_tf     - total count in CDR regions\n");
    fprintf(stderr, "  cdr_df     - number of unique CDR regions containing k-mer\n");
    fprintf(stderr, "  flank_tf   - total count in flanking regions\n");
    fprintf(stderr, "  flank_df   - number of unique flanking regions containing k-mer\n");
    fprintf(stderr, "  total_tf   - total count\n");
    fprintf(stderr, "  log2_fc    - log2 fold change (CDR/Flanking)\n");
    fprintf(stderr, "  cdr_df_pct - %% of CDR regions containing k-mer\n");
    fprintf(stderr, "  flank_df_pct - %% of flanking regions containing k-mer\n");
}

int main(int argc, char *argv[]) {
    static struct option long_options[] = {
        {"genome",    required_argument, 0, 'g'},
        {"cdr",       required_argument, 0, 'c'},
        {"flanking",  required_argument, 0, 'f'},
        {"output",    required_argument, 0, 'o'},
        {"max-k",     required_argument, 0, 'k'},
        {"min-count", required_argument, 0, 'm'},
        {"threads",   required_argument, 0, 't'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "g:c:f:o:k:m:t:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'g': genome_file = optarg; break;
            case 'c': cdr_bed = optarg; break;
            case 'f': centromere_bed = optarg; break;
            case 'o': output_file = optarg; break;
            case 'k': max_k = atoi(optarg); break;
            case 'm': min_count = atoi(optarg); break;
            case 't': num_threads = atoi(optarg); break;
            case 'h': print_usage(argv[0]); return 0;
            default:  print_usage(argv[0]); return 1;
        }
    }

    if (!genome_file || !cdr_bed || !centromere_bed || !output_file) {
        fprintf(stderr, "Error: Missing required arguments\n\n");
        print_usage(argv[0]);
        return 1;
    }

    if (max_k < 1 || max_k > MAX_KMER) {
        fprintf(stderr, "Error: max-k must be between 1 and %d\n", MAX_KMER);
        return 1;
    }

    fprintf(stderr, "K-mer Enrichment Analysis (v2 with TF-DF)\n");
    fprintf(stderr, "==========================================\n");
    fprintf(stderr, "Genome:     %s\n", genome_file);
    fprintf(stderr, "CDR BED:    %s\n", cdr_bed);
    fprintf(stderr, "Centro BED: %s\n", centromere_bed);
    fprintf(stderr, "Output:     %s\n", output_file);
    fprintf(stderr, "Max k:      %d\n", max_k);
    fprintf(stderr, "Min count:  %d\n", min_count);
    fprintf(stderr, "Threads:    %d\n\n", num_threads);

    // Read genome
    if (read_genome(genome_file) < 0) return 1;

    // Read BED files
    if (read_bed(cdr_bed, &cdr_regions, &num_cdr_regions) < 0) return 1;
    if (read_bed(centromere_bed, &centromere_regions, &num_centromere_regions) < 0) return 1;

    // Compute flanking regions
    compute_flanking_regions();

    // Calculate bitmap sizes
    cdr_bitmap_words = (num_cdr_regions + 63) / 64;
    flank_bitmap_words = (num_flank_regions + 63) / 64;

    fprintf(stderr, "CDR regions: %d, Flanking regions: %d\n", num_cdr_regions, num_flank_regions);
    fprintf(stderr, "Bitmap sizes: CDR=%d words, Flank=%d words\n\n", cdr_bitmap_words, flank_bitmap_words);

    // Initialize k-mer tables
    fprintf(stderr, "Initializing k-mer tables...\n");
    for (int k = 1; k <= max_k; k++) {
        init_kmer_table(k);
    }

    // Count k-mers in CDR regions
    fprintf(stderr, "\nCounting k-mers in CDR regions...\n");
    pthread_t *threads = malloc(sizeof(pthread_t) * num_threads);
    ThreadArg *thread_args = malloc(sizeof(ThreadArg) * num_threads);

    for (int i = 0; i < num_threads; i++) {
        thread_args[i].thread_id = i;
        thread_args[i].num_threads = num_threads;
        thread_args[i].is_cdr = 1;
        thread_args[i].min_k = 1;
        thread_args[i].max_k = max_k;
        pthread_create(&threads[i], NULL, count_kmers_cdr_thread, &thread_args[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Count k-mers in flanking regions
    fprintf(stderr, "Counting k-mers in flanking regions...\n");
    for (int i = 0; i < num_threads; i++) {
        thread_args[i].is_cdr = 0;
        pthread_create(&threads[i], NULL, count_kmers_flank_thread, &thread_args[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Write results
    fprintf(stderr, "\nWriting results...\n");
    write_results(output_file);

    // Cleanup
    free(threads);
    free(thread_args);
    for (int i = 0; i < num_chromosomes; i++) {
        if (chromosomes[i].seq) free(chromosomes[i].seq);
    }
    free(chromosomes);
    free(cdr_regions);
    free(centromere_regions);
    free(flank_regions);

    for (int k = 1; k <= max_k; k++) {
        KmerTable *table = &kmer_tables[k];
        for (uint64_t i = 0; i < table->capacity; i++) {
            if (table->entries[i].cdr_bitmap) free(table->entries[i].cdr_bitmap);
            if (table->entries[i].flank_bitmap) free(table->entries[i].flank_bitmap);
        }
        free(table->entries);
        pthread_mutex_destroy(&table_mutexes[k]);
    }

    fprintf(stderr, "Done!\n");
    return 0;
}
