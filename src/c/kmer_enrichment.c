/*
 * K-mer Enrichment Analysis for Centromeric Regions
 *
 * Compares k-mer frequencies between CDR (deep centromere) and flanking regions.
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
} Region;

// Chromosome sequence
typedef struct {
    char name[64];
    char *seq;
    uint64_t length;
    uint64_t offset;  // offset in file for fai
} ChromSeq;

// K-mer count entry
typedef struct {
    uint64_t kmer;      // encoded k-mer
    uint64_t cdr_count;
    uint64_t flank_count;
    int k;
} KmerCount;

// Global data
ChromSeq *chromosomes = NULL;
int num_chromosomes = 0;

Region *cdr_regions = NULL;
int num_cdr_regions = 0;

Region *centromere_regions = NULL;
int num_centromere_regions = 0;

// K-mer hash tables (one per k value)
// Using simple arrays for k <= 14, hash tables for larger
typedef struct {
    KmerCount *entries;
    uint64_t capacity;
    uint64_t count;
} KmerTable;

KmerTable kmer_tables[MAX_KMER + 1];
pthread_mutex_t table_mutexes[MAX_KMER + 1];

// Thread arguments
typedef struct {
    int thread_id;
    int num_threads;
    int region_type;  // 0 = CDR, 1 = flanking
    int min_k;
    int max_k;
} ThreadArg;

// Parameters
int max_k = 20;
int min_count = 40;
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

// Encode k-mer to uint64
static inline uint64_t encode_kmer(const char *seq, int k) {
    uint64_t kmer = 0;
    for (int i = 0; i < k; i++) {
        int code = base_to_code(seq[i]);
        if (code < 0) return UINT64_MAX;  // invalid
        kmer = (kmer << 2) | code;
    }
    return kmer;
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
        rc = (rc << 2) | (3 - (kmer & 3));  // complement and reverse
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
    // MurmurHash-like mixing
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
        capacity = 1ULL << (2 * k);  // exact size for small k
    } else {
        capacity = 1ULL << 28;  // ~256M entries for large k
    }

    kmer_tables[k].entries = calloc(capacity, sizeof(KmerCount));
    kmer_tables[k].capacity = capacity;
    kmer_tables[k].count = 0;
    pthread_mutex_init(&table_mutexes[k], NULL);
}

// Add or update k-mer count
void add_kmer(int k, uint64_t kmer, int is_cdr) {
    KmerTable *table = &kmer_tables[k];
    uint64_t idx;

    if (k <= 14) {
        // Direct indexing for small k
        idx = kmer;
    } else {
        // Hash table lookup for large k
        idx = hash_kmer(kmer, table->capacity);

        // Linear probing
        while (table->entries[idx].k != 0 && table->entries[idx].kmer != kmer) {
            idx = (idx + 1) % table->capacity;
        }
    }

    pthread_mutex_lock(&table_mutexes[k]);

    if (table->entries[idx].k == 0) {
        table->entries[idx].kmer = kmer;
        table->entries[idx].k = k;
        table->count++;
    }

    if (is_cdr) {
        table->entries[idx].cdr_count++;
    } else {
        table->entries[idx].flank_count++;
    }

    pthread_mutex_unlock(&table_mutexes[k]);
}

// Read FASTA index
int read_fai(const char *fai_path) {
    FILE *fp = fopen(fai_path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open FAI file: %s\n", fai_path);
        return -1;
    }

    chromosomes = malloc(sizeof(ChromSeq) * 1000);
    char line[MAX_LINE];

    while (fgets(line, MAX_LINE, fp)) {
        ChromSeq *chr = &chromosomes[num_chromosomes];
        uint64_t line_bases, line_bytes;

        if (sscanf(line, "%63s\t%llu\t%llu\t%llu\t%llu",
                   chr->name, (unsigned long long*)&chr->length, (unsigned long long*)&chr->offset,
                   (unsigned long long*)&line_bases, (unsigned long long*)&line_bytes) >= 3) {
            chr->seq = NULL;
            num_chromosomes++;
        }
    }

    fclose(fp);
    fprintf(stderr, "Loaded %d chromosomes from FAI\n", num_chromosomes);
    return 0;
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
            // New chromosome
            if (current_chr >= 0) {
                chromosomes[current_chr].length = seq_pos;
            }

            current_chr = num_chromosomes++;
            ChromSeq *chr = &chromosomes[current_chr];

            // Parse chromosome name
            char *name_end = strchr(line + 1, ' ');
            if (!name_end) name_end = strchr(line + 1, '\n');
            if (!name_end) name_end = strchr(line + 1, '\t');

            int name_len = name_end ? (name_end - line - 1) : strlen(line + 1);
            if (name_len > 63) name_len = 63;
            strncpy(chr->name, line + 1, name_len);
            chr->name[name_len] = '\0';

            // Remove trailing whitespace
            while (name_len > 0 && (chr->name[name_len-1] == '\n' || chr->name[name_len-1] == '\r')) {
                chr->name[--name_len] = '\0';
            }

            // Allocate sequence buffer
            seq_capacity = 300000000;  // 300MB initial
            chr->seq = malloc(seq_capacity);
            seq_pos = 0;

            fprintf(stderr, "  Reading %s...\n", chr->name);
        } else if (current_chr >= 0) {
            // Sequence line
            ChromSeq *chr = &chromosomes[current_chr];
            size_t len = strlen(line);

            // Remove newline
            while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
                line[--len] = '\0';
            }

            // Expand buffer if needed
            if (seq_pos + len >= seq_capacity) {
                seq_capacity *= 2;
                chr->seq = realloc(chr->seq, seq_capacity);
            }

            // Convert to uppercase and copy
            for (size_t i = 0; i < len; i++) {
                chr->seq[seq_pos++] = toupper(line[i]);
            }
        }
    }

    // Finalize last chromosome
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
        if (sscanf(line, "%63s\t%llu\t%llu", r->chrom, (unsigned long long*)&r->start, (unsigned long long*)&r->end) >= 3) {
            (*num_regions)++;
        }
    }

    fclose(fp);
    fprintf(stderr, "Loaded %d regions from %s\n", *num_regions, bed_path);
    return 0;
}

// Check if position is in CDR region
int is_in_cdr(const char *chrom, uint64_t pos) {
    for (int i = 0; i < num_cdr_regions; i++) {
        if (strcmp(cdr_regions[i].chrom, chrom) == 0 &&
            pos >= cdr_regions[i].start && pos < cdr_regions[i].end) {
            return 1;
        }
    }
    return 0;
}

// Process sequence for k-mers
void process_sequence(const char *seq, uint64_t len, int is_cdr, int min_k, int max_k_local) {
    if (len < (uint64_t)min_k) return;

    for (int k = min_k; k <= max_k_local && k <= (int)len; k++) {
        // Sliding window
        uint64_t valid_pos = 0;
        uint64_t current_kmer = 0;

        for (uint64_t i = 0; i < len; i++) {
            int code = base_to_code(seq[i]);

            if (code < 0) {
                // Reset on N
                valid_pos = 0;
                current_kmer = 0;
                continue;
            }

            current_kmer = ((current_kmer << 2) | code) & ((1ULL << (2*k)) - 1);
            valid_pos++;

            if (valid_pos >= (uint64_t)k) {
                uint64_t canon = canonical_kmer(current_kmer, k);
                add_kmer(k, canon, is_cdr);
            }
        }
    }
}

// Thread worker function
void* count_kmers_thread(void *arg) {
    ThreadArg *targ = (ThreadArg*)arg;
    int is_cdr = (targ->region_type == 0);

    // Select regions
    Region *regions = is_cdr ? cdr_regions : centromere_regions;
    int num_regions_local = is_cdr ? num_cdr_regions : num_centromere_regions;

    for (int i = targ->thread_id; i < num_regions_local; i += targ->num_threads) {
        Region *r = &regions[i];
        ChromSeq *chr = find_chromosome(r->chrom);

        if (!chr || !chr->seq) {
            fprintf(stderr, "Warning: Chromosome %s not found\n", r->chrom);
            continue;
        }

        uint64_t start = r->start;
        uint64_t end = r->end;
        if (end > chr->length) end = chr->length;
        if (start >= end) continue;

        if (is_cdr) {
            // Process CDR region directly
            process_sequence(chr->seq + start, end - start, 1, targ->min_k, targ->max_k);
        } else {
            // For centromere regions, process flanking parts (exclude CDR)
            // Find overlapping CDR regions and process gaps
            uint64_t current_pos = start;

            for (int j = 0; j < num_cdr_regions; j++) {
                if (strcmp(cdr_regions[j].chrom, r->chrom) != 0) continue;

                uint64_t cdr_start = cdr_regions[j].start;
                uint64_t cdr_end = cdr_regions[j].end;

                // Check if CDR overlaps with this centromere region
                if (cdr_end <= start || cdr_start >= end) continue;

                // Clip CDR to centromere boundaries
                if (cdr_start < start) cdr_start = start;
                if (cdr_end > end) cdr_end = end;

                // Process flanking before CDR
                if (current_pos < cdr_start) {
                    process_sequence(chr->seq + current_pos, cdr_start - current_pos, 0, targ->min_k, targ->max_k);
                }

                current_pos = cdr_end;
            }

            // Process remaining flanking after last CDR
            if (current_pos < end) {
                process_sequence(chr->seq + current_pos, end - current_pos, 0, targ->min_k, targ->max_k);
            }
        }
    }

    return NULL;
}

// Comparison function for sorting by enrichment
int compare_enrichment(const void *a, const void *b) {
    const KmerCount *ka = (const KmerCount*)a;
    const KmerCount *kb = (const KmerCount*)b;

    // Calculate log2 fold change
    double fc_a = (ka->cdr_count + 1.0) / (ka->flank_count + 1.0);
    double fc_b = (kb->cdr_count + 1.0) / (kb->flank_count + 1.0);

    if (fc_b > fc_a) return 1;
    if (fc_b < fc_a) return -1;
    return 0;
}

// Write results
void write_results(const char *output_path) {
    FILE *fp = fopen(output_path, "w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open output file: %s\n", output_path);
        return;
    }

    // Header
    fprintf(fp, "k\trank\tkmer\tcdr_count\tflank_count\ttotal\tlog2_fc\tenrichment\n");

    // Process each k
    for (int k = 1; k <= max_k; k++) {
        KmerTable *table = &kmer_tables[k];

        // Collect significant k-mers
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

        // Sort by enrichment (CDR/flanking ratio)
        qsort(significant, num_significant, sizeof(KmerCount), compare_enrichment);

        // Output top 20
        int top_n = num_significant < 20 ? num_significant : 20;
        char kmer_str[MAX_KMER + 1];

        for (int i = 0; i < top_n; i++) {
            KmerCount *kc = &significant[i];
            decode_kmer(kc->kmer, k, kmer_str);

            uint64_t total = kc->cdr_count + kc->flank_count;
            double fc = (kc->cdr_count + 1.0) / (kc->flank_count + 1.0);
            double log2_fc = log2(fc);
            const char *enrichment = (log2_fc > 0.5) ? "CDR" : (log2_fc < -0.5) ? "Flanking" : "Neutral";

            fprintf(fp, "%d\t%d\t%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.3f\t%s\n",
                    k, i + 1, kmer_str, kc->cdr_count, kc->flank_count, total, log2_fc, enrichment);
        }

        // Also output top 20 flanking-enriched (reverse order)
        for (int i = num_significant - 1; i >= 0 && i >= num_significant - 20; i--) {
            KmerCount *kc = &significant[i];
            double fc = (kc->cdr_count + 1.0) / (kc->flank_count + 1.0);
            double log2_fc = log2(fc);

            if (log2_fc >= -0.5) break;  // Stop if not flanking-enriched

            decode_kmer(kc->kmer, k, kmer_str);
            uint64_t total = kc->cdr_count + kc->flank_count;

            fprintf(fp, "%d\t%d\t%s\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%.3f\tFlanking\n",
                    k, -(num_significant - i), kmer_str, kc->cdr_count, kc->flank_count, total, log2_fc);
        }

        free(significant);
        fprintf(stderr, "k=%d: %d significant k-mers\n", k, num_significant);
    }

    fclose(fp);
    fprintf(stderr, "Results written to %s\n", output_path);
}

void print_usage(const char *prog) {
    fprintf(stderr, "K-mer Enrichment Analysis for Centromeric Regions\n\n");
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
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "  TSV file with columns: k, rank, kmer, cdr_count, flank_count, total, log2_fc, enrichment\n");
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

    // Validate parameters
    if (!genome_file || !cdr_bed || !centromere_bed || !output_file) {
        fprintf(stderr, "Error: Missing required arguments\n\n");
        print_usage(argv[0]);
        return 1;
    }

    if (max_k < 1 || max_k > MAX_KMER) {
        fprintf(stderr, "Error: max-k must be between 1 and %d\n", MAX_KMER);
        return 1;
    }

    fprintf(stderr, "K-mer Enrichment Analysis\n");
    fprintf(stderr, "=========================\n");
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

    // Initialize k-mer tables
    fprintf(stderr, "\nInitializing k-mer tables...\n");
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
        thread_args[i].region_type = 0;  // CDR
        thread_args[i].min_k = 1;
        thread_args[i].max_k = max_k;
        pthread_create(&threads[i], NULL, count_kmers_thread, &thread_args[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Count k-mers in flanking regions
    fprintf(stderr, "Counting k-mers in flanking regions...\n");
    for (int i = 0; i < num_threads; i++) {
        thread_args[i].region_type = 1;  // Flanking
        pthread_create(&threads[i], NULL, count_kmers_thread, &thread_args[i]);
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
    for (int k = 1; k <= max_k; k++) {
        free(kmer_tables[k].entries);
        pthread_mutex_destroy(&table_mutexes[k]);
    }

    fprintf(stderr, "Done!\n");
    return 0;
}
