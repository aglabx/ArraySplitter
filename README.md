# ArraySplitter: De Novo Decomposition of Satellite DNA Arrays

Decomposes satellite DNA arrays into monomers within telomere-to-telomere (T2T) assemblies. Ideal for analyzing centromeric and pericentromeric regions on monomeric level.

**Status:** Production ready. Successfully handles arrays from kilobase to megabase scale.

**Key Features:**
- De novo monomer identification without prior knowledge
- Autocorrelation-based period detection for robust periodicity analysis
- Automatic orientation to canonical form (A>T, C>G)
- Deterministic output sorted by genomic coordinates
- Multi-threaded processing

**Performance:** CHM13v2.0 assembly (~1300 alpha satellite arrays) processes in ~3.5 minutes (16 threads)

## Installation

```bash
pip install arraysplitter
```

Or build from source:
```bash
cd src/rust/arraysplitter
cargo build --release
```

## Quick Start

```bash
# Basic decomposition
arraysplitter -i arrays.fa -o output_prefix -t 16

# With predefined cut sequences
arraysplitter -i arrays.fa -o output_prefix -c ATG,CGCG -t 16

# Show version
arraysplitter --version
```

## Output Files

All output is **deterministically sorted by chromosome and genomic position** (chr1 → chr22 → chrX → chrY → chrM).

| File | Description |
|------|-------------|
| `.decomposed.fasta` | Monomers with orientation info in headers |
| `.hors.tsv` | HOR-level decomposition with metrics per HOR monomer |
| `.monomers.tsv` | Base-level monomers from recursive HOR decomposition |
| `.lengths` | Fragment lengths for each array |

### HORs TSV Columns (`.hors.tsv`)

Contains the primary decomposition into HOR (Higher Order Repeat) monomers.

| Column | Description |
|--------|-------------|
| `array_id` | Array identifier (chr_start_end_len_period_type) |
| `type` | `pred_array`, `monomer`, `flank`, `array`, `consensus` |
| `idx` | Monomer index within array |
| `length` | Sequence length |
| `source` | Detection method (`anchor`, `split_2x`, etc.) |
| `ed_tmpl` | Edit distance to consensus template |
| `ed_prev` | Edit distance to previous monomer |
| `ed_next` | Edit distance to next monomer |
| `period` | Detected repeat period |
| `autocorr` | Autocorrelation value at period |
| `cut_sequence` | Anchor sequence used for splitting |
| `orientation` | `fwd` or `rev` (reverse complemented) |
| `sequence` | Actual DNA sequence |

### Monomers TSV Columns (`.monomers.tsv`)

Contains base-level monomers after recursive decomposition of HORs. Each HOR is recursively decomposed until no further periodicity is detected (autocorrelation ≤ 0.5) or minimum length (5bp) is reached.

| Column | Description |
|--------|-------------|
| `array_id` | Array identifier |
| `hor_idx` | Index of parent HOR from primary decomposition |
| `sub_idx` | Index within parent HOR (hierarchical for nested decomposition) |
| `level` | Recursion depth (1 = direct child of HOR) |
| `length` | Sequence length |
| `period` | Detected period at this level (0 if base monomer) |
| `autocorr` | Autocorrelation value at detected period |
| `source` | `recursive_anchor`, `recursive_split`, `base`, `recursive_flank` |
| `sequence` | Actual DNA sequence |

## Algorithm

ArraySplitter employs an autocorrelation-based algorithm for detecting repeat periods and decomposing satellite DNA arrays.

### 1. Canonical Orientation

Arrays are oriented to canonical form:
- Primary rule: A > T (more A's than T's)
- Secondary rule: C > G (if A=T)
- Non-canonical arrays are reverse complemented

### 2. Period Detection via Autocorrelation

The algorithm computes sequence autocorrelation to detect periodicity:

```
autocorr(offset) = matches / comparisons
```

Where `matches` counts identical nucleotides at positions `i` and `i + offset`.

**Key innovations:**
- **Random expectation correction:** Subtracts expected random match rate based on nucleotide composition
- **Refined period search:** Uses FFT-like peak detection to find true period vs harmonics
- **Confidence scoring:** Autocorrelation excess over random indicates detection confidence

### 3. Anchor Selection

For the detected period, finds optimal anchor (cut sequence) using:

1. **K-mer enumeration:** Extract all k-mers (k=10 by default) from the sequence
2. **Position analysis:** For each k-mer, record all occurrence positions
3. **Scoring metrics:**
   - **Uniqueness:** Fraction of occurrences exactly `period` apart
   - **Regularity:** How evenly spaced the occurrences are
4. **Combined score:** `uniqueness × regularity`
5. **Deterministic selection:** K-mers sorted lexicographically for reproducible tie-breaking

### 4. Array Decomposition

Using the selected anchor:
1. Split array at all anchor occurrences
2. First fragment → left flank (if < 70% of period)
3. Middle fragments → monomers
4. Last fragment → right flank (if < 70% of period)
5. Apply heuristics for multiplet splitting (doublets, triplets, etc.)

### 5. Output Generation

Results are:
- Sorted by chromosome (natural order: 1, 2, ..., 22, X, Y, M)
- Within chromosome, sorted by start position
- Fully deterministic across runs

## Methods

### `autocorr` (Default)

Uses autocorrelation for period detection. Best for:
- Regular tandem repeats
- Alpha satellite arrays
- HOR (Higher Order Repeat) structures

### `classic`

Uses frequency suffix tree approach. Better for:
- Irregular or degenerate repeats
- Very short arrays
- Arrays with high mutation rates

### `both`

Tries autocorrelation first, falls back to classic if autocorr fails.

## Command Line Options

```
arraysplitter --help

Options:
  -i, --input <FILE>       Input FASTA file
  -o, --output <PREFIX>    Output prefix
  -t, --threads <N>        Number of threads [default: all cores]
  -c, --cuts <SEQ,SEQ>     Predefined cut sequences (comma-separated)
  -d, --depth <N>          Max depth for cut search [default: 100]
  --method <METHOD>        Detection method: autocorr, classic, both [default: autocorr]
  --max-ed-len <N>         Max monomer length for edit distance [default: 10000]
  --stats                  Print detailed statistics
  --top-outliers <N>       Number of outliers to show [default: 10]
  -V, --version            Print version
```

## Citation

If you use ArraySplitter in your research, please cite:
[Publication pending]

## Contact

For questions or support: ad3002@gmail.com
