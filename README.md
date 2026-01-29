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

## Understanding Key Metrics

### `orientation` — Strand Direction
- **`fwd`** = sequence kept as-is (canonical form: more A's than T's)
- **`rev`** = sequence was reverse-complemented to canonical form

Why it matters: Satellite arrays on opposite strands appear as different sequences. Canonical orientation ensures consistent comparison across the genome.

### `ed_tmpl` — Edit Distance to Template (Consensus)
Measures how different each monomer is from the **consensus sequence** of all monomers in the array.

- **Low ed_tmpl (0-5)** = monomer is very similar to consensus → typical/canonical monomer
- **High ed_tmpl (>20)** = monomer is divergent → possibly a variant, mutation hotspot, or misalignment

Use case: Filter out divergent monomers, identify conserved vs variable positions.

### `ed_prev` / `ed_next` — Edit Distance to Neighbors
Measures how different each monomer is from its **adjacent monomers** (previous and next).

- **Low ed_prev/ed_next** = adjacent monomers are similar → homogeneous region
- **High ed_prev/ed_next** = sudden change → possible HOR boundary, structural variant, or sequence transition

Use case: Detect HOR structure boundaries, identify recombination breakpoints.

### `autocorr` — Autocorrelation Score
Measures periodicity strength (0.0 to 1.0):

- **High autocorr (>0.7)** = strong, regular periodicity → well-conserved tandem repeat
- **Medium autocorr (0.5-0.7)** = detectable but irregular periodicity
- **Low autocorr (<0.5)** = weak or no periodicity → degenerate or non-repetitive

### `cv` — Coefficient of Variation
Standard deviation of monomer lengths divided by mean length:

- **Low cv (<0.1)** = uniform monomer sizes → canonical repeat structure
- **High cv (>0.2)** = variable sizes → insertions/deletions, heterogeneous array

### `ed_per_bp` — Normalized Edit Distance
Edit distance divided by sequence length. Allows comparison across different monomer sizes:

- **0.00-0.02** = highly conserved
- **0.02-0.05** = moderately divergent
- **>0.05** = significantly divergent

## Output Files

All output is **deterministically sorted by chromosome and genomic position** (chr1 → chr22 → chrX → chrY → chrM).

| File | Description |
|------|-------------|
| `.decomposed.fasta` | Monomers with orientation info in headers |
| `.hors.tsv` | HOR-level decomposition (16 columns) |
| `.monomers.tsv` | Base-level monomers from recursive decomposition (17 columns) |
| `.summary.tsv` | One-row-per-array summary with HOR and monomer statistics (23 columns) |
| `.lengths` | Fragment lengths for each array |

### Summary TSV Columns (`.summary.tsv`)

One row per array combining HOR-level and monomer-level statistics. Useful for overview analysis.

| Column | Description |
|--------|-------------|
| `array_id` | Array identifier (chr_start_end_len_period_type) |
| `array_length` | Total array length in bp |
| `orientation` | `fwd` or `rev` (reverse complemented to canonical) |
| `method` | Detection method used (`autocorr`, `classic`) |
| **HOR-level stats** | |
| `hor_period` | Detected HOR period in bp |
| `hor_autocorr` | Autocorrelation at HOR period |
| `hor_n_monomers` | Number of HOR-level monomers |
| `hor_mean_ed_tmpl` | Mean edit distance to HOR consensus |
| `hor_mean_ed_prev` | Mean edit distance between adjacent HORs |
| `hor_cv` | Coefficient of variation for HOR lengths |
| `hor_consensus` | Consensus sequence at HOR level |
| `hor_iupac` | IUPAC ambiguity codes (bases ≥20% frequency) |
| `hor_quality` | Per-position support (digit 0-9, 9=90-100%) |
| **Monomer-level stats** | |
| `mono_period` | Median base monomer period |
| `mono_autocorr` | Mean autocorrelation at monomer level |
| `mono_n_monomers` | Total number of base monomers |
| `mono_mean_ed_tmpl` | Mean edit distance to monomer consensus |
| `mono_mean_ed_prev` | Mean edit distance between adjacent monomers |
| `mono_cv` | Mean coefficient of variation |
| `mono_consensus` | Consensus sequence at monomer level |
| `mono_iupac` | IUPAC ambiguity codes |
| `mono_quality` | Per-position support |
| `cut_sequence` | Anchor k-mer used for splitting |

### HORs TSV Columns (`.hors.tsv`)

Contains the primary decomposition into HOR (Higher Order Repeat) monomers. Multiple rows per array.

**Row types** (in order):
1. `pred_array` - Array-level prediction/header row
2. `flank` - Terminal fragments <70% of period
3. `monomer` - Full HOR monomers (sorted by idx)
4. `array` - Summary statistics row
5. `consensus` - Consensus sequence row

| Column | Description |
|--------|-------------|
| `array_id` | Array identifier (chr_start_end_len_period_type) |
| `type` | `pred_array`, `monomer`, `flank`, `array`, `consensus` |
| `idx` | Monomer index within array (0-based) |
| `length` | Sequence length in bp |
| `source` | Detection method: `anchor`, `split_2x`, `split_3x`, `left_flank`, `right_flank` |
| `ed_tmpl` | Edit distance to consensus template |
| `ed_prev` | Edit distance to previous monomer |
| `ed_next` | Edit distance to next monomer |
| `period` | Detected repeat period in bp |
| `autocorr` | Autocorrelation value at detected period |
| `n_expected` | Expected count of monomers (array_len / period) |
| `ed_per_bp` | Normalized edit distance (ed / length) |
| `cv` | Coefficient of variation for lengths |
| `cut_sequence` | Anchor sequence used for splitting |
| `orientation` | `fwd` or `rev` (reverse complemented) |
| `sequence` | Actual DNA sequence (or `-` for pred_array/array rows) |

### Monomers TSV Columns (`.monomers.tsv`)

Contains base-level monomers after recursive HOR decomposition. **Unified format** matching `.hors.tsv` plus `parent_idx`.

Each HOR is recursively decomposed until:
- No further periodicity detected (autocorrelation ≤ 0.5)
- Minimum length (5bp) reached

**Row types** (in order):
1. `pred_array` - Array-level summary row
2. `base_monomer` - Base-level monomers from recursive decomposition
3. `monomer` - Non-decomposable monomers (e.g., telomeres)

| Column | Description |
|--------|-------------|
| `array_id` | Array identifier |
| `type` | `pred_array`, `base_monomer`, `monomer` |
| `idx` | Global index within array (0-based) |
| `length` | Sequence length in bp |
| `source` | `recursive_anchor`, `recursive_split`, `base`, `recursive_flank` |
| `ed_tmpl` | Edit distance to submonomer consensus |
| `ed_prev` | Edit distance to previous base monomer |
| `ed_next` | Edit distance to next base monomer |
| `period` | Detected period at this level (0 if base) |
| `autocorr` | Autocorrelation value |
| `n_expected` | Always 1 for individual monomers |
| `ed_per_bp` | Normalized edit distance |
| `cv` | Coefficient of variation within parent group |
| `cut_sequence` | Inherited anchor sequence |
| `orientation` | Inherited from array (`fwd`/`rev`) |
| `parent_idx` | Index of parent HOR from `.hors.tsv` |
| `sequence` | Actual DNA sequence |

### Example: α-satellite HOR Decomposition

For a typical α-satellite HOR (512bp → 3×171bp monomers):

**`.hors.tsv`** - 10 HOR monomers (~512bp each):
```
array_id                type        idx  length  period  ...
chr1_centromere         pred_array  10   5120    512     ...
chr1_centromere         monomer     0    512     512     ...
chr1_centromere         monomer     1    512     512     ...
...
chr1_centromere         array       10   5120    512     ...
chr1_centromere         consensus   10   512     512     ... [consensus seq]
```

**`.monomers.tsv`** - 30 base monomers (~171bp each):
```
array_id                type          idx  length  parent_idx  ...
chr1_centromere         pred_array    30   5120    -           ...
chr1_centromere         base_monomer  0    171     0           ...
chr1_centromere         base_monomer  1    171     0           ...
chr1_centromere         base_monomer  2    170     0           ...
chr1_centromere         base_monomer  3    171     1           ...
...
```

**`.summary.tsv`** - Single row with both levels:
```
array_id         length  hor_period  hor_n_monomers  mono_period  mono_n_monomers  ...
chr1_centromere  5120    512         10              171          30               ...
```

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
