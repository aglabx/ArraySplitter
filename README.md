# ArraySplitter: De Novo Decomposition of Satellite DNA Arrays

Decomposes satellite DNA arrays into monomers within telomere-to-telomere (T2T) assemblies. Ideal for analyzing centromeric and pericentromeric regions on monomeric level.

**Status:** Production ready. Successfully handles arrays from kilobase to megabase scale.

**Key Features:**
- De novo monomer identification without prior knowledge
- Automatic orientation to canonical form (A>T, C>G)
- Integrated monomer rotation for standardized comparison
- Classification of arrays into repeat families
- Structural importance scoring to identify centromeric and other functional regions

**Performance:** CHM13v2.0 assembly (13K arrays > 1Kb) processes in ~3 hours (single-threaded)

## Installation

**Prerequisites**

* Python 3.6 or later

**Installation with pip:**

```bash
pip install arraysplitter
```

## Tool Overview

ArraySplitter provides a unified command-line interface with multiple subcommands:

### Main Command: `arraysplitter`

```bash
arraysplitter <command> [options]
```

Available commands:
- `split` - Decompose arrays into monomers
- `classify` - Group arrays into families based on patterns
- `rotate` - Normalize monomer orientation
- `extract` - Extract and count unique monomers

### 1. `arraysplitter split` - Array Decomposition

Performs de novo decomposition of satellite DNA arrays into individual monomers.

**Usage:**
```bash
# Basic decomposition
arraysplitter split -i arrays.fa -o output_prefix

# With predefined cut sequences
arraysplitter split -i arrays.fa -o output_prefix -c ATG,CGCG
```

**Output files:**
- `.decomposed.fasta` - Monomers with orientation and rotation applied
- `.monomers.tsv` - Detailed table with monomer information
- `.lengths` - Fragment lengths for each array

**Features:**
- Automatic canonical orientation (A>T, C>G)
- Integrated rotation (monomers start with cut sequence)
- Flank detection and labeling
- Perfect reconstruction guarantee

### 2. `arraysplitter classify` - Family Classification

Groups arrays into families based on cut sequences and decomposition patterns. Identifies structurally important regions.

**Usage:**
```bash
# Classify using lengths file from split step
arraysplitter classify -i output_prefix.lengths -o classification

# With custom similarity threshold
arraysplitter classify -i output_prefix.lengths -o classification -s 0.9
```

**Output files:**
- `.families.tsv` - Array assignments with stability metrics
- `.family_stats.tsv` - Basic statistics per family
- `.family_summary.tsv` - Detailed family analysis with structural importance scores
- `.features.json` - Complete feature data

**Key metrics:**
- **Step stability**: Consistency of monomer length differences
- **Local variability**: Sliding window analysis (10kb windows)
- **Structural importance score**: Combined metric (0-1) identifying potential functional regions

### 3. `arraysplitter rotate` - Monomer Normalization

Rotates monomers to start with the same sequence pattern.

**Usage:**
```bash
# Auto-detect best rotation
arraysplitter rotate -i monomers.fa -o rotated.fa

# With specific starting sequence
arraysplitter rotate -i monomers.fa -o rotated.fa -s ATTCC
```

### 4. `arraysplitter extract` - Monomer Analysis

Extracts unique monomers and calculates frequencies.

**Usage:**
```bash
arraysplitter extract -i monomers.fa -o stats_prefix
```

**Output format:**
```
<length> <frequency> <sequence>
```

## Complete Workflow Example

```bash
# Step 1: Decompose arrays into monomers
arraysplitter split -i centromere_arrays.fa -o centromere

# Step 2: Classify arrays into families and identify structural regions
arraysplitter classify -i centromere.lengths -o centromere_families

# Step 3: Extract and count unique monomers (optional)
arraysplitter extract -i centromere.decomposed.fasta -o centromere_monomers

# Step 4: Additional rotation if needed (optional)
arraysplitter rotate -i centromere.decomposed.fasta -o centromere.rotated.fa -s ATTCC
```

### Output Analysis

After classification, examine the `.family_summary.tsv` file to identify:
- Families with high structural importance scores (>0.8) - potential centromeric regions
- Families with low step variability - structurally conserved regions
- Most stable arrays within each family - best representatives for further study

## Algorithm Description

ArraySplitter employs a novel de novo algorithm for decomposing satellite DNA arrays into constituent monomers without prior knowledge of the monomer sequences. The algorithm is specifically designed to handle the challenges of centromeric and pericentromeric regions in telomere-to-telomere assemblies.

### Overview

The algorithm uses a frequency suffix tree (fs_tree) approach to identify optimal cut sequences that split tandem repeat arrays into individual monomer units. It handles variable-length monomers (polymorphic repeats) common in biological sequences and automatically orients sequences to canonical form for comparability.

### Detailed Algorithm Steps

#### 1. Canonical Orientation
Arrays are first oriented to canonical form using the rules:
- Primary: A > T (arrays with more A's than T's are kept as-is)
- Secondary: C > G (if A=T, arrays with more C's than G's are kept as-is)
- Arrays not in canonical form are reverse complemented

#### 2. All Nucleotide Analysis
The algorithm analyzes all nucleotides (A, C, T, G) to avoid bias from single nucleotide selection. Each nucleotide serves as an anchor point for building frequency suffix trees.

#### 3. Frequency Suffix Tree Construction
The core data structure is a frequency suffix tree that efficiently identifies repetitive patterns:

- **Starting positions**: All positions containing each nucleotide become root nodes
- **Iterative extension**: The tree grows by extending sequences one nucleotide at a time (A/C/G/T)
- **Frequency filtering**: Only branches exceeding a dynamic cutoff threshold are retained:
  - Arrays > 1MB: cutoff = 1000
  - Arrays > 100KB: cutoff = 250
  - Arrays > 10KB: cutoff = 10
  - Arrays ≤ 10KB: cutoff = 3
- **Heap-based optimization**: Uses a priority queue to efficiently process high-frequency patterns first

#### 4. Candidate Cut Sequence Generation
From the frequency suffix tree, the algorithm extracts potential cut sequences:
- For each sequence length (up to a configurable depth, default 100)
- Identifies the sequence with maximum coverage (highest frequency)
- Generates a ranked list of candidate cut sequences

#### 5. Optimal Cut Selection
The algorithm evaluates each candidate cut sequence by:
- Splitting the array at all occurrences of the cut sequence
- Calculating period distribution between cuts
- Handling perfect/near-perfect repeats (≥80% empty parts between cuts)
- Computing scores with tie-breaking rules:
  1. Prefer cuts producing fewer segments
  2. Prefer smaller fundamental periods (using GCD)
- The cut sequence becomes the START of each monomer (biologically relevant)

#### 6. Monomer Construction
- First fragment (before any cut) is treated as left flank
- Each monomer is built as: cut_sequence + following_part
- Last fragment <70% of average length is treated as right flank
- All monomers are automatically rotated to start with the cut sequence
- Ensures all output monomers are directly comparable

### Key Features

1. **De novo approach**: No prior knowledge of monomer sequences required
2. **Variable-length monomers**: Handles polymorphic repeats naturally
3. **Canonical orientation**: All arrays oriented consistently (A>T, C>G)
4. **Integrated rotation**: Monomers automatically aligned to cut sequence
5. **Family classification**: Groups arrays by cut patterns and variability
6. **Structural scoring**: Identifies centromeric and functional regions

### Performance Characteristics

- **Optimized for**: 100Kb scale arrays
- **Scalability**: Successfully handles megabase-scale arrays (largest ~5 minutes)
- **Benchmarking**: CHM13v20 assembly (13K arrays > 1Kb) processes in ~3 hours
- **Current limitation**: Single-threaded (parallelization planned)

### Classification Algorithm

ArraySplitter includes a classification system that groups arrays into families based on:

1. **Cut sequence identity**: Arrays with different cuts belong to different families
2. **Pattern similarity**: Within same cut, arrays are clustered by:
   - Mean monomer length
   - Length variability
   - Step stability (consistency between consecutive monomers)
3. **Structural importance scoring**:
   - High step stability (consistent monomer spacing)
   - Low local variability (stable regions in 10kb windows)
   - Consistency across arrays in family
   - Score 0-1 (higher = more structurally important)

### Applications

The algorithm is particularly well-suited for:
- Analyzing centromeric and pericentromeric regions in T2T assemblies
- Identifying structurally important genomic regions
- Studying satellite DNA evolution and variation
- Discovering novel tandem repeat families
- Quantifying monomer composition and variability
- Finding stable regions for functional studies

## Contact

For questions or support: ad3002@gmail.com
