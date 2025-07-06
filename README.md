# ArraySplitter: De Novo Decomposition of Satellite DNA Arrays

Decomposes satellite DNA arrays into monomers within telomere-to-telomere (T2T) assemblies. Ideal for analyzing centromeric and pericentromeric regions on monomeric level.

**Status:** In development. Optimized for 100Kb scale arrays; longer arrays will work but may take longer to process. Signigicanlty longer time.

**Update:** From 1.1.6, ArraySplitter now successfully decomposes arrays on megabase scale. Largest arrays takes around 5 minutes to process. Fortunatelly, there are only 41 arrays large 1 Mb in CHM13v20 assembly. And I'm going to add parallel processing to speed up singificantly the process. Currently, it is single-threaded.

**Update:** Monomers are required some polising of borders, I am working on it.

**Update:** To test ArraySplitter, I used CHM13v20 assembly, it takes around 3 hours, to decompose all arrays longer than 1 Kb (13K arrays).

## Installation

**Prerequisites**

* Python 3.6 or later

**Installation with pip:**

```bash
pip install arraysplitter
```

## Tool Overview

ArraySplitter provides three complementary tools for satellite DNA analysis:

### 1. `arraysplitter` - Main Decomposition Tool
The core tool that performs de novo decomposition of satellite DNA arrays into individual monomers.

**Usage:**
```bash
arraysplitter -i chr1.arrays.fa -o chr1.arrays
```

**Purpose:**
- Takes tandem repeat arrays as input
- Identifies optimal cut sequences using frequency suffix tree algorithm
- Splits arrays into constituent monomers
- Outputs monomers separated by spaces within each array

**Example:**
```
Input:  >array1
        CAGCAGCAGCAGCAG
Output: >array1
        CAG CAG CAG CAG CAG
```

### 2. `arraysplitter_rotate` - Monomer Normalization Tool
Rotates monomers to start with the same sequence pattern for standardized comparison.

**Usage:**
```bash
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa
# Or with specific starting sequence:
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa -s ATTCC
```

**Purpose:**
- Normalizes monomer orientation across different arrays
- Essential for comparing monomers from different decompositions
- Finds optimal rotation to match specified pattern
- If no pattern given, automatically selects most common start

**Example:**
```
Input:  >monomer1
        CAGCAGCAG
        >monomer2
        GCAGCAGCA
        >monomer3
        AGCAGCAGC
Output: >monomer1
        CAGCAGCAG
        >monomer2
        CAGCAGCAG
        >monomer3
        CAGCAGCAG
```

### 3. `arraysplitter_extract` - Monomer Analysis Tool
Extracts unique monomers and calculates their frequencies across all arrays.

**Usage:**
```bash
arraysplitter_extract -i arrays.norm.fa -o arrays.stats
```

**Purpose:**
- Identifies all unique monomer sequences
- Counts frequency of each monomer
- Groups by monomer length
- Outputs statistics sorted by frequency

**Output format:**
```
<length> <frequency> <sequence>
```

**Example:**
```
171     1250    ATTCCATTCCATTCC...
171     890     ATTCCATTCTATTCC...
171     245     ATTCCGTTCCATTCC...
169     122     ATTCCATTCCATT...
```

## Complete Workflow Example

```bash
# Step 1: Decompose arrays into monomers
arraysplitter -i centromere_arrays.fa -o centromere_monomers

# Step 2: Normalize monomer orientation
arraysplitter_rotate -i centromere_monomers.fa -o centromere_monomers.normalized.fa

# Step 3: Extract and count unique monomers
arraysplitter_extract -i centromere_monomers.normalized.fa -o centromere_monomers.stats
```

This three-step process provides a complete analysis pipeline from raw tandem arrays to monomer frequency statistics.

## Rotating monomers to start with the same sequence

We found that different arrays of the same repeat family can be decomposed sligtly differently. To make them comparable, ArraySplitter can rotate monomers to start with the same sequence. 

```bash
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa
```

And you can give the sequence to start with:

```bash
arraysplitter_rotate -i arrays.fa -o arrays.norm.fa -s TTTC
```

**Explanation**

* **`-i arrays.fa`:**  FASTA file of monomers.
* **`-o arrays.norm.fa`:** Output FASTA file with rotated monomers.

## Extracting and counting monomers

And finally, you can extract and count monomers from the arrays:

```bash
arraysplitter_extract -i arrays.norm.fa -o arrays.norm
```

It will create a file with monomer length, monomer frequency, and monomer sequence (ordered by frequency). For example, for the arrays.norm.fa file above, the output will be like this:

```bash
514     10      ATCCCATTCC
514     10      GATTGGAGTG
514     6       TCCTTT
514     5       TGCTG
514     10      ATTGAATGGA
514     10      ATGCAATGGA
514     5       TCCTA
```

## Algorithm Description

ArraySplitter employs a novel de novo algorithm for decomposing satellite DNA arrays into constituent monomers without prior knowledge of the monomer sequences. The algorithm is specifically designed to handle the challenges of centromeric and pericentromeric regions in telomere-to-telomere assemblies.

### Overview

The algorithm uses a frequency suffix tree (fs_tree) approach to identify optimal cut sequences that split tandem repeat arrays into individual monomer units. It operates in multiple stages: building a frequency-based suffix tree, identifying candidate cut sequences, selecting the optimal cut, and iteratively refining the decomposition.

### Detailed Algorithm Steps

#### 1. Most Frequent Nucleotide Selection
The algorithm begins by identifying the most frequent nucleotide (A, C, T, or G) in the input array. This nucleotide serves as an anchor point for building the frequency suffix tree, reducing the search space while maintaining effectiveness.

#### 2. Frequency Suffix Tree Construction
The core data structure is a frequency suffix tree that efficiently identifies repetitive patterns:

- **Starting positions**: All positions containing the most frequent nucleotide become root nodes
- **Iterative extension**: The tree grows by extending sequences one nucleotide at a time (A/C/G/T)
- **Frequency filtering**: Only branches exceeding a dynamic cutoff threshold are retained:
  - Arrays > 1MB: cutoff = 1000
  - Arrays > 100KB: cutoff = 250
  - Arrays > 10KB: cutoff = 10
  - Arrays ≤ 10KB: cutoff = 3
- **Heap-based optimization**: Uses a priority queue to efficiently process high-frequency patterns first

#### 3. Candidate Cut Sequence Generation
From the frequency suffix tree, the algorithm extracts potential cut sequences:
- For each sequence length (up to a configurable depth, default 100)
- Identifies the sequence with maximum coverage (highest frequency)
- Generates a ranked list of candidate cut sequences

#### 4. Optimal Cut Selection
The algorithm evaluates each candidate cut sequence by:
- Splitting the array at all occurrences of the cut sequence
- Calculating the period length distribution (period = segment length + cut length)
- Computing a uniformity score: fraction of segments matching the most common period
- Selecting the cut sequence that produces the most uniform period distribution

#### 5. Two-Stage Decomposition Process

**Stage 1 - Initial Decomposition**:
- Splits the array using the selected cut sequence
- Handles edge cases for first and last segments
- For segments matching the expected period: adds directly as monomers
- For longer segments: attempts even splitting based on expected period length

**Stage 2 - Iterative Refinement**:
- Identifies irregular segments (> 1.3× expected period length)
- Uses the most common monomer as a reference template
- For each irregular segment:
  - Tests cut positions at beginning (positions 0-4) and end (last 5 positions)
  - Evaluates cuts using edit distance to the reference monomer
  - Accepts cuts where edit distance < 50% of monomer length
- Iterates until no further improvements are found

### Key Features

1. **De novo approach**: No prior knowledge of monomer sequences required
2. **Adaptive thresholds**: Cutoff values scale with array size for optimal performance
3. **Iterative refinement**: Multiple passes improve boundary precision
4. **Edit distance validation**: Ensures biological relevance of identified monomers

### Performance Characteristics

- **Optimized for**: 100Kb scale arrays
- **Scalability**: Successfully handles megabase-scale arrays (largest ~5 minutes)
- **Benchmarking**: CHM13v20 assembly (13K arrays > 1Kb) processes in ~3 hours
- **Current limitation**: Single-threaded (parallelization planned)

### Applications

The algorithm is particularly well-suited for:
- Analyzing centromeric and pericentromeric regions in T2T assemblies
- Studying satellite DNA evolution and variation
- Identifying novel tandem repeat families
- Quantifying monomer composition in complex arrays

## Contact

For questions or support: ad3002@gmail.com
