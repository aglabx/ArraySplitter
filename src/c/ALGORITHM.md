# ArraySplitter C - Algorithm Documentation

## Overview

ArraySplitter decomposes satellite DNA arrays into monomers using a frequency-based suffix tree approach. The algorithm finds optimal "cut sequences" that appear regularly throughout the array, then splits the array at these positions.

## Pipeline

```
Input FASTA
    │
    ▼
┌─────────────────────┐
│ 1. Read & Orient    │  Canonical orientation (A>T, C>G)
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│ 2. Build FS-Tree    │  Find frequent subsequences (hints)
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│ 3. Evaluate Cuts    │  Score each hint as potential cut
└─────────────────────┘
    │
    ▼
┌─────────────────────┐
│ 4. Decompose        │  Split array at cut positions
└─────────────────────┘
    │
    ▼
Output Files (.decomposed.fasta, .monomers.tsv, .lengths)
```

## Algorithm Details

### 1. Canonical Orientation

Before processing, each array is oriented to a canonical form:
- Count nucleotides A, T, C, G
- If T > A: reverse complement the sequence
- If T == A and G > C: reverse complement

This ensures consistent results regardless of input strand.

```c
bool seq_is_canonical(const char *seq, size_t len) {
    // Count A, T, C, G
    // Return true if A > T, or (A == T and C >= G)
}
```

### 2. Frequency Suffix Tree (FS-Tree)

The FS-tree identifies frequent subsequences that are candidates for cut sequences.

**Concept:**
- Start from each nucleotide position
- Extend one nucleotide at a time
- Track how many positions share the same prefix
- When count drops below cutoff, stop extending

**Data Structure:**
```
HeapNode {
    level: current depth (subsequence length)
    names: starting positions of matching sequences
    positions: current positions in sequence
    count: number of matching positions
}
```

**Algorithm (BFS with min-heap by level):**
```
1. For each starting nucleotide (A, C, G, T):
   - Find all positions where it occurs
   - Push to heap with level=1

2. While heap not empty:
   - Pop node with smallest level
   - If level changed, emit best hint from previous level
   - Group positions by next nucleotide (A, C, G, T)
   - For each group with count > cutoff:
     - Push new node with level+1
```

**Cutoff Calculation:**
```c
size_t cutoff = len > 1000000 ? 1000 :
                len > 100000  ? 250  :
                len > 10000   ? 10   : 3;
```

**Self-Repeating Pattern Detection:**
If a hint like "ATAT" is found, check if it's composed of smaller repeats ("AT").
Use the minimal unit as the actual hint.

### 3. Cut Evaluation

Each hint is evaluated as a potential cut sequence by:

1. **Split the array** at all occurrences of the hint
2. **Calculate periods** (length of each segment + cut length)
3. **Find mode period** (most common segment length)
4. **Score calculation:**

```c
base_score = mode_count / total_segments;  // 0.0 to 1.0

// Fragmentation penalty (short fragments indicate overcutting)
short_threshold = mode_period * 0.5;
fragmentation = short_fragments / total_segments;

adjusted_score = base_score * (1.0 - fragmentation * 0.5);
```

**Selection:**
- Candidates within 0.05 of best base_score are considered
- Sort by adjusted_score (descending), then by period (ascending)
- Select the best candidate

### 4. Decomposition

Once the best cut is selected:

1. Find all non-overlapping positions of cut sequence
2. Split array at these positions
3. Prepend cut sequence to each segment (except left flank)

**Handling Overlapping Matches:**
```c
// Filter overlapping positions
for each position:
    if position >= last_end:
        keep position
        last_end = position + cut_length
```

**Output Structure:**
```
[LEFT_FLANK] [CUT+SEG1] [CUT+SEG2] ... [CUT+SEGn]
```

- First segment (before first cut) = LEFT_FLANK
- Last segment (if short) = RIGHT_FLANK
- All others = MONOMER

## Multithreading

The C implementation uses pthreads for parallel processing:

```
┌─────────────┐
│ Main Thread │
└─────────────┘
       │
       ▼
┌─────────────────────────────────┐
│ Create Work Queue (all arrays)  │
└─────────────────────────────────┘
       │
       ▼
┌─────────────────────────────────┐
│ Spawn N Worker Threads          │
└─────────────────────────────────┘
       │
       ├──► Thread 1 ──► Process array[i]
       ├──► Thread 2 ──► Process array[j]
       ├──► Thread 3 ──► Process array[k]
       └──► ...
       │
       ▼
┌─────────────────────────────────┐
│ Join Threads, Write Results     │
└─────────────────────────────────┘
```

**Thread-Safety:**
- Each array is processed independently
- Work queue protected by mutex
- Results stored in pre-allocated array (no conflicts)

## Complexity Analysis

For an array of length N with M hints:

| Step | Time Complexity | Space Complexity |
|------|-----------------|------------------|
| Orientation | O(N) | O(N) for revcomp |
| FS-Tree | O(N × D) | O(N) positions |
| Cut Evaluation | O(M × N) | O(N) per candidate |
| Decomposition | O(N) | O(N) for result |

Where D = depth parameter (default 100)

**Total:** O(N × D + M × N) ≈ O(M × N) for typical parameters

## Output Files

### .decomposed.fasta
```
>header cut=SEQUENCE orientation=fwd/rev n_monomers=N range=MIN-MAX avg=AVG
MONO1 MONO2 MONO3 ... (space-separated)
```

### .monomers.tsv
```
sequence_id  orientation  index  type        length  is_flank  sequence
header       fwd          0      LEFT_FLANK  123     TRUE      ACGT...
header       fwd          1      MONOMER     200     FALSE     ACGT...
```

### .lengths
```
>header cut=SEQUENCE orientation=fwd/rev n_monomers=N range=MIN-MAX avg=AVG
123 200 200 198 45 (space-separated lengths)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-t, --threads` | 4 | Number of parallel threads |
| `-d, --depth` | 100 | Maximum depth for FS-tree exploration |
| `-v, --verbose` | off | Print detailed progress |

## References

The algorithm is based on the frequency suffix tree approach for finding conserved regions in tandem repeats, adapted for satellite DNA decomposition.

Key insight: In satellite arrays, the monomer boundary sequence appears at regular intervals. By finding the most frequent subsequence with consistent spacing, we identify the optimal cut point.
