# Rotation Algorithm Improvements

## Current Problems

1. **Orientation Logic**: Current code uses simple nucleotide counting which doesn't guarantee canonical form
2. **K-mer Selection**: Looks for k-mers appearing exactly once per monomer, which may not find the most conserved regions
3. **No Integration with Decomposition**: Doesn't use the cut sequences from decomposition, which are proven conserved regions

## Proposed Improvements

### 1. Canonical Orientation (A>T, C>G)

```python
def get_canonical_orientation(sequence):
    """Determine if sequence needs reverse complement for canonical form."""
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    
    # Primary rule: A > T
    if a_count != t_count:
        return a_count > t_count
    
    # Secondary rule: C > G (when A == T)
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    return c_count > g_count
```

### 2. Use Cut Sequences as Anchors

Since ArraySplitter already identifies conserved cut sequences, we should use them:

```python
# In the decomposition output, save cut sequences
>array1 cut=ATGATG n_monomers=10 range=165-175 avg=171.2
ATGATG... ATGATG... ATGATG...

# When rotating, use the cut sequence as anchor
arraysplitter_rotate -i decomposed.fasta -o rotated.fasta --use-cuts
```

### 3. Better Conserved Region Detection

Instead of F1 score based on single occurrence, find truly conserved regions:
- Present in >80% of arrays
- Appears ~1 time per monomer on average
- Longer k-mers (8bp instead of 5bp) for better specificity

## Integration Plan

### Option 1: Modify Existing rotate.py

Add parameters:
- `--canonical`: Force canonical orientation
- `--use-cuts`: Extract and use cut sequences from headers
- `--min-conservation`: Minimum fraction of arrays containing k-mer

### Option 2: New Rotation Mode

Create a post-processing mode that:
1. Reads decomposed FASTA with cut info in headers
2. Orients all sequences canonically
3. Uses cut sequences for rotation
4. Outputs aligned sequences

## Usage Examples

### Current:
```bash
arraysplitter_rotate -i arrays.fa -o rotated.fa -s ATGAT
```

### Proposed:
```bash
# Use cut sequences from decomposition
arraysplitter_rotate -i decomposed.fa -o rotated.fa --use-cuts --canonical

# Or specify conservation threshold
arraysplitter_rotate -i arrays.fa -o rotated.fa --min-conservation 0.8
```

## Benefits

1. **Consistency**: All sequences in same orientation
2. **Biological Relevance**: Uses proven conserved regions (cuts)
3. **Better Alignment**: More reliable rotation anchors
4. **Integration**: Works seamlessly with decomposition output

## Implementation Steps

1. Add canonical orientation function
2. Add header parsing to extract cut sequences
3. Modify k-mer selection to prioritize conservation
4. Add command-line options
5. Update documentation

This would make the rotation step more reliable and better integrated with the overall ArraySplitter workflow.