# ArraySplitter C

High-performance C implementation of ArraySplitter with multithreading support.

## Features

- **Parallel processing** - Each array processed in separate thread
- **Same algorithm** - Identical to Python version (FS-tree, compute_cuts)
- **Compatible output** - Same file formats as Python version
- **Gzip support** - Read .fa.gz files directly

## Building

```bash
cd src/c
make            # Build optimized binary
make debug      # Build with debug symbols
make clean      # Remove build artifacts
```

### Requirements

- GCC with C11 support
- pthreads
- zlib (for gzip support)

On macOS, these are available via Xcode Command Line Tools.
On Linux: `apt install build-essential zlib1g-dev`

## Installation

```bash
make install        # Copy to ../../bin/
make install-user   # Copy to ~/.local/bin/ or ~/bin/
```

## Usage

```bash
arraysplitter_c -i input.fasta -o output_prefix [options]

Required:
  -i, --input FILE     Input FASTA file (supports .gz)
  -o, --output PREFIX  Output file prefix

Options:
  -t, --threads N      Number of threads (default: 4)
  -d, --depth N        FS-tree depth (default: 100)
  -v, --verbose        Verbose output
  -h, --help           Show help
```

### Examples

```bash
# Basic usage with 8 threads
arraysplitter_c -i arrays.fa -o result -t 8

# Gzipped input
arraysplitter_c -i arrays.fa.gz -o result

# Verbose mode for debugging
arraysplitter_c -i arrays.fa -o result -v
```

## Output Files

| File | Description |
|------|-------------|
| `PREFIX.decomposed.fasta` | Monomers in FASTA format (space-separated) |
| `PREFIX.monomers.tsv` | Detailed TSV with type, length, sequence |
| `PREFIX.lengths` | Monomer lengths only |

### FASTA Header Format

```
>original_header cut=CUTSEQ orientation=fwd/rev n_monomers=N range=MIN-MAX avg=AVG
```

### TSV Columns

```
sequence_id  orientation  index  type  length  is_flank  sequence
```

Types: `LEFT_FLANK`, `MONOMER`, `RIGHT_FLANK`

## Performance

Benchmark on Apple M1 (8 cores):

| Dataset | Arrays | Size | Time | CPU Usage |
|---------|--------|------|------|-----------|
| rotation.fa | 20 | 0.5 MB | <1s | 400% |
| zebra_finch | 429 | 60 MB | 16s | 590% |
| satellome | 875 | - | - | - |

## Algorithm

See [ALGORITHM.md](ALGORITHM.md) for detailed algorithm documentation.

**Summary:**
1. Orient sequence to canonical form (A>T, C>G)
2. Build Frequency Suffix Tree to find frequent subsequences
3. Evaluate each candidate as potential cut sequence
4. Select best cut based on regularity score
5. Split array at cut positions

## Differences from Python Version

| Aspect | Python | C |
|--------|--------|---|
| Threading | Single-threaded | Multi-threaded (pthreads) |
| Memory | Higher (Python overhead) | Lower (direct allocation) |
| Speed | Baseline | ~5-10x faster |
| Anchor Graph | Full implementation | Simplified (direct cut) |

The C version uses a simplified decomposition approach without the full Anchor Graph analysis. For most satellite arrays, results are equivalent.

## File Structure

```
src/c/
├── arraysplitter.h   # Header with structures and declarations
├── main.c            # CLI, threading, output writing
├── fasta.c           # FASTA reader (with gzip)
├── sequence.c        # Sequence utilities
├── fstree.c          # Frequency Suffix Tree
├── decompose.c       # Decomposition logic
├── Makefile          # Build system
├── README.md         # This file
└── ALGORITHM.md      # Algorithm documentation
```

## License

Same as main ArraySplitter project.
