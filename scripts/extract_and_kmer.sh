#!/bin/bash
# Extract one array and compute k-mer spectrum step by step
# Usage: ./extract_and_kmer.sh

INPUT="/mnt/data/claude/phase1_alpha_extraction/all_alpha_satellite.fasta"
ARRAY_ID="chr7_60370542_60382239_11698_172_MONO"
WORKDIR="/mnt/data/claude/arraysplitter_rust/kmer_experiment"

mkdir -p "$WORKDIR"

# Step 1: Extract the array
echo "=== Step 1: Extracting $ARRAY_ID ==="
awk -v id="$ARRAY_ID" '
    /^>/ { printing = (substr($1, 2) == id) }
    printing { print }
' "$INPUT" > "$WORKDIR/$ARRAY_ID.fasta"

echo "Extracted: $(grep -c "^>" "$WORKDIR/$ARRAY_ID.fasta") sequence(s)"
echo "Length: $(grep -v "^>" "$WORKDIR/$ARRAY_ID.fasta" | tr -d '\n' | wc -c) bp"

# Step 2: Compute k-mer frequencies for each k
echo ""
echo "=== Step 2: K-mer frequency spectrum ==="
python3 - "$WORKDIR/$ARRAY_ID.fasta" "$WORKDIR" << 'PYTHON'
import sys
from collections import Counter

fasta_file = sys.argv[1]
outdir = sys.argv[2]

# Read sequence
seq = ""
with open(fasta_file) as f:
    for line in f:
        if not line.startswith(">"):
            seq += line.strip().upper()

print(f"Sequence length: {len(seq)} bp")
print(f"Expected period: 172 bp")
print(f"Expected N monomers: {len(seq) // 172}")
print()

# For each k: count k-mers, find max frequency
print(f"{'k':>5} {'max_freq':>10} {'n_distinct':>12} {'est_period':>12} {'top_kmer':>25}")
print(f"{'-'*5} {'-'*10} {'-'*12} {'-'*12} {'-'*25}")

results = []
for k in list(range(3, 50)) + list(range(50, 200, 2)) + list(range(200, min(500, len(seq)//2), 5)):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        counts[seq[i:i+k]] += 1

    if not counts:
        break

    most_common = counts.most_common(1)[0]
    max_freq = most_common[1]
    top_kmer = most_common[0]
    n_distinct = len(counts)
    est_period = len(seq) / max_freq

    display_kmer = top_kmer[:20] + "..." if len(top_kmer) > 20 else top_kmer
    print(f"{k:>5} {max_freq:>10} {n_distinct:>12} {est_period:>12.1f} {display_kmer:>25}")

    results.append((k, max_freq, n_distinct, est_period, top_kmer))

# Save full results
with open(f"{outdir}/kmer_spectrum.tsv", "w") as f:
    f.write("k\tmax_freq\tn_distinct\test_period\ttop_kmer\n")
    for k, mf, nd, ep, tk in results:
        f.write(f"{k}\t{mf}\t{nd}\t{ep:.1f}\t{tk}\n")

print(f"\nFull results saved to {outdir}/kmer_spectrum.tsv")
PYTHON

echo ""
echo "Done!"
