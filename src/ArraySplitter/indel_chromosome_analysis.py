#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Indel Chromosome Specificity Analysis.

Analyzes how chromosome-specific indels are in MSA alignments.
For each indel, calculates:
- Which chromosomes contain sequences with this indel
- What fraction of sequences from each chromosome have the indel
- Overall specificity score
"""

import argparse
from collections import defaultdict
from typing import List, Tuple, Dict, Set
from dataclasses import dataclass
import re
import json


@dataclass
class Indel:
    start: int
    end: int
    length: int

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
    """Parse FASTA file."""
    sequences = []
    with open(filepath) as f:
        name = ""
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append((name, seq))
                name = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        if seq:
            sequences.append((name, seq))
    return sequences


def extract_chromosome(name: str) -> str:
    """Extract chromosome name from sequence header."""
    # Pattern: chr1_mat:123-456|... or chr1A_pat:...
    match = re.match(r'(chr[^:]+)', name)
    if match:
        return match.group(1)
    return "unknown"


def find_indels_in_sequence(seq: str) -> List[Indel]:
    """Find all indel regions in a sequence."""
    indels = []
    in_gap = False
    gap_start = 0

    for i, base in enumerate(seq):
        if base == '-':
            if not in_gap:
                in_gap = True
                gap_start = i
        else:
            if in_gap:
                indels.append(Indel(start=gap_start, end=i, length=i - gap_start))
                in_gap = False

    if in_gap:
        indels.append(Indel(start=gap_start, end=len(seq), length=len(seq) - gap_start))

    return indels


def analyze_chromosome_specificity(input_fasta: str, output_tsv: str, output_json: str,
                                    min_count: int = 3, verbose: bool = False) -> Dict:
    """
    Analyze chromosome specificity of indels.

    Args:
        input_fasta: Aligned FASTA file
        output_tsv: Output TSV with detailed analysis
        output_json: Output JSON for HTML report
        min_count: Minimum sequence count for indel to be included
        verbose: Print progress
    """
    sequences = parse_fasta(input_fasta)
    n_seqs = len(sequences)

    if verbose:
        print(f"Loaded {n_seqs} sequences")

    # Count sequences per chromosome
    chr_seq_counts: Dict[str, int] = defaultdict(int)
    for name, _ in sequences:
        chrom = extract_chromosome(name)
        chr_seq_counts[chrom] += 1

    if verbose:
        print(f"Found {len(chr_seq_counts)} chromosomes")

    # Map indels to chromosomes and sequences
    indel_to_chr_seqs: Dict[Indel, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))

    for name, seq in sequences:
        chrom = extract_chromosome(name)
        indels = find_indels_in_sequence(seq)
        for indel in indels:
            indel_to_chr_seqs[indel][chrom].append(name)

    # Analyze each indel
    results = []

    for indel, chr_seqs in indel_to_chr_seqs.items():
        total_count = sum(len(seqs) for seqs in chr_seqs.values())

        # Skip rare indels
        if total_count < min_count:
            continue

        # Calculate per-chromosome statistics
        chr_stats = []
        for chrom, seq_names in chr_seqs.items():
            count_in_chr = len(seq_names)
            total_in_chr = chr_seq_counts[chrom]
            fraction = count_in_chr / total_in_chr if total_in_chr > 0 else 0
            chr_stats.append({
                "chr": chrom,
                "count": count_in_chr,
                "total": total_in_chr,
                "fraction": fraction
            })

        # Sort by fraction descending
        chr_stats.sort(key=lambda x: -x["fraction"])

        # Calculate specificity metrics
        n_chromosomes = len(chr_stats)
        max_fraction = chr_stats[0]["fraction"] if chr_stats else 0

        # Specificity score: high if indel is concentrated in few chromosomes
        # If indel appears in all sequences of 1-2 chromosomes, it's highly specific
        top_chr_contribution = sum(s["count"] for s in chr_stats[:2]) / total_count if total_count > 0 else 0

        # Is this indel chromosome-specific?
        # Criteria: >80% of occurrences from top 2 chromosomes AND top chr has >50% of its sequences with indel
        is_chr_specific = top_chr_contribution > 0.8 and max_fraction > 0.5

        results.append({
            "start": indel.start,
            "end": indel.end,
            "length": indel.length,
            "total_count": total_count,
            "frequency": total_count / n_seqs,
            "n_chromosomes": n_chromosomes,
            "max_fraction": max_fraction,
            "top_chr_contribution": top_chr_contribution,
            "is_chr_specific": is_chr_specific,
            "chr_stats": chr_stats
        })

    # Sort by total count descending
    results.sort(key=lambda x: -x["total_count"])

    # Summary statistics
    n_total_indels = len(results)
    n_chr_specific = sum(1 for r in results if r["is_chr_specific"])

    # Write TSV
    with open(output_tsv, 'w') as f:
        f.write("# Indel Chromosome Specificity Analysis\n")
        f.write(f"# Input: {input_fasta}\n")
        f.write(f"# Total sequences: {n_seqs}\n")
        f.write(f"# Chromosomes: {len(chr_seq_counts)}\n")
        f.write(f"# Indels analyzed (count >= {min_count}): {n_total_indels}\n")
        f.write(f"# Chromosome-specific indels: {n_chr_specific} ({n_chr_specific/n_total_indels*100:.1f}%)\n")
        f.write("#\n")
        f.write("start\tend\tlength\ttotal_count\tfrequency\tn_chromosomes\tmax_chr_fraction\ttop2_contribution\tis_specific\ttop_chromosomes\n")

        for r in results:
            top_chrs = ", ".join(f"{s['chr']}({s['count']}/{s['total']}={s['fraction']*100:.0f}%)"
                                  for s in r["chr_stats"][:5])
            f.write(f"{r['start']}\t{r['end']}\t{r['length']}\t{r['total_count']}\t"
                    f"{r['frequency']*100:.2f}%\t{r['n_chromosomes']}\t{r['max_fraction']*100:.1f}%\t"
                    f"{r['top_chr_contribution']*100:.1f}%\t{'YES' if r['is_chr_specific'] else 'no'}\t{top_chrs}\n")

    # Write JSON for HTML report
    json_data = {
        "summary": {
            "total_sequences": n_seqs,
            "n_chromosomes": len(chr_seq_counts),
            "n_indels_analyzed": n_total_indels,
            "n_chr_specific": n_chr_specific,
            "pct_chr_specific": round(n_chr_specific / n_total_indels * 100, 1) if n_total_indels > 0 else 0
        },
        "chromosome_counts": dict(chr_seq_counts),
        "indels": results[:50]  # Top 50 for HTML
    }

    with open(output_json, 'w') as f:
        json.dump(json_data, f, indent=2)

    if verbose:
        print(f"\nResults:")
        print(f"  Indels analyzed: {n_total_indels}")
        print(f"  Chromosome-specific: {n_chr_specific} ({n_chr_specific/n_total_indels*100:.1f}%)")
        print(f"\nOutput: {output_tsv}, {output_json}")

    return json_data


def main():
    parser = argparse.ArgumentParser(description="Analyze chromosome specificity of indels")
    parser.add_argument("-i", "--input", required=True, help="Aligned FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output prefix")
    parser.add_argument("-m", "--min-count", type=int, default=3,
                        help="Minimum count for indel to be included (default: 3)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    analyze_chromosome_specificity(
        input_fasta=args.input,
        output_tsv=f"{args.output}.chr_specificity.tsv",
        output_json=f"{args.output}.chr_specificity.json",
        min_count=args.min_count,
        verbose=args.verbose
    )


if __name__ == "__main__":
    main()
