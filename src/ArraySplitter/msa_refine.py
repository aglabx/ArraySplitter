#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MSA Refine - Remove sequences with rare indels from multiple sequence alignment.

This module analyzes indel patterns in MSA and removes sequences containing
rare indels (below frequency threshold). This helps clean up alignments
before final consensus generation.

Algorithm:
1. Parse aligned FASTA and identify all indel patterns (contiguous gaps)
2. Calculate frequency of each indel pattern across all sequences
3. Mark sequences containing rare indels (frequency < threshold)
4. Output: cleaned sequences, removed sequences, statistics report
"""

import argparse
from collections import Counter, defaultdict
from typing import List, Tuple, Dict, Set
from dataclasses import dataclass
import re


@dataclass
class Indel:
    """Represents an indel (gap region) in the alignment."""
    start: int
    end: int
    length: int

    def __hash__(self):
        return hash((self.start, self.end))

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end


@dataclass
class IndelStats:
    """Statistics for an indel pattern."""
    indel: Indel
    count: int
    frequency: float
    sequences: List[str]  # sequence names containing this indel


def parse_fasta(filepath: str) -> List[Tuple[str, str]]:
    """Parse FASTA file, return list of (name, sequence) tuples."""
    sequences = []
    with open(filepath) as f:
        name = ""
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append((name, seq))
                name = line[1:].split()[0]  # Take first word as name
                seq = ""
            else:
                seq += line
        if seq:
            sequences.append((name, seq))
    return sequences


def find_indels_in_sequence(seq: str) -> List[Indel]:
    """Find all indel (gap) regions in a sequence."""
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

    # Handle gap at end
    if in_gap:
        indels.append(Indel(start=gap_start, end=len(seq), length=len(seq) - gap_start))

    return indels


def analyze_indels(sequences: List[Tuple[str, str]]) -> Dict[Indel, IndelStats]:
    """Analyze all indels across sequences and compute statistics."""
    indel_to_seqs: Dict[Indel, List[str]] = defaultdict(list)
    n_seqs = len(sequences)

    for name, seq in sequences:
        indels = find_indels_in_sequence(seq)
        for indel in indels:
            indel_to_seqs[indel].append(name)

    # Compute statistics
    stats = {}
    for indel, seq_names in indel_to_seqs.items():
        count = len(seq_names)
        stats[indel] = IndelStats(
            indel=indel,
            count=count,
            frequency=count / n_seqs,
            sequences=seq_names
        )

    return stats


def refine_msa(
    input_fasta: str,
    output_clean: str,
    output_removed: str,
    output_stats: str,
    min_freq: float = 0.05,
    verbose: bool = False
) -> Dict:
    """
    Refine MSA by removing sequences with rare indels.

    Args:
        input_fasta: Input aligned FASTA file
        output_clean: Output FASTA with cleaned sequences
        output_removed: Output FASTA with removed sequences
        output_stats: Output TSV with indel statistics
        min_freq: Minimum frequency threshold for indels (default 5%)
        verbose: Print progress

    Returns:
        Statistics dictionary
    """
    # Parse sequences
    sequences = parse_fasta(input_fasta)
    n_total = len(sequences)

    if verbose:
        print(f"Loaded {n_total} aligned sequences")

    # Analyze indels
    indel_stats = analyze_indels(sequences)

    if verbose:
        print(f"Found {len(indel_stats)} unique indel patterns")

    # Find rare indels and sequences to remove
    rare_indels = {indel: stats for indel, stats in indel_stats.items()
                   if stats.frequency < min_freq}

    # Collect sequences to remove (those with ANY rare indel)
    seqs_to_remove: Set[str] = set()
    for stats in rare_indels.values():
        seqs_to_remove.update(stats.sequences)

    if verbose:
        print(f"Found {len(rare_indels)} rare indels (freq < {min_freq*100:.1f}%)")
        print(f"Sequences to remove: {len(seqs_to_remove)}")

    # Split sequences
    clean_seqs = []
    removed_seqs = []

    for name, seq in sequences:
        if name in seqs_to_remove:
            removed_seqs.append((name, seq))
        else:
            clean_seqs.append((name, seq))

    # Write clean sequences (remove columns that are now all gaps)
    if clean_seqs:
        # Find columns that have at least one non-gap
        align_len = len(clean_seqs[0][1])
        non_gap_cols = []
        for i in range(align_len):
            has_base = any(seq[i] != '-' for _, seq in clean_seqs)
            if has_base:
                non_gap_cols.append(i)

        # Write cleaned sequences with gap-only columns removed
        with open(output_clean, 'w') as f:
            for name, seq in clean_seqs:
                cleaned_seq = ''.join(seq[i] for i in non_gap_cols)
                f.write(f">{name}\n{cleaned_seq}\n")

        cols_removed = align_len - len(non_gap_cols)
        if verbose and cols_removed > 0:
            print(f"Removed {cols_removed} gap-only columns")

    # Write removed sequences
    with open(output_removed, 'w') as f:
        for name, seq in removed_seqs:
            f.write(f">{name}\n{seq}\n")

    # Write statistics
    with open(output_stats, 'w') as f:
        f.write("# MSA Refinement Statistics\n")
        f.write(f"# Input: {input_fasta}\n")
        f.write(f"# Total sequences: {n_total}\n")
        f.write(f"# Clean sequences: {len(clean_seqs)}\n")
        f.write(f"# Removed sequences: {len(removed_seqs)}\n")
        f.write(f"# Min frequency threshold: {min_freq*100:.1f}%\n")
        f.write(f"# Unique indel patterns: {len(indel_stats)}\n")
        f.write(f"# Rare indels (below threshold): {len(rare_indels)}\n")
        f.write("#\n")
        f.write("# Indel details (sorted by frequency):\n")
        f.write("start\tend\tlength\tcount\tfrequency\tis_rare\texample_sequences\n")

        # Sort by frequency descending
        sorted_indels = sorted(indel_stats.items(),
                               key=lambda x: (-x[1].frequency, x[0].start))

        for indel, stats in sorted_indels:
            is_rare = "RARE" if stats.frequency < min_freq else ""
            example_seqs = ",".join(stats.sequences[:3])
            if len(stats.sequences) > 3:
                example_seqs += f"... (+{len(stats.sequences)-3} more)"
            f.write(f"{indel.start}\t{indel.end}\t{indel.length}\t"
                    f"{stats.count}\t{stats.frequency*100:.2f}%\t{is_rare}\t{example_seqs}\n")

    if verbose:
        print(f"\nOutput files:")
        print(f"  Clean sequences: {output_clean}")
        print(f"  Removed sequences: {output_removed}")
        print(f"  Statistics: {output_stats}")

    return {
        "total": n_total,
        "clean": len(clean_seqs),
        "removed": len(removed_seqs),
        "total_indels": len(indel_stats),
        "rare_indels": len(rare_indels),
        "min_freq": min_freq
    }


def main():
    parser = argparse.ArgumentParser(
        description="Refine MSA by removing sequences with rare indels"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input aligned FASTA file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix for cleaned files")
    parser.add_argument("-f", "--min-freq", type=float, default=0.05,
                        help="Minimum indel frequency threshold (default: 0.05 = 5%%)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    output_clean = f"{args.output}.clean.fasta"
    output_removed = f"{args.output}.removed.fasta"
    output_stats = f"{args.output}.indel_stats.tsv"

    stats = refine_msa(
        input_fasta=args.input,
        output_clean=output_clean,
        output_removed=output_removed,
        output_stats=output_stats,
        min_freq=args.min_freq,
        verbose=args.verbose
    )

    print(f"\nRefinement complete:")
    print(f"  Total: {stats['total']} → Clean: {stats['clean']}, Removed: {stats['removed']}")
    print(f"  Indels: {stats['total_indels']} total, {stats['rare_indels']} rare (< {stats['min_freq']*100:.1f}%)")


if __name__ == "__main__":
    main()
