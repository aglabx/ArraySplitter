#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ArraySplitter Monomer Rotation - Rotate monomers to common starting position.

This module finds the optimal rotation k-mer and rotates all monomers to start
from that position. Important for tandem repeats: the sequence is duplicated
before searching to handle cases where the k-mer spans the boundary.

Algorithm:
1. Find k-mer that appears exactly ONCE in the maximum number of monomers
2. Duplicate each sequence (seq+seq) to handle boundary cases
3. Find best match position using edit distance
4. Rotate sequence to start from that position
5. Ensure all sequences are on the same strand
"""

import argparse
from collections import Counter
from typing import List, Tuple, Dict, Optional
import editdistance as ed


def get_revcomp(seq: str) -> str:
    """Get reverse complement of sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(b, b) for b in reversed(seq))


def find_best_rotation_kmer(sequences: List[str], k: int = 16) -> Tuple[Optional[str], float]:
    """
    Find k-mer that appears exactly once in the maximum number of sequences.

    Args:
        sequences: List of DNA sequences (monomers)
        k: K-mer length (default 16)

    Returns:
        (best_kmer, fraction_unique) - k-mer and fraction of sequences with unique occurrence
    """
    kmer_unique_count: Counter = Counter()

    for seq in sequences:
        # Count k-mers in this sequence
        seq_kmers: Counter = Counter()
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            seq_kmers[kmer] += 1

        # Record k-mers that appear exactly once
        for kmer, count in seq_kmers.items():
            if count == 1:
                kmer_unique_count[kmer] += 1

    if not kmer_unique_count:
        return None, 0.0

    best_kmer, best_count = kmer_unique_count.most_common(1)[0]
    fraction = best_count / len(sequences)

    return best_kmer, fraction


def find_kmer_position_in_doubled(seq: str, kmer: str, max_ed: int = 5) -> Tuple[int, int, str]:
    """
    Find best position of k-mer in doubled sequence (seq+seq) to handle boundary cases.
    Checks both strands.

    Args:
        seq: Original sequence
        kmer: K-mer to find
        max_ed: Maximum edit distance allowed

    Returns:
        (position, edit_distance, strand) where position is in original [0, len(seq))
    """
    k = len(kmer)
    seq_len = len(seq)

    # Double the sequence to handle boundary
    doubled_fwd = seq + seq
    doubled_rev = get_revcomp(seq) + get_revcomp(seq)

    best_pos = -1
    best_dist = max_ed + 1
    best_strand = "unk"

    # Search in forward strand (doubled)
    for i in range(seq_len):  # Only search first copy positions
        window = doubled_fwd[i:i+k]
        dist = ed.eval(kmer, window)
        if dist < best_dist:
            best_dist = dist
            best_pos = i
            best_strand = "fwd"

    # Search in reverse complement (doubled)
    for i in range(seq_len):
        window = doubled_rev[i:i+k]
        dist = ed.eval(kmer, window)
        if dist < best_dist:
            best_dist = dist
            best_pos = i
            best_strand = "rev"

    if best_dist <= max_ed:
        return best_pos, best_dist, best_strand

    return -1, best_dist, "unk"


def rotate_sequence(seq: str, position: int, strand: str) -> str:
    """
    Rotate sequence to start from given position.

    Args:
        seq: Original sequence
        position: Position to start from
        strand: "fwd" or "rev" - which strand to use

    Returns:
        Rotated sequence
    """
    if strand == "rev":
        seq = get_revcomp(seq)

    return seq[position:] + seq[:position]


def rotate_monomers(input_fasta: str, output_fasta: str,
                    kmer: Optional[str] = None, k: int = 16,
                    max_ed: int = 5, verbose: bool = False) -> Dict:
    """
    Main rotation function.

    Args:
        input_fasta: Input FASTA file with monomers
        output_fasta: Output FASTA file with rotated monomers
        kmer: Optional rotation k-mer (auto-detected if not provided)
        k: K-mer length for auto-detection
        max_ed: Maximum edit distance for rotation
        verbose: Print progress

    Returns:
        Statistics dictionary
    """
    # Read sequences
    sequences = []
    headers = []

    with open(input_fasta) as f:
        seq = ""
        header = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    headers.append(header)
                seq = ""
                header = line.strip()
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
            headers.append(header)

    if verbose:
        print(f"Loaded {len(sequences)} sequences")

    # Auto-detect rotation k-mer if not provided
    if kmer is None:
        # Filter to normal-length monomers for k-mer detection
        normal_seqs = [s for s in sequences if 600 <= len(s) <= 800]
        if not normal_seqs:
            normal_seqs = sequences

        kmer, fraction = find_best_rotation_kmer(normal_seqs, k)

        if verbose:
            print(f"Auto-detected rotation k-mer: {kmer}")
            print(f"  Unique in {fraction*100:.1f}% of monomers")

    if kmer is None:
        raise ValueError("Could not detect rotation k-mer")

    # Rotate all sequences
    rotated = []
    stats = {
        "total": len(sequences),
        "fwd": 0,
        "rev": 0,
        "failed": 0,
        "ed_distribution": Counter()
    }

    for header, seq in zip(headers, sequences):
        pos, dist, strand = find_kmer_position_in_doubled(seq, kmer, max_ed)

        if pos >= 0:
            rotated_seq = rotate_sequence(seq, pos, strand)
            stats[strand] += 1
            stats["ed_distribution"][dist] += 1
        else:
            # Keep original if rotation failed
            rotated_seq = seq
            strand = "failed"
            stats["failed"] += 1

        rotated.append((header, rotated_seq, strand, dist))

    # Write output
    with open(output_fasta, 'w') as f:
        for header, seq, strand, dist in rotated:
            new_header = f"{header}|rot_{strand}|ed{dist}"
            f.write(f"{new_header}\n{seq}\n")

    if verbose:
        print(f"\nRotation results:")
        print(f"  Forward: {stats['fwd']}")
        print(f"  Reverse: {stats['rev']}")
        print(f"  Failed: {stats['failed']}")
        print(f"\nEdit distance distribution:")
        for d, c in sorted(stats["ed_distribution"].items()):
            print(f"  ED={d}: {c}")

    stats["kmer"] = kmer
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Rotate monomers to common starting position"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file with monomers")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA file with rotated monomers")
    parser.add_argument("-k", "--kmer", default=None,
                        help="Rotation k-mer (auto-detected if not provided)")
    parser.add_argument("-l", "--kmer-length", type=int, default=16,
                        help="K-mer length for auto-detection (default: 16)")
    parser.add_argument("-e", "--max-ed", type=int, default=5,
                        help="Maximum edit distance for rotation (default: 5)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    stats = rotate_monomers(
        input_fasta=args.input,
        output_fasta=args.output,
        kmer=args.kmer,
        k=args.kmer_length,
        max_ed=args.max_ed,
        verbose=args.verbose
    )

    print(f"Rotation complete. K-mer used: {stats['kmer']}")


if __name__ == "__main__":
    main()
