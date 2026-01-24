#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ArraySplitter Refine - Post-processing of decomposed monomers.

This module refines the .monomers.tsv file by:
1. Splitting merged monomers (doubled, tripled, etc.) that have mutations in cut sites
2. Merging fragmented half-monomers with their neighbors

Input/Output: .monomers.tsv format
"""

import argparse
from collections import Counter
from typing import List, Dict, Optional, Tuple
import editdistance as ed


def read_monomers_tsv(filepath: str) -> Tuple[List[str], List[Dict]]:
    """
    Read .monomers.tsv file.
    Returns (header_fields, list of row dicts).
    """
    rows = []
    header_fields = []

    with open(filepath, 'r') as f:
        header = f.readline().strip()
        header_fields = header.split('\t')

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            row = {}
            for i, field in enumerate(header_fields):
                if i < len(parts):
                    row[field] = parts[i]
                else:
                    row[field] = ''

            # Convert numeric fields
            if 'length' in row:
                row['length'] = int(row['length']) if row['length'] else 0
            if 'index' in row:
                row['index'] = int(row['index']) if row['index'] else 0

            rows.append(row)

    return header_fields, rows


def write_monomers_tsv(filepath: str, header_fields: List[str], rows: List[Dict]):
    """Write rows to .monomers.tsv file."""
    with open(filepath, 'w') as f:
        f.write('\t'.join(header_fields) + '\n')
        for row in rows:
            values = []
            for field in header_fields:
                val = row.get(field, '')
                values.append(str(val))
            f.write('\t'.join(values) + '\n')


def find_cut_sequence(rows: List[Dict], kmer_len: int = 16) -> Optional[str]:
    """
    Discover the most common starting k-mer from monomers.
    """
    kmer_counts: Counter = Counter()

    for row in rows:
        if row.get('type') == 'MONOMER':
            seq = row.get('sequence', '')
            if len(seq) >= kmer_len:
                kmer = seq[:kmer_len]
                kmer_counts[kmer] += 1

    if not kmer_counts:
        return None

    return kmer_counts.most_common(1)[0][0]


def get_expected_lengths(rows: List[Dict]) -> Tuple[List[int], int]:
    """
    Analyze length distribution to find expected monomer lengths.
    Returns (common_lengths, modal_length).
    """
    lengths = []
    for row in rows:
        if row.get('type') == 'MONOMER':
            lengths.append(row.get('length', 0))

    if not lengths:
        return [], 700

    length_counts = Counter(lengths)
    common_lengths = [l for l, c in length_counts.most_common(10) if c >= 5]

    if not common_lengths:
        common_lengths = [length_counts.most_common(1)[0][0]]

    modal_length = common_lengths[0]
    return common_lengths, modal_length


def find_approximate_cut_positions(sequence: str, cut_seq: str,
                                    max_mismatches: int = 3) -> List[int]:
    """
    Find positions where cut_seq appears with up to max_mismatches.
    """
    positions = []
    cut_len = len(cut_seq)

    for i in range(len(sequence) - cut_len + 1):
        window = sequence[i:i + cut_len]
        dist = ed.eval(cut_seq, window)
        if dist <= max_mismatches:
            positions.append(i)

    return positions


def split_sequence(sequence: str, cut_seq: str, expected_length: int,
                   max_mismatches: int = 3) -> List[str]:
    """
    Split a merged sequence into individual monomers.
    """
    seq_len = len(sequence)
    num_monomers = round(seq_len / expected_length)

    if num_monomers < 2:
        return [sequence]

    cut_positions = find_approximate_cut_positions(sequence, cut_seq, max_mismatches)

    if not cut_positions:
        return [sequence]

    # Filter positions - need minimum distance between cuts
    min_distance = expected_length * 0.6
    filtered_positions = []
    last_pos = 0

    for pos in sorted(cut_positions):
        if pos > 0 and pos - last_pos >= min_distance:
            filtered_positions.append(pos)
            last_pos = pos

    if not filtered_positions:
        return [sequence]

    # Split at found positions
    result = []
    prev_pos = 0

    for pos in filtered_positions:
        part = sequence[prev_pos:pos]
        if part:
            result.append(part)
        prev_pos = pos

    last_part = sequence[prev_pos:]
    if last_part:
        result.append(last_part)

    return result


def refine_monomers_tsv(input_tsv: str, output_tsv: str,
                        expected_length: Optional[int] = None,
                        max_mismatches: int = 3,
                        min_length_ratio: float = 0.5,
                        max_length_ratio: float = 1.5,
                        verbose: bool = False) -> Dict:
    """
    Main refinement function for TSV files.
    """
    # Read input
    header_fields, rows = read_monomers_tsv(input_tsv)

    if verbose:
        print(f"Read {len(rows)} rows from {input_tsv}")

    # Find cut sequence
    cut_seq = find_cut_sequence(rows)
    if cut_seq is None:
        print("Warning: Could not detect cut sequence")
    if verbose:
        print(f"Detected cut sequence: {cut_seq}")

    # Analyze length distribution
    common_lengths, modal_length = get_expected_lengths(rows)

    if expected_length is None:
        expected_length = modal_length

    if verbose:
        print(f"Expected monomer length: {expected_length} bp")
        print(f"Common lengths: {common_lengths[:5]}")

    # Define thresholds
    min_length = int(expected_length * min_length_ratio)
    max_length = int(expected_length * max_length_ratio)

    if verbose:
        print(f"Normal range: {min_length}-{max_length} bp")

    # Process rows
    refined_rows = []
    stats = {
        'total_input': len(rows),
        'normal': 0,
        'split': 0,
        'merged': 0,
        'flanks': 0,
        'split_into': 0,
        'merged_from': 0,
    }

    i = 0
    while i < len(rows):
        row = rows[i]
        row_type = row.get('type', '')
        seq = row.get('sequence', '')
        seq_len = row.get('length', len(seq))

        # Skip flanks - keep as is
        if 'FLANK' in row_type:
            refined_rows.append(row)
            stats['flanks'] += 1
            i += 1
            continue

        # Normal length - keep as is
        if min_length <= seq_len <= max_length:
            refined_rows.append(row)
            stats['normal'] += 1
            i += 1
            continue

        # Merged monomer (too long) - try to split
        if seq_len > max_length:
            if cut_seq is not None:
                parts = split_sequence(seq, cut_seq, expected_length, max_mismatches)

                if len(parts) > 1:
                    stats['split'] += 1
                    stats['split_into'] += len(parts)

                    if verbose:
                        print(f"  Split {seq_len}bp -> {[len(p) for p in parts]}")

                    # Create new rows for each part
                    for j, part in enumerate(parts):
                        new_row = row.copy()
                        new_row['sequence'] = part
                        new_row['length'] = len(part)
                        new_row['type'] = f"MONOMER|split_{j+1}of{len(parts)}"
                        refined_rows.append(new_row)

                    i += 1
                    continue

            # Could not split
            row_copy = row.copy()
            row_copy['type'] = row_type + '|unsplit'
            refined_rows.append(row_copy)
            i += 1
            continue

        # Short monomer - try to merge with neighbor
        if seq_len < min_length:
            merged = False

            # Look at next row
            if i + 1 < len(rows):
                next_row = rows[i + 1]
                next_type = next_row.get('type', '')
                next_seq = next_row.get('sequence', '')
                next_len = next_row.get('length', len(next_seq))

                # Check if same array and both are short monomers
                if (row.get('sequence_id') == next_row.get('sequence_id') and
                    'FLANK' not in next_type and
                    next_len < min_length):

                    combined_len = seq_len + next_len
                    # Check if combined is reasonable
                    if min_length <= combined_len <= max_length * 1.1:
                        # Merge
                        merged_seq = seq + next_seq
                        merged_row = row.copy()
                        merged_row['sequence'] = merged_seq
                        merged_row['length'] = len(merged_seq)
                        merged_row['type'] = 'MONOMER|merged'
                        refined_rows.append(merged_row)

                        stats['merged'] += 1
                        stats['merged_from'] += 2

                        if verbose:
                            print(f"  Merged {seq_len}bp + {next_len}bp = {len(merged_seq)}bp")

                        i += 2
                        merged = True

            if not merged:
                # Keep as is with tag
                row_copy = row.copy()
                row_copy['type'] = row_type + '|short'
                refined_rows.append(row_copy)
                i += 1
            continue

        # Fallback
        refined_rows.append(row)
        i += 1

    stats['total_output'] = len(refined_rows)

    # Renumber indices within each sequence_id
    current_seq_id = None
    current_idx = 0
    for row in refined_rows:
        seq_id = row.get('sequence_id')
        if seq_id != current_seq_id:
            current_seq_id = seq_id
            current_idx = 0
        row['index'] = current_idx
        current_idx += 1

    # Write output
    write_monomers_tsv(output_tsv, header_fields, refined_rows)

    if verbose:
        print(f"\nRefinement statistics:")
        print(f"  Input rows: {stats['total_input']}")
        print(f"  Output rows: {stats['total_output']}")
        print(f"  Normal (unchanged): {stats['normal']}")
        print(f"  Split merged: {stats['split']} -> {stats['split_into']} parts")
        print(f"  Merged short: {stats['merged']} (from {stats['merged_from']} fragments)")
        print(f"  Flanks (unchanged): {stats['flanks']}")

    return stats


def main(input_tsv: str, output_tsv: str,
         expected_length: Optional[int] = None,
         max_mismatches: int = 3,
         verbose: bool = False) -> Dict:
    """Main entry point for refine command."""

    stats = refine_monomers_tsv(
        input_tsv=input_tsv,
        output_tsv=output_tsv,
        expected_length=expected_length,
        max_mismatches=max_mismatches,
        verbose=verbose
    )

    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Refine ArraySplitter .monomers.tsv by splitting merged monomers and joining fragments"
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input .monomers.tsv file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output .monomers.tsv file (refined)")
    parser.add_argument("-l", "--length", type=int, default=None,
                        help="Expected monomer length (auto-detected if not specified)")
    parser.add_argument("-m", "--mismatches", type=int, default=3,
                        help="Maximum mismatches for cut sequence matching (default: 3)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    main(
        input_tsv=args.input,
        output_tsv=args.output,
        expected_length=args.length,
        max_mismatches=args.mismatches,
        verbose=args.verbose
    )
