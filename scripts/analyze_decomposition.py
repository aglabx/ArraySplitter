#!/usr/bin/env python3
"""
Analyze ArraySplitter decomposition results vs expected monomer sizes.

Parses array names like: chr14_9840269_9854538_14270_171_MONO
Where 171 is the expected monomer size predicted by the upstream tool.

Compares with actual decomposition results to identify discrepancies.
"""

import sys
import re
import argparse
from collections import defaultdict
from pathlib import Path


def parse_array_name(name: str) -> dict:
    """
    Parse array name to extract metadata.
    Format: chr{N}_{start}_{end}_{length}_{expected_monomer}_MONO
    Example: chr14_9840269_9854538_14270_171_MONO
    """
    # Try standard format
    match = re.match(r'^(chr\w+)_(\d+)_(\d+)_(\d+)_(\d+)_MONO', name)
    if match:
        return {
            'chrom': match.group(1),
            'start': int(match.group(2)),
            'end': int(match.group(3)),
            'array_len': int(match.group(4)),
            'expected_monomer': int(match.group(5)),
        }
    return None


def parse_monomers_tsv(filepath: str) -> dict:
    """Parse monomers.tsv file and collect stats per array."""
    arrays = defaultdict(lambda: {
        'lengths': [],
        'eds': [],
        'n_monomers': 0,
        'n_flanks': 0,
        'cut_seq': None,
        'orientation': None,
    })

    with open(filepath) as f:
        header = f.readline().strip().split('\t')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            seq_id = parts[0]
            orientation = parts[1]
            idx = int(parts[2])
            piece_type = parts[3]
            length = int(parts[4])
            ed_str = parts[5]

            arrays[seq_id]['orientation'] = orientation

            if piece_type in ('LEFT_FLANK', 'RIGHT_FLANK'):
                arrays[seq_id]['n_flanks'] += 1
            else:
                arrays[seq_id]['lengths'].append(length)
                arrays[seq_id]['n_monomers'] += 1

                if ed_str != 'N/A' and ed_str != '':
                    try:
                        arrays[seq_id]['eds'].append(int(ed_str))
                    except ValueError:
                        pass

    return dict(arrays)


def parse_decomposed_fasta(filepath: str) -> dict:
    """Parse decomposed.fasta to get cut sequences."""
    arrays = {}

    with open(filepath) as f:
        for line in f:
            if line.startswith('>'):
                # Parse header like: >chr14_... cut=GTTTCT orientation=fwd n_monomers=27 range=169-173 avg=171.0
                parts = line[1:].strip().split()
                array_id = parts[0]

                cut_seq = None
                for part in parts:
                    if part.startswith('cut='):
                        cut_seq = part.split('=')[1]
                        break

                arrays[array_id] = {'cut_seq': cut_seq}

    return arrays


def calculate_stats(lengths: list) -> dict:
    """Calculate statistics for a list of lengths."""
    if not lengths:
        return {'n': 0, 'median': 0, 'mean': 0, 'min': 0, 'max': 0, 'cv': float('inf'), 'rcv': float('inf')}

    n = len(lengths)
    sorted_lens = sorted(lengths)
    median = sorted_lens[n // 2]
    mean = sum(lengths) / n
    min_len = min(lengths)
    max_len = max(lengths)

    # Standard CV
    if mean > 0:
        variance = sum((x - mean) ** 2 for x in lengths) / n
        std = variance ** 0.5
        cv = std / mean
    else:
        cv = float('inf')

    # Robust CV (MAD / median)
    if median > 0:
        deviations = sorted(abs(x - median) for x in lengths)
        mad = deviations[n // 2]
        rcv = mad / median
    else:
        rcv = float('inf')

    return {
        'n': n,
        'median': median,
        'mean': mean,
        'min': min_len,
        'max': max_len,
        'cv': cv,
        'rcv': rcv,
    }


def analyze_results(monomers_file: str, fasta_file: str, output_file: str = None):
    """Main analysis function."""

    # Parse files
    print(f"Parsing {monomers_file}...", file=sys.stderr)
    arrays_data = parse_monomers_tsv(monomers_file)

    print(f"Parsing {fasta_file}...", file=sys.stderr)
    fasta_data = parse_decomposed_fasta(fasta_file)

    # Merge data
    for array_id, fasta_info in fasta_data.items():
        if array_id in arrays_data:
            arrays_data[array_id]['cut_seq'] = fasta_info.get('cut_seq')

    # Analyze each array
    results = []

    for array_id, data in arrays_data.items():
        meta = parse_array_name(array_id)
        if not meta:
            continue

        stats = calculate_stats(data['lengths'])

        expected = meta['expected_monomer']
        actual_median = stats['median']
        actual_mean = stats['mean']

        # Calculate deviation
        if expected > 0:
            deviation_pct = abs(actual_median - expected) / expected * 100
            ratio = actual_median / expected
        else:
            deviation_pct = 0
            ratio = 1.0

        # Classify match quality
        if deviation_pct < 5:
            match_quality = 'EXACT'
        elif deviation_pct < 15:
            match_quality = 'CLOSE'
        elif 0.45 < ratio < 0.55 or 1.9 < ratio < 2.1:
            match_quality = 'HALF/DOUBLE'
        elif 0.3 < ratio < 0.35 or 2.9 < ratio < 3.1:
            match_quality = 'THIRD/TRIPLE'
        else:
            match_quality = 'MISMATCH'

        # Mean ED per bp
        mean_ed = sum(data['eds']) / len(data['eds']) if data['eds'] else 0
        ed_per_bp = mean_ed / actual_median if actual_median > 0 else 0

        results.append({
            'array_id': array_id,
            'chrom': meta['chrom'],
            'array_len': meta['array_len'],
            'expected_monomer': expected,
            'actual_median': actual_median,
            'actual_mean': round(actual_mean, 1),
            'n_monomers': stats['n'],
            'range': f"{stats['min']}-{stats['max']}",
            'cv': round(stats['cv'], 3),
            'rcv': round(stats['rcv'], 3),
            'mean_ed': round(mean_ed, 1),
            'ed_per_bp': round(ed_per_bp, 3),
            'deviation_pct': round(deviation_pct, 1),
            'ratio': round(ratio, 2),
            'match': match_quality,
            'cut_seq': data.get('cut_seq', ''),
            'orientation': data.get('orientation', ''),
        })

    # Sort by deviation
    results.sort(key=lambda x: -x['deviation_pct'])

    # Output
    out = open(output_file, 'w') if output_file else sys.stdout

    # Header
    cols = ['array_id', 'array_len', 'expected', 'actual', 'n_mono', 'range',
            'rcv', 'mean_ed', 'ed/bp', 'dev%', 'ratio', 'match', 'cut_seq']
    print('\t'.join(cols), file=out)

    # Data rows
    for r in results:
        row = [
            r['array_id'],
            str(r['array_len']),
            str(r['expected_monomer']),
            str(r['actual_median']),
            str(r['n_monomers']),
            r['range'],
            f"{r['rcv']:.3f}",
            f"{r['mean_ed']:.1f}",
            f"{r['ed_per_bp']:.3f}",
            f"{r['deviation_pct']:.1f}",
            f"{r['ratio']:.2f}",
            r['match'],
            r['cut_seq'] or '',
        ]
        print('\t'.join(row), file=out)

    if output_file:
        out.close()

    # Summary stats
    print(f"\n=== Summary ===", file=sys.stderr)
    print(f"Total arrays: {len(results)}", file=sys.stderr)

    by_match = defaultdict(int)
    for r in results:
        by_match[r['match']] += 1

    for match_type in ['EXACT', 'CLOSE', 'HALF/DOUBLE', 'THIRD/TRIPLE', 'MISMATCH']:
        count = by_match.get(match_type, 0)
        pct = count / len(results) * 100 if results else 0
        print(f"  {match_type:12s}: {count:4d} ({pct:5.1f}%)", file=sys.stderr)

    # Worst mismatches
    mismatches = [r for r in results if r['match'] == 'MISMATCH']
    if mismatches:
        print(f"\n=== Top 10 Mismatches ===", file=sys.stderr)
        print(f"{'Array':<50s} {'Exp':>5s} {'Act':>5s} {'Ratio':>6s} {'RCV':>6s}", file=sys.stderr)
        for r in mismatches[:10]:
            print(f"{r['array_id']:<50s} {r['expected_monomer']:>5d} {r['actual_median']:>5d} {r['ratio']:>6.2f} {r['rcv']:>6.3f}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Analyze ArraySplitter decomposition results')
    parser.add_argument('-m', '--monomers', required=True, help='Input .monomers.tsv file')
    parser.add_argument('-f', '--fasta', required=True, help='Input .decomposed.fasta file')
    parser.add_argument('-o', '--output', help='Output TSV file (default: stdout)')

    args = parser.parse_args()

    analyze_results(args.monomers, args.fasta, args.output)


if __name__ == '__main__':
    main()
