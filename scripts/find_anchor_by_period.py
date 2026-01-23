#!/usr/bin/env python3
"""
Find anchor (cut sequence) for a diverged tandem repeat array.

Strategy:
1. Period P is known (from autocorrelation or metadata)
2. Expected copy number N = seq_len / P
3. Find k-mers with frequency ≈ N that are unique within each P-window

A good anchor appears exactly once per monomer.
"""

import sys
from collections import Counter


def read_fasta_single(filepath):
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def find_kmer_positions(seq, kmer):
    """Find all positions of kmer in seq."""
    positions = []
    k = len(kmer)
    for i in range(len(seq) - k + 1):
        if seq[i:i+k] == kmer:
            positions.append(i)
    return positions


def check_unique_in_window(positions, period, seq_len):
    """
    Check if k-mer is unique within each period-sized window.
    Returns: (n_windows_with_exactly_1, n_windows_total, max_per_window)

    For each occurrence, check how many other occurrences fall within ±period/2.
    """
    n = len(positions)
    if n == 0:
        return 0, 0, 0

    max_per_window = 0
    violations = 0

    for i in range(n):
        # Count how many occurrences are within [pos - period/2, pos + period/2]
        count_in_window = 0
        for j in range(n):
            if abs(positions[j] - positions[i]) < period and i != j:
                count_in_window += 1
        if count_in_window > 0:
            violations += 1
        max_per_window = max(max_per_window, count_in_window + 1)

    unique_fraction = (n - violations) / n if n > 0 else 0
    return n - violations, n, max_per_window


def check_spacing_regularity(positions, period):
    """
    Check how regular the spacing between consecutive occurrences is.
    Ideal: all gaps are close to `period`.
    Returns: (mean_gap, cv_gap, n_gaps_near_period)
    """
    if len(positions) < 2:
        return 0, 0, 0

    gaps = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]
    mean_gap = sum(gaps) / len(gaps)

    if mean_gap == 0:
        return 0, 0, 0

    variance = sum((g - mean_gap) ** 2 for g in gaps) / len(gaps)
    std_gap = variance ** 0.5
    cv_gap = std_gap / mean_gap

    # Count gaps within ±20% of period
    tolerance = period * 0.3
    near_period = sum(1 for g in gaps if abs(g - period) <= tolerance)

    return mean_gap, cv_gap, near_period


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.fasta> [period] [min_k] [max_k]")
        sys.exit(1)

    seq = read_fasta_single(sys.argv[1])
    period = int(sys.argv[2]) if len(sys.argv) > 2 else 172
    min_k = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    max_k = int(sys.argv[4]) if len(sys.argv) > 4 else 20

    seq_len = len(seq)
    n_expected = seq_len / period
    freq_min = int(n_expected * 0.5)
    freq_max = int(n_expected * 1.5)

    print(f"Sequence: {seq_len} bp")
    print(f"Period: {period} bp")
    print(f"Expected copy number N: {n_expected:.1f}")
    print(f"Frequency filter: [{freq_min}, {freq_max}]")
    print()

    all_candidates = []

    for k in range(min_k, max_k + 1):
        # Count all k-mers
        counts = Counter()
        for i in range(seq_len - k + 1):
            counts[seq[i:i+k]] += 1

        # Filter by frequency
        candidates = [(kmer, freq) for kmer, freq in counts.items()
                      if freq_min <= freq <= freq_max]

        if not candidates:
            print(f"k={k:>2}: no candidates in frequency range [{freq_min}, {freq_max}]")
            continue

        print(f"k={k:>2}: {len(candidates)} candidates with freq in [{freq_min}, {freq_max}]")

        # For each candidate, check uniqueness and spacing
        for kmer, freq in sorted(candidates, key=lambda x: -x[1]):
            positions = find_kmer_positions(seq, kmer)
            n_unique, n_total, max_per_window = check_unique_in_window(positions, period, seq_len)
            mean_gap, cv_gap, near_period = check_spacing_regularity(positions, period)

            score = (n_unique / n_total) * (near_period / max(n_total - 1, 1))

            all_candidates.append({
                "k": k,
                "kmer": kmer,
                "freq": freq,
                "n_unique": n_unique,
                "max_per_window": max_per_window,
                "mean_gap": mean_gap,
                "cv_gap": cv_gap,
                "near_period": near_period,
                "score": score,
            })

    # Sort by score (uniqueness × regularity)
    all_candidates.sort(key=lambda x: -x["score"])

    print(f"\n{'='*100}")
    print(f"Top 30 anchor candidates (sorted by score = uniqueness × regularity):")
    print(f"{'='*100}")
    print(f"{'k':>3} {'kmer':>15} {'freq':>5} {'unique':>7} {'max/win':>7} {'mean_gap':>8} {'CV':>6} {'near_P':>6} {'score':>6}")
    print(f"{'-'*3} {'-'*15} {'-'*5} {'-'*7} {'-'*7} {'-'*8} {'-'*6} {'-'*6} {'-'*6}")

    for c in all_candidates[:30]:
        print(f"{c['k']:>3} {c['kmer']:>15} {c['freq']:>5} "
              f"{c['n_unique']:>3}/{c['freq']:<3} {c['max_per_window']:>7} "
              f"{c['mean_gap']:>8.1f} {c['cv_gap']:>6.3f} "
              f"{c['near_period']:>3}/{c['freq']-1:<3} {c['score']:>6.3f}")

    # Show positions for top 5
    print(f"\n{'='*100}")
    print(f"Positions of top 5 candidates:")
    print(f"{'='*100}")
    for c in all_candidates[:5]:
        positions = find_kmer_positions(seq, c["kmer"])
        gaps = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]
        print(f"\n  {c['kmer']} (k={c['k']}, freq={c['freq']}, score={c['score']:.3f}):")
        print(f"  Gaps: {gaps[:20]}{'...' if len(gaps) > 20 else ''}")
        print(f"  Positions: {positions[:20]}{'...' if len(positions) > 20 else ''}")


if __name__ == "__main__":
    main()
