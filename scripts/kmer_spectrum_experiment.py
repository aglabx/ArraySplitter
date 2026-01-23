#!/usr/bin/env python3
"""
K-mer frequency spectrum experiment.

Hypothesis: For a tandem repeat array with period P, the k-mer frequency
spectrum should show a characteristic signal at k ≈ P.

For each array, we compute:
- max_freq(k): maximum frequency of any k-mer for given k
- n_distinct(k): number of distinct k-mers
- n_frequent(k): number of k-mers with frequency >= N/2 (where N ≈ array_len/period)

We look for transitions, plateaus, or peaks that correlate with the monomer period.
"""

import sys
from collections import Counter


def read_fasta(filepath):
    """Read FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_id = ""
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line.upper())
    if current_id:
        sequences.append((current_id, "".join(current_seq)))

    return sequences


def count_kmers(sequence, k):
    """Count all k-mers in sequence."""
    counts = Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        counts[kmer] += 1
    return counts


def compute_kmer_spectrum(sequence, k_values):
    """
    For each k, compute:
    - max_freq: frequency of the most common k-mer
    - n_distinct: number of distinct k-mers observed
    - top3_freq: frequencies of top 3 k-mers
    - top_kmer: the most frequent k-mer sequence
    """
    results = []
    seq_len = len(sequence)

    for k in k_values:
        if k >= seq_len:
            break

        counts = count_kmers(sequence, k)

        if not counts:
            break

        most_common = counts.most_common(3)
        max_freq = most_common[0][1]
        top_kmer = most_common[0][0]
        top3_freq = [x[1] for x in most_common]
        n_distinct = len(counts)

        results.append({
            "k": k,
            "max_freq": max_freq,
            "n_distinct": n_distinct,
            "top3_freq": top3_freq,
            "top_kmer": top_kmer if k <= 20 else top_kmer[:10] + "...",
            "total_kmers": seq_len - k + 1,
        })

    return results


def find_period_signal(results, seq_len):
    """
    Look for the period signal in the k-mer spectrum.

    Key insight: for a repeat of period P with N copies:
    - max_freq ≈ N for k < conserved_core_length
    - estimated_period = seq_len / max_freq should be ~P

    We look for:
    1. Plateaus where max_freq stays constant (ratio ≈ 1.0)
    2. Steps where max_freq drops to a new plateau
    3. The estimated_period = seq_len / max_freq at each plateau
    """
    if len(results) < 3:
        return [], []

    # Detect plateaus: sequences of k values where max_freq is stable
    plateaus = []  # (start_k, end_k, freq, est_period)
    plateau_start = 0
    plateau_freq = results[0]["max_freq"]

    for i in range(1, len(results)):
        curr_freq = results[i]["max_freq"]
        # Plateau continues if freq is within 5% of plateau freq
        if abs(curr_freq - plateau_freq) / max(plateau_freq, 1) <= 0.05:
            continue
        else:
            # End of plateau
            if i - plateau_start >= 2:  # At least 2 consecutive k values
                est_period = seq_len / plateau_freq
                plateaus.append({
                    "start_k": results[plateau_start]["k"],
                    "end_k": results[i-1]["k"],
                    "freq": plateau_freq,
                    "est_period": est_period,
                    "n_points": i - plateau_start,
                })
            plateau_start = i
            plateau_freq = curr_freq

    # Don't forget last plateau
    if len(results) - plateau_start >= 2:
        est_period = seq_len / plateau_freq
        plateaus.append({
            "start_k": results[plateau_start]["k"],
            "end_k": results[-1]["k"],
            "freq": plateau_freq,
            "est_period": est_period,
            "n_points": len(results) - plateau_start,
        })

    # Also compute estimated_period for each k
    period_estimates = []
    for r in results:
        if r["max_freq"] > 0:
            period_estimates.append({
                "k": r["k"],
                "max_freq": r["max_freq"],
                "est_period": seq_len / r["max_freq"],
            })

    return plateaus, period_estimates


def analyze_array(sequence, max_k=300):
    """Full analysis of one array."""
    seq_len = len(sequence)

    # Determine k range
    k_max = min(max_k, seq_len // 2)

    # Use step=1 for small k, larger step for big k
    k_values = list(range(3, min(50, k_max), 1))
    k_values += list(range(50, min(200, k_max), 2))
    k_values += list(range(200, k_max, 5))
    k_values = sorted(set(k_values))

    results = compute_kmer_spectrum(sequence, k_values)
    plateaus, period_estimates = find_period_signal(results, seq_len)

    return results, plateaus, period_estimates


def print_spectrum(header, sequence, results, plateaus, known_period=None):
    """Print the spectrum for one array."""
    seq_len = len(sequence)

    print(f"\n{'='*80}")
    print(f"Array: {header} ({seq_len:,} bp)")
    if known_period:
        expected_n = seq_len // known_period
        print(f"Known period: {known_period} bp, expected N monomers: ~{expected_n}")
    print(f"{'='*80}")

    print(f"\n{'k':>5} {'max_freq':>10} {'n_distinct':>12} {'ratio':>8} {'est_period':>12} {'top_kmer':>22}")
    print(f"{'-'*5} {'-'*10} {'-'*12} {'-'*8} {'-'*12} {'-'*22}")

    prev_freq = None
    for r in results:
        k = r["k"]
        mf = r["max_freq"]
        nd = r["n_distinct"]
        ratio = mf / prev_freq if prev_freq and prev_freq > 0 else 1.0
        top = r["top_kmer"]
        est_p = seq_len / mf if mf > 0 else 0

        # Mark interesting points
        marker = ""
        if known_period and abs(k - known_period) <= 2:
            marker = " <-- PERIOD"
        elif prev_freq and ratio < 0.7:
            marker = " <-- DROP"
        elif prev_freq and ratio >= 0.99:
            marker = " plateau"

        print(f"{k:>5} {mf:>10,} {nd:>12,} {ratio:>8.3f} {est_p:>12.1f} {top:>22}{marker}")
        prev_freq = mf

    # Print plateaus
    if plateaus:
        print(f"\nDetected plateaus (stable max_freq regions):")
        print(f"  {'k_range':>12} {'freq':>8} {'est_period':>12} {'n_points':>8}")
        for p in plateaus:
            period_str = f"{p['est_period']:.1f}"
            match = ""
            if known_period and abs(p['est_period'] - known_period) / known_period < 0.15:
                match = " <-- MATCH"
            print(f"  {p['start_k']:>4}-{p['end_k']:<5} {p['freq']:>8} {period_str:>12} {p['n_points']:>8}{match}")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.fasta> [max_arrays] [max_k]")
        print(f"  Computes k-mer frequency spectrum for each array")
        print(f"  Looks for period signal in the spectrum")
        sys.exit(1)

    input_file = sys.argv[1]
    max_arrays = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    max_k = int(sys.argv[3]) if len(sys.argv) > 3 else 300

    sequences = read_fasta(input_file)
    print(f"Read {len(sequences)} arrays from {input_file}")

    # Parse known periods from headers (format: chrN_start_end_len_period_TYPE)
    for header, sequence in sequences[:max_arrays]:
        parts = header.split("_")
        known_period = None
        if len(parts) >= 6:
            try:
                known_period = int(parts[4])
            except ValueError:
                pass

        results, plateaus, _ = analyze_array(sequence, max_k)
        print_spectrum(header, sequence, results, plateaus, known_period)

    # Summary: for each array, print the best plateau estimate
    print(f"\n\n{'='*80}")
    print(f"SUMMARY: Period detection from k-mer spectrum plateaus")
    print(f"{'='*80}")
    print(f"{'Array':>40} {'SeqLen':>8} {'Known_P':>8} {'1st_plateau':>12} {'est_period':>12} {'Match':>6}")

    for header, sequence in sequences[:max_arrays]:
        parts = header.split("_")
        known_period = None
        if len(parts) >= 6:
            try:
                known_period = int(parts[4])
            except ValueError:
                pass

        _, plateaus, _ = analyze_array(sequence, max_k)

        # Use the first significant plateau (freq > 10, spanning > 3 k values)
        sig_plateaus = [p for p in plateaus if p["freq"] > 10 and p["n_points"] >= 3]
        if sig_plateaus:
            best = sig_plateaus[0]
            est_p = best["est_period"]
            match = "YES" if known_period and abs(est_p - known_period) / known_period < 0.15 else "no"
            plateau_range = f"{best['start_k']}-{best['end_k']}"
        else:
            est_p = 0
            match = "?"
            plateau_range = "none"

        short_name = header[:38] if len(header) > 38 else header
        print(f"{short_name:>40} {len(sequence):>8,} {known_period or '?':>8} {plateau_range:>12} {est_p:>12.1f} {match:>6}")


if __name__ == "__main__":
    main()
