#!/usr/bin/env python3
"""
Find period of tandem repeat using nucleotide autocorrelation.

For each offset d, compute fraction of positions where seq[i] == seq[i+d].
At d = period, this fraction should peak.
"""

import sys


def read_fasta_single(filepath):
    """Read first sequence from FASTA."""
    seq = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def autocorrelation(seq, d):
    """Fraction of positions where seq[i] == seq[i+d]."""
    matches = 0
    total = len(seq) - d
    for i in range(total):
        if seq[i] == seq[i + d]:
            matches += 1
    return matches / total if total > 0 else 0


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.fasta> [min_d] [max_d]")
        sys.exit(1)

    seq = read_fasta_single(sys.argv[1])
    min_d = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    max_d = int(sys.argv[3]) if len(sys.argv) > 3 else 500

    print(f"Sequence: {len(seq)} bp")
    print(f"Scanning offsets {min_d}..{max_d}")
    print()

    # Random expectation for 4-letter alphabet: 0.25
    # But real sequences have nucleotide bias
    a_frac = seq.count('A') / len(seq)
    c_frac = seq.count('C') / len(seq)
    g_frac = seq.count('G') / len(seq)
    t_frac = seq.count('T') / len(seq)
    random_expect = a_frac**2 + c_frac**2 + g_frac**2 + t_frac**2
    print(f"Composition: A={a_frac:.3f} C={c_frac:.3f} G={g_frac:.3f} T={t_frac:.3f}")
    print(f"Random expectation: {random_expect:.4f}")
    print()

    results = []
    for d in range(min_d, min(max_d + 1, len(seq) // 2)):
        ac = autocorrelation(seq, d)
        results.append((d, ac))

    # Sort by autocorrelation to find peaks
    sorted_results = sorted(results, key=lambda x: -x[1])

    # Print all results
    print(f"{'d':>5} {'autocorr':>10} {'excess':>10} {'bar'}")
    print(f"{'-'*5} {'-'*10} {'-'*10} {'-'*40}")
    for d, ac in results:
        excess = ac - random_expect
        bar_len = int(max(0, excess) * 200)
        bar = "#" * min(bar_len, 50)
        marker = ""
        if d in [172, 344]:  # known period and 2x
            marker = " <-- PERIOD" if d == 172 else " <-- 2x"
        print(f"{d:>5} {ac:>10.4f} {excess:>10.4f} {bar}{marker}")

    # Top 10 peaks
    print(f"\nTop 10 peaks:")
    for d, ac in sorted_results[:10]:
        excess = ac - random_expect
        print(f"  d={d:>4}: autocorr={ac:.4f}, excess={excess:.4f}")

    # Check if multiples of top peak
    if sorted_results:
        best_d = sorted_results[0][0]
        print(f"\nBest period candidate: d={best_d}")
        print(f"Multiples:")
        for mult in range(1, 6):
            target = best_d * mult
            if target < max_d:
                for d, ac in results:
                    if d == target:
                        print(f"  {mult}x = {target}: autocorr={ac:.4f}")
                        break


if __name__ == "__main__":
    main()
