#!/usr/bin/env python3
"""T2: is the 98,714 bp 'HOR period' reported on chr17_MATERNAL_43624480_43998424 a real
large-scale unit, or a harmonic of a much shorter real period?

Method: FFT-based composition-corrected autocorrelation excess for ALL lags in O(N log N).
- excess(d) = matches(d)/(N-d) - sum_c f_c^2, where matches(d)=sum_c (indicator_c autocorr)[d].
- A clean tandem repeat of period p shows a HARMONIC COMB: peaks at p, 2p, 3p, ... of similar
  excess. A non-folding detector can then report a large multiple (e.g. ~98 kb) though the real
  period is short. If excess(claimed) is comparable to excess(short p) AND claimed ~ k*p, the large
  period is a harmonic artifact, not an independent large unit.

Reports: strongest period in the monomer band, HOR band, and large band; excess at the claimed
period; and whether the claimed period is an integer multiple of a strong shorter period.
"""
import argparse, os
import numpy as np

MAP = {65:0, 67:1, 71:2, 84:3}

def read_one(path, target):
    name, buf, keep = None, [], False
    with open(path) as fh:
        for line in fh:
            if line[0] == '>':
                if keep:
                    return ''.join(buf)
                name = line[1:].split()[0]
                keep = (name == target) or name.startswith(target)
                buf = []
            elif keep:
                buf.append(line.strip())
    return ''.join(buf) if keep else None

def excess_profile(seq, dmax):
    N = len(seq)
    a = np.frombuffer(seq.encode().upper(), dtype=np.uint8)
    code = np.full(N, 4, np.int8)
    for b, v in MAP.items():
        code[a == b] = v
    # FFT size
    L = 1
    while L < 2 * N:
        L <<= 1
    match = np.zeros(N, dtype=np.float64)
    counts = np.zeros(4)
    for c in range(4):
        x = (code == c).astype(np.float64)
        counts[c] = x.sum()
        X = np.fft.rfft(x, L)
        r = np.fft.irfft(X * np.conj(X), L)[:N]
        match += r
    f = counts / counts.sum()
    base = (f**2).sum()
    d = np.arange(N)
    valid = np.arange(N, 0, -1).astype(np.float64)  # N-d pairs (approx; ambiguous ignored)
    rate = np.zeros(N)
    rate[:dmax+1] = match[:dmax+1] / valid[:dmax+1]
    exc = rate - base
    return exc, base

def band_top(exc, lo, hi, k=6):
    hi = min(hi, len(exc)-1)
    if lo >= hi: return []
    idx = np.arange(lo, hi)
    order = idx[np.argsort(-exc[lo:hi])][:k]
    return [(int(i), float(exc[i])) for i in order]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--array', required=True)
    ap.add_argument('--claimed', type=int, default=98714)
    ap.add_argument('--out', default='results')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    seq = read_one(a.fasta, a.array)
    if not seq:
        print("array not found"); return
    N = len(seq)
    dmax = N // 3
    print(f"array {a.array}  len={N} bp  dmax={dmax}")
    exc, base = excess_profile(seq, dmax)
    print(f"composition baseline sum f^2 = {base:.4f}\n")

    print("=== strongest periods by band ===")
    mono = band_top(exc, 100, 400)
    hor  = band_top(exc, 1000, 8000)
    large= band_top(exc, 40000, dmax)
    print("  monomer band [100-400]:", ", ".join(f"{d}:{e:.3f}" for d,e in mono))
    print("  HOR band    [1k-8k]  :", ", ".join(f"{d}:{e:.3f}" for d,e in hor))
    print("  large band  [40k-max]:", ", ".join(f"{d}:{e:.3f}" for d,e in large))

    p_mono = mono[0][0] if mono else 0
    p_hor  = hor[0][0]  if hor  else 0
    cl = a.claimed
    print(f"\n=== claimed period {cl} ===")
    if cl < len(exc):
        print(f"  excess({cl}) = {exc[cl]:.4f}")
    # harmonic check vs strong short periods
    for label, p in [("monomer", p_mono), ("HOR", p_hor)]:
        if p:
            k = round(cl / p)
            print(f"  claimed / {label} {p} = {cl/p:.2f}  (nearest integer {k}; {label}×{k} = {k*p}, "
                  f"excess there = {exc[k*p]:.4f})" if k*p < len(exc) else f"  claimed/{label} {p} = {cl/p:.2f}")
    # harmonic comb of the strongest HOR period: excess at 1x..Nx
    if p_hor:
        print(f"\n=== harmonic comb of the strongest HOR period {p_hor} bp ===")
        combs = []
        kmax = min(dmax // p_hor, 60)
        for k in range(1, kmax+1):
            d = k * p_hor
            if d < len(exc): combs.append((k, d, exc[d]))
        # show first few + the ones near claimed
        for k, d, e in combs[:8]:
            print(f"  {k}x = {d:>7} bp  excess {e:.3f}")
        near = [c for c in combs if abs(c[1]-cl) < p_hor]
        for k, d, e in near:
            print(f"  --> near claimed: {k}x = {d} bp  excess {e:.3f}  (claimed={cl})")
    np.save(os.path.join(a.out, 'excess_profile.npy'), exc)
    print(f"\nsaved excess profile to {a.out}/excess_profile.npy")

if __name__ == '__main__':
    main()
