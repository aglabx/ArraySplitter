#!/usr/bin/env python3
"""Prototype: FFT autocorrelation + 'comb modulation' rule for fundamental period and HOR order.
Not a replacement implementation -- a feasibility check for a review suggestion."""
import numpy as np, math

def excess_profile(seq, dmax):
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    N = len(a); L = 1
    while L < 2 * N: L <<= 1
    match = np.zeros(N); valid = np.zeros(N, dtype=bool); comp = []
    for c in b'ACGT':
        x = (a == c).astype(np.float64); comp.append(x.sum())
        X = np.fft.rfft(x, L); match += np.fft.irfft(X * np.conj(X), L)[:N]
    comp = np.array(comp); f = comp / comp.sum(); base = float((f ** 2).sum())
    d = np.arange(N); n = (N - d).astype(float)
    ac = np.zeros(N); ac[1:] = match[1:] / n[1:]
    dmax = min(dmax, N // 2)
    return ac[:dmax + 1] - base, base, N

def local_peak(E, d, w):
    lo, hi = max(1, d - w), min(len(E) - 1, d + w)
    if lo > hi: return None, 0.0
    j = lo + int(np.argmax(E[lo:hi + 1])); return j, float(E[j])

def call(seq, min_d=5, dmax=None, floor=0.05):
    N = len(seq); dmax = dmax or N // 3
    E, base, N = excess_profile(seq, dmax)
    if len(E) <= min_d + 1: return None
    dstar = min_d + int(np.argmax(E[min_d:])); Estar = float(E[dstar])
    if Estar < floor: return None
    # ---- fundamental: smallest P such that dstar ~ k*P and P carries a significant local peak + a comb
    P = dstar
    for k in range(min(dstar // min_d, 400), 1, -1):
        cand = int(round(dstar / k))
        if cand < min_d: continue
        w = max(1, int(0.02 * cand))
        p, e = local_peak(E, cand, w)
        if p is None or p < min_d: continue
        se = math.sqrt(0.25 / max(1, N - p))
        if e < max(floor, 6 * se) or e < 0.2 * Estar: continue
        # peak must dominate its neighbourhood (not a shoulder of low-complexity decay)
        lo, hi = max(min_d, p - max(3, p // 3)), min(len(E) - 1, p + max(3, p // 3))
        bg = float(np.median(E[lo:hi + 1]))
        if e - bg < 0.5 * e: continue
        # comb: first few multiples must also be peaks
        K = min(8, (len(E) - 1) // p)
        if K >= 3:
            ok = 0
            for j in range(2, K + 1):
                _, ej = local_peak(E, j * p, max(1, int(0.02 * j * p)))
                if ej >= 0.5 * e: ok += 1
            if ok < 0.6 * (K - 1): continue
        P = p; break
    # ---- HOR order from modulation of comb heights
    K = (len(E) - 1) // P
    h = np.array([local_peak(E, k * P, max(1, int(0.015 * k * P)))[1] for k in range(1, K + 1)])
    best_n, best_score = 1, 0.0
    se = math.sqrt(0.25 / max(1, N - dstar))
    for n in range(2, K // 2 + 1):
        on = h[n - 1::n]; mask = np.ones(K, bool); mask[n - 1::n] = False
        off = h[mask][: len(on) * (n - 1)] if n > 1 else h[mask]
        if len(on) < 2 or len(off) == 0: continue
        score = float(np.median(on[:6]) - np.median(off[: 6 * (n - 1)]))
        if score > best_score * 1.25 or (best_n == 1 and score > 0):
            if score > max(0.03, 8 * se): best_n, best_score = n, score
    hor_p = local_peak(E, best_n * P, max(1, int(0.015 * best_n * P)))[0] if best_n > 1 else P
    return dict(mono=P, hor=hor_p, n=best_n, score=round(best_score, 4), raw_argmax=dstar)
