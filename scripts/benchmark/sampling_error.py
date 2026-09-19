"""Re-implement v1.8.1 generate_sample_positions (SplitMix64 stratified jitter) and measure
|sampled - exhaustive| autocorrelation over all scanned lags (5..5000) on >50 kb arrays."""
import sys, numpy as np
M = (1 << 64) - 1
def positions(L, n=10000):
    st = 0x9e3779b97f4a7c15; out = []
    for i in range(n):
        st = (st + 0x9e3779b97f4a7c15) & M
        z = st
        z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & M
        z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & M
        r = z ^ (z >> 31)
        b0 = i * L // n; b1 = (i + 1) * L // n; bl = b1 - b0
        out.append(b0 + (r % bl if bl > 0 else 0))
    return np.array(out)
fa = {}; name = None
for l in open(sys.argv[1]):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; fa[name] = []
    else: fa[name].append(l.upper())
print(f"{'array':36s} {'L':>8s} {'max|err|':>9s} {'95th pct':>9s} {'lags with |err|>0.003':>22s}")
for k, v in fa.items():
    s = ''.join(v); L = len(s)
    if L <= 50000: continue
    a = np.frombuffer(s.encode(), dtype=np.uint8); pos = positions(L)
    errs = []
    for d in range(5, 5001):
        full = float((a[:-d] == a[d:]).mean())
        p = pos[pos + d < L]; samp = float((a[p] == a[p + d]).mean())
        errs.append(abs(samp - full))
    errs = np.array(errs)
    print(f"{k[:36]:36s} {L:8d} {errs.max():9.4f} {np.percentile(errs,95):9.4f} {int((errs>0.003).sum()):>10d} of {len(errs)}")
