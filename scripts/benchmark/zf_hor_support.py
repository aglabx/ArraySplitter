import sys, csv, collections, math
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from comb import excess_profile, local_peak
csv.field_size_limit(10**9)
PANEL, PREFIX = sys.argv[1], sys.argv[2]   # zf.fa, output prefix of an arraysplitter run
fa = {}; name = None
for l in open(PANEL):
    l = l.strip()
    if l.startswith('>'): name = l[1:].split()[0]; fa[name] = []
    else: fa[name].append(l.upper())
fa = {k: ''.join(v) for k, v in fa.items()}
rows = list(csv.DictReader(open(PREFIX + '.summary.tsv'), delimiter='\t'))
hor = [r for r in rows if int(r['mono_period']) > 0 and int(r['hor_period']) >= 1.5 * int(r['mono_period'])]
out = []
for r in hor:
    s = fa[r['array_id']]; L = len(s); P = int(r['mono_period']); H = int(r['hor_period'])
    if r['orientation'] == 'rev': s = s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    dmax = min(L // 3, max(3 * H, 40 * P), 60000)
    E, base, N = excess_profile(s, dmax)
    # exact excess at reported HOR lag and at the base-monomer lag (whole array)
    eH = local_peak(E, H, max(1, int(0.01 * H)))[1]
    eP = local_peak(E, P, max(1, int(0.02 * P)))[1] if P < len(E) else float('nan')
    ratio = H / P; n = round(ratio); is_int = abs(ratio - n) < 0.08 * max(1, n) ** 0.5 and n >= 2
    # comb heights at multiples of H itself: is 2H, 3H as high (harmonic family of H)? and at sub-multiples H/j?
    subs = []
    for j in (2, 3, 4, 5, 6, 7, 8):
        d = int(round(H / j))
        if d >= 5: subs.append((j, local_peak(E, d, max(1, int(0.02 * d)))[1]))
    best_sub = max(subs, key=lambda x: x[1]) if subs else (0, 0.0)
    # neighbours among multiples of P: k = n-1, n+1 (if P is meaningful)
    neigh = []
    if P >= 5 and n >= 2:
        for k in (n - 1, n + 1):
            d = k * P
            if 5 <= d < len(E): neigh.append(local_peak(E, d, max(1, int(0.015 * d)))[1])
    out.append(dict(id=r['array_id'], L=L, P=P, H=H, ratio=ratio, eH=eH, eP=eP, best_sub=best_sub, neigh=neigh, ncopies=L / H))
# classify
def cls(o):
    if o['best_sub'][1] >= o['eH'] - 0.01: return 'HARMONIC-LIKE: a shorter period H/j is as strong as H'
    if o['neigh'] and o['eH'] - max(o['neigh']) > 0.03: return 'SUPPORTED: H stands out vs neighbouring multiples of mono period'
    if o['eH'] - o['best_sub'][1] > 0.03: return 'SUPPORTED: H clearly stronger than any H/j'
    return 'AMBIGUOUS'
c = collections.Counter(cls(o) for o in out)
print(f'{len(out)} arrays flagged HOR (hor >= 1.5 x mono) on the zebra finch panel, v1.8.1 defaults\n')
for k, v in c.most_common(): print(f'  {v:3d}  {k}')
print('\nper-array detail (sorted by H/P):')
print(f"{'array':34s} {'L':>8s} {'mono':>5s} {'HOR':>6s} {'H/P':>6s} {'copies':>7s} {'exc(H)':>7s} {'exc(P)':>7s} {'best H/j':>14s} {'class'}")
for o in sorted(out, key=lambda x: -x['ratio']):
    print(f"{o['id'][:34]:34s} {o['L']:8d} {o['P']:5d} {o['H']:6d} {o['ratio']:6.1f} {o['ncopies']:7.1f} {o['eH']:7.3f} {o['eP']:7.3f}   j={o['best_sub'][0]} {o['best_sub'][1]:6.3f}   {cls(o).split(':')[0]}")
