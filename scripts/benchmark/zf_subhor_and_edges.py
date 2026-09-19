"""Real-panel checks: (1) byte-identical duplicate records, (2) what the sub_hor rows are,
(3) partial monomers at the edges of top-level HOR units.  Usage: python3 zf_subhor_and_edges.py <panel.fa> <run_prefix>"""
import os, sys, csv, collections, hashlib, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from asparse import load
csv.field_size_limit(10**9)
PANEL, PREFIX = sys.argv[1], sys.argv[2]
rows = {r['array_id']: r for r in csv.DictReader(open(PREFIX + '.summary.tsv'), delimiter='\t')}
fa = {}; n = None
for l in open(PANEL):
    l = l.strip()
    if l.startswith('>'): n = l[1:].split()[0]; fa[n] = []
    else: fa[n].append(l.upper())
h = collections.defaultdict(list)
for k, v in fa.items(): h[hashlib.md5(''.join(v).encode()).hexdigest()].append(k)
dups = [v for v in h.values() if len(v) > 1]
hap = collections.Counter('mat' if '_mat_' in k else ('pat' if '_pat_' in k else '?') for k in fa)
print(f"panel: {len(fa)} records; haplotype labels {dict(hap)}; byte-identical sequences shared by >1 record: {len(dups)} groups / {sum(len(v) for v in dups)} records")
sub = collections.defaultdict(list); lvl = collections.Counter()
for r in csv.DictReader(open(PREFIX + '.hors.tsv'), delimiter='\t'):
    if r['type'] == 'sub_hor':
        sub[r['array_id']].append((int(r['length']), int(r['period']), int(r['n_expected']), int(r['level']))); lvl[int(r['level'])] += 1
allr = [x for v in sub.values() for x in v]; tot = len(allr)
print(f"\nsub_hor rows: {tot} in {len(sub)} arrays; by level {dict(sorted(lvl.items()))}")
cmpn = [l - p for l, p, c, lv in allr]
print("compared positions behind each sub_hor call (length - period): median=%d; <10: %.1f%%; <20: %.1f%%; <30: %.1f%%; >=100: %.1f%%" % (statistics.median(cmpn), *[100 * sum(1 for c in cmpn if c < t) / tot for t in (10, 20, 30)], 100 * sum(1 for c in cmpn if c >= 100) / tot))
print("sub_hor unit length: median=%d bp; <=30 bp: %.1f%%; <=60 bp: %.1f%%; >=500 bp: %.1f%%" % (statistics.median([x[0] for x in allr]), *[100 * sum(1 for x in allr if x[0] <= t) / tot for t in (30, 60)], 100 * sum(1 for x in allr if x[0] >= 500) / tot))
print("children per sub_hor: ==2: %.1f%%; <=3: %.1f%%" % (100 * sum(1 for x in allr if x[2] == 2) / tot, 100 * sum(1 for x in allr if x[2] <= 3) / tot))
t191 = [k for k in sub if 186 <= int(rows[k]['hor_period']) <= 196]
print(f"arrays whose TOP-LEVEL period is 186-196 bp (the 191-bp satellite): {len(t191)}; they hold {sum(len(sub[k]) for k in t191)} of {tot} sub_hor rows ({100*sum(len(sub[k]) for k in t191)/tot:.1f}%)")
print("top arrays by sub_hor rows:")
for k, v in sorted(sub.items(), key=lambda kv: -len(kv[1]))[:8]: print(f"   {k:34s} rows={len(v):5d} mono_period={rows[k]['mono_period']:>4s} hor_period={rows[k]['hor_period']:>5s}")
res = load(PREFIX); tu = pe = lv_ = fr = na = 0
for aid, r in res.items():
    P = int(r['summary']['mono_period']); H = int(r['summary']['hor_period'])
    if P < 98 or H < 1.5 * P: continue
    na += 1
    for a, b, t in r['top']:
        if t != 'monomer': continue
        inner = [(s, e) for s, e in r['base'] if a <= s and e <= b]
        if len(inner) < 2: continue
        tu += 1; lv_ += len(inner); fr += sum(1 for s, e in inner if (e - s) < 0.7 * P)
        if (inner[0][1] - inner[0][0]) < 0.7 * P or (inner[-1][1] - inner[-1][0]) < 0.7 * P: pe += 1
print(f"\nHOR arrays with monomer >= 98 bp: {na}; top-level units split further: {tu}")
print(f"  units whose first or last leaf is a partial monomer (<70% of the monomer period): {pe} ({100*pe/max(1,tu):.1f}%)")
print(f"  leaves shorter than 70% of the monomer period: {fr} of {lv_} ({100*fr/max(1,lv_):.1f}%)")
