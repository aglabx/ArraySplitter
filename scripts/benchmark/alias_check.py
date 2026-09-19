"""Aliasing of the evenly spaced 10,000-position subsampling (arrays > 50 kb) with the repeat period.
When step = L // 10000 is a multiple of the period, every sampled position has the same phase and the
autocorrelation scan reports a wrong period.   Usage: python3 alias_check.py <workdir> <repo>"""
import os, sys, random, subprocess, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from simlib import flat_array, mutate, write_fasta
W, REPO = sys.argv[1], sys.argv[2]; os.makedirs(W + '/alias', exist_ok=True)
rng = random.Random(99); recs = []
for name, P, n in [('alias_P171_step171', 171, 10000), ('control_P171_step173', 171, 10130),
                   ('alias_P191_step191', 191, 10000), ('alias_P340_step340', 340, 10000)]:
    seq, _ = flat_array(P, n, rng, sub=0.05, indel=0.0, flank=100); recs.append((f'{name}_L{len(seq)}', seq))
rng = random.Random(5)
for L in (48000, 60000, 66000, 69996, 75000):
    seq = ''.join(mutate('TTAGGG', rng, sub=0.02) for _ in range(L // 6))[:L]
    step = len(seq) // 10000 if len(seq) > 50000 else 0   # 0 = no subsampling (<= 50 kb)
    recs.append((f'telo_TTAGGG_L{len(seq)}_step{step}', seq))
write_fasta(W + '/alias/alias.fa', recs)
bin_path = REPO + '/src/rust/arraysplitter/target/release/arraysplitter'
if not os.path.isfile(bin_path):
    import shutil
    bin_path = shutil.which('arraysplitter') or bin_path
subprocess.run([bin_path, '-i', W + '/alias/alias.fa', '-o', W + '/alias/out', '-t', '2'],
               check=True, stderr=subprocess.DEVNULL)
csv.field_size_limit(10**9)
import numpy as np
seqs = dict(recs)
def true_ac(seq, d, rev):
    if rev: seq = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    a = np.frombuffer(seq.encode(), dtype=np.uint8); return float((a[:-d] == a[d:]).mean())
print(f"{'array':34s} {'hor_period':>10s} {'reported hor_autocorr':>22s} {'TRUE full autocorr at that lag':>31s} {'mono_period':>11s}")
for r in csv.DictReader(open(W + '/alias/out.summary.tsv'), delimiter='\t'):
    d = int(r['hor_period']); t = true_ac(seqs[r['array_id']], d, r['orientation'] == 'rev')
    print(f"{r['array_id']:34s} {d:10d} {float(r['hor_autocorr']):22.4f} {t:31.4f} {r['mono_period']:>11s}")
