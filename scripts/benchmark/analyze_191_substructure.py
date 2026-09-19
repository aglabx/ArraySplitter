#!/usr/bin/env python3
"""191 bp sub-period forensic: is the short (6-14 bp) signal a real sub-HOR of the
191 bp unit, or over-decomposition of a genuine ~191 bp alpha-sat monomer?

Tests, per array on the 6 chr16/34/35 arrays that reported a short mono_period:
  (A) composition-corrected autocorrelation excess profile d=2..260 -> where are the
      REAL peaks? (expect strong at 191; if 6/10 are weak/broad -> spurious)
  (B) cut monomers at the best period near 191, MSA a sample WITHIN the array
      -> intra-array conservation
  (C) MSA a sample ACROSS all arrays -> cross-chromosome conservation. A real 10 bp
      microsatellite CANNOT produce a specific 191 bp sequence conserved across
      chr16/chr34/chr35; conservation proves a genuine alpha-sat monomer.
  (D) cut the same monomers at period=10 and measure mutual identity -> if the "10 bp
      sub-unit" is real it is conserved; if it is a substring artifact it is not.

Reproducible; writes to results/191_substructure/.
"""
import sys, os, subprocess, argparse, random
import numpy as np

MAP = {65:0,67:1,71:2,84:3}  # A C G T ; others -> 4 (ambiguous)

def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    keep = False
    with open(path) as fh:
        for line in fh:
            if line[0] == '>':
                if keep and name: seqs[name] = ''.join(buf)
                name = line[1:].split()[0]
                keep = any(name.startswith(w) for w in wanted)
                buf = []
            elif keep:
                buf.append(line.strip())
        if keep and name: seqs[name] = ''.join(buf)
    return seqs

def encode(s):
    a = np.frombuffer(s.encode().upper(), dtype=np.uint8)
    out = np.full(a.shape, 4, dtype=np.int8)
    for b,v in MAP.items(): out[a==b] = v
    return out

def autocorr_excess(code, dmax=260):
    valid = code < 4
    # baseline: sum f_i^2 over composition
    counts = np.array([(code==i).sum() for i in range(4)], dtype=float)
    f = counts / counts.sum()
    base = (f**2).sum()
    prof = []
    for d in range(2, dmax+1):
        a, b = code[:-d], code[d:]
        m = valid[:-d] & valid[d:]
        nd = m.sum()
        if nd == 0: prof.append((d, 0.0)); continue
        match = ((a==b) & m).sum()
        prof.append((d, match/nd - base))
    return prof

def top_peaks(prof, k=8):
    return sorted(prof, key=lambda x: -x[1])[:k]

def cut_monomers(seq, period, n, start=0):
    mons = []
    i = start
    while i+period <= len(seq) and len(mons) < n:
        mons.append(seq[i:i+period]); i += period
    return mons

def mafft_identity(seqs, tag, outdir):
    if len(seqs) < 3: return None, None
    fa = os.path.join(outdir, f"{tag}.fa")
    al = os.path.join(outdir, f"{tag}.aln")
    with open(fa,'w') as f:
        for i,s in enumerate(seqs): f.write(f">{tag}_{i}\n{s}\n")
    with open(al,'w') as f:
        subprocess.run(["mafft","--auto","--quiet","--thread","4",fa], stdout=f, check=True)
    aln, nm, bufd = {}, None, []
    with open(al) as fh:
        for line in fh:
            if line[0]=='>':
                if nm: aln[nm]=''.join(bufd)
                nm=line[1:].strip(); bufd=[]
            else: bufd.append(line.strip())
        if nm: aln[nm]=''.join(bufd)
    rows = list(aln.values())
    if len(rows)<2: return None, None
    L = len(rows[0]); import itertools
    ids=[]
    for x,y in itertools.combinations(rows, 2):
        pos = [(1 if a==b else 0) for a,b in zip(x,y) if a!='-' and b!='-']
        if pos: ids.append(sum(pos)/len(pos))
    # consensus
    cons=[]
    for c in range(L):
        col=[r[c] for r in rows if r[c]!='-']
        if col: cons.append(max(set(col), key=col.count))
    return (np.mean(ids) if ids else None), ''.join(cons)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--out', default='results/191_substructure')
    ap.add_argument('--arrays', nargs='+', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    random.seed(42)
    seqs = read_fasta(a.fasta, a.arrays)
    print(f"loaded {len(seqs)} arrays\n")
    across191 = []
    for aid in a.arrays:
        seq = seqs.get(aid) or next((v for k,v in seqs.items() if k.startswith(aid)), None)
        if not seq: print(f"!! {aid} not found"); continue
        code = encode(seq)
        prof = autocorr_excess(code)
        peaks = top_peaks(prof)
        # best period in the 180-200 (alpha) band and in the 4-16 (short) band
        alpha = max((p for p in prof if 180<=p[0]<=200), key=lambda x:x[1])
        short = max((p for p in prof if 4<=p[0]<=16),  key=lambda x:x[1])
        print(f"=== {aid}  len={len(seq)} ===")
        print("  top-8 autocorr-excess peaks (period:excess):",
              ", ".join(f"{d}:{e:.3f}" for d,e in peaks))
        print(f"  alpha band best = {alpha[0]}bp (excess {alpha[1]:.3f}) | "
              f"short band best = {short[0]}bp (excess {short[1]:.3f}) | "
              f"ratio alpha/short = {alpha[1]/short[1]:.2f}" if short[1]>0 else "")
        # cut ~191 monomers from a clean interior window, sample 30
        p191 = alpha[0]
        mons = cut_monomers(seq, p191, 30, start=len(seq)//4)
        idn, cons = mafft_identity(mons, f"{aid}_m191", a.out)
        print(f"  intra-array 191bp-monomer mean pairwise identity = "
              f"{idn:.3f}" if idn else "  (too few monomers)")
        across191.append((aid, mons[0] if mons else None))
        # cut at 10bp, sample 30, mutual identity of the '10bp sub-unit'
        m10 = cut_monomers(seq, 10, 30, start=len(seq)//4)
        id10,_ = mafft_identity(m10, f"{aid}_m10", a.out)
        print(f"  '10bp sub-unit' mutual identity = {id10:.3f} "
              f"(1.0 would mean a real conserved 10bp repeat)\n" if id10 else "\n")
    # cross-array: one 191 monomer from each array, MSA together
    xseqs = [m for _,m in across191 if m]
    idx, consx = mafft_identity(xseqs, "cross_array_191", a.out)
    print("=== CROSS-ARRAY (one 191bp monomer per array) ===")
    print(f"  mean pairwise identity across chr16/34/35 = {idx:.3f}" if idx else "  n/a")
    print(f"  cross-array 191bp consensus:\n  {consx}" if consx else "")

if __name__ == '__main__':
    main()
