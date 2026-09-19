#!/usr/bin/env python3
"""T5: factorial synthetic benchmark for ArraySplitter (replaces the 12 hand-made fixtures).

Constructs synthetic tandem / HOR / nested / negative-control arrays with EXACT ground truth
(monomer period, HOR period, per-monomer boundary positions, nesting depth, periodicity flag),
swept over the axes the reviewer asked for: monomer length, HOR cardinality, sequence divergence,
indel rate, copy number, GC, flank length, nesting depth, ambiguous bases, and negatives.

All fixtures go into ONE multi-FASTA (run ArraySplitter once); truth in a companion TSV.
Deterministic: one fixed master seed; every fixture derives its own seed from its id.
"""
import argparse, random, hashlib

BASES = "ACGT"

def seed_for(fid):
    return int(hashlib.md5(fid.encode()).hexdigest()[:8], 16)

def rand_seq(n, gc, rng):
    # gc = P(G or C); split evenly G/C and A/T
    out = []
    for _ in range(n):
        if rng.random() < gc:
            out.append(rng.choice("GC"))
        else:
            out.append(rng.choice("AT"))
    return "".join(out)

def mutate(seq, sub_rate, rng):
    if sub_rate <= 0:
        return seq
    s = list(seq)
    for i in range(len(s)):
        if rng.random() < sub_rate:
            s[i] = rng.choice([b for b in BASES if b != s[i]])
    return "".join(s)

def add_indels(seq, indel_rate, rng):
    if indel_rate <= 0:
        return seq
    out = []
    for ch in seq:
        r = rng.random()
        if r < indel_rate / 2:
            continue  # deletion
        elif r < indel_rate:
            out.append(rng.choice(BASES))  # insertion before
            out.append(ch)
        else:
            out.append(ch)
    return "".join(out)

def build_array(P, C, K, div, indel, gc, rng, variant_div=0.18):
    """Return (array_seq, boundaries) where boundaries are the 0-based cut positions
    BETWEEN monomers within the array (not incl 0 or the array end). C=1 -> flat tandem;
    C>=2 -> HOR of C fixed distinct variants repeated K times, each copy noised by div."""
    base = rand_seq(P, gc, rng)
    variants = [base] if C == 1 else [mutate(base, variant_div, rng) for _ in range(C)]
    parts, bounds, pos = [], [], 0
    for _ in range(K):
        for i in range(C):
            mono = add_indels(mutate(variants[i], div, rng), indel, rng)
            parts.append(mono)
            pos += len(mono)
            bounds.append(pos)
    seq = "".join(parts)
    return seq, bounds[:-1]  # drop final (array end)

def build_nested(P, C1, C2, K, div, gc, rng):
    """Two-level: sub-HOR = C1 variants; top HOR = C2 distinct sub-HORs; repeated K times."""
    base = rand_seq(P, gc, rng)
    # C2 sub-HORs, each a fixed set of C1 variants
    subhors = []
    for _ in range(C2):
        subhors.append([mutate(base, 0.18, rng) for _ in range(C1)])
    parts, bounds, pos = [], [], 0
    for _ in range(K):
        for sh in subhors:
            for i in range(C1):
                mono = mutate(sh[i], div, rng)
                parts.append(mono); pos += len(mono); bounds.append(pos)
    return "".join(parts), bounds[:-1]

def with_flank(seq, flank, gc, rng):
    if flank <= 0:
        return seq, 0
    lf = rand_seq(flank, gc, rng); rf = rand_seq(flank, gc, rng)
    return lf + seq + rf, flank

def emit(fixtures, fid, klass, P, hor, depth, periodic, seq, bounds, flank):
    fixtures.append(dict(fid=fid, klass=klass, P=P, hor=hor, depth=depth,
                         periodic=int(periodic), n_mono=len(bounds) + 1,
                         flank=flank, seq=seq, bounds=bounds))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-fasta', required=True)
    ap.add_argument('--out-truth', required=True)
    a = ap.parse_args()
    fixtures = []

    # Baseline: P=171, C=12, K=20, div=0.02, indel=0, gc=0.4, flank=0
    def add(fid, klass, P=171, C=12, K=20, div=0.02, indel=0.0, gc=0.4, flank=0,
            nested=None):
        rng = random.Random(seed_for(fid))
        if klass == "negative_random":
            seq = rand_seq(P * C * K, gc, rng)  # P*C*K just a length knob
            emit(fixtures, fid, klass, 0, 0, 0, False, seq, [], 0)
            return
        if klass == "negative_lowcomplexity":
            unit = rng.choice(["A", "AT", "AAAT", "N"])
            seq = (unit * (P * C * K // len(unit) + 1))[:P * C * K]
            emit(fixtures, fid, klass, 0, 0, 0, False, seq, [], 0)
            return
        if nested:
            C1, C2 = nested
            seq, b = build_nested(P, C1, C2, K, div, gc, rng)
            seq, off = with_flank(seq, flank, gc, rng)
            b = [x + off for x in b]
            emit(fixtures, fid, "nested", P, C1 * P, 2, True, seq, b, off)
            return
        seq, b = build_array(P, C, K, div, indel, gc, rng)
        if klass == "ambiguous":
            # inject two N-runs of 40 bp inside the array (should not break period)
            s = list(seq)
            for st in (len(s) // 3, 2 * len(s) // 3):
                for j in range(40):
                    if st + j < len(s):
                        s[st + j] = "N"
            seq = "".join(s)
        seq, off = with_flank(seq, flank, gc, rng)
        b = [x + off for x in b]
        hor = C * P if C >= 2 else 0
        emit(fixtures, fid, klass, P, hor, 1, True, seq, b, off)

    # --- one-axis sweeps around the baseline (3 replicate seeds each) ---
    for rep in range(3):
        for Pv in (6, 42, 100, 171, 340, 500):
            add(f"P{Pv}_r{rep}", "sweep_monolen", P=Pv)
        for Cv in (1, 2, 4, 6, 12, 20):
            add(f"C{Cv}_r{rep}", "sweep_cardinality", C=Cv)
        for dv in (0.0, 0.02, 0.05, 0.10, 0.15):
            add(f"div{int(dv*100)}_r{rep}", "sweep_divergence", div=dv)
        for iv in (0.0, 0.005, 0.01, 0.02):
            add(f"indel{int(iv*1000)}_r{rep}", "sweep_indel", indel=iv)
        for Kv in (3, 5, 10, 20, 60):
            add(f"K{Kv}_r{rep}", "sweep_copies", K=Kv)
        for gv in (0.3, 0.4, 0.5, 0.6):
            add(f"gc{int(gv*100)}_r{rep}", "sweep_gc", gc=gv)
        for fv in (0, 100, 500):
            add(f"flank{fv}_r{rep}", "sweep_flank", flank=fv)

    # --- special classes ---
    for rep in range(3):
        add(f"micro6_r{rep}", "microsatellite", P=6, C=1, K=200)
        add(f"micro2_r{rep}", "microsatellite", P=2, C=1, K=400)
        add(f"ambiguous_r{rep}", "ambiguous")
        add(f"nested2_r{rep}", "nested", P=100, K=10, nested=(4, 3))
        add(f"nested3_r{rep}", "nested", P=80, K=8, nested=(3, 4))  # scored as depth>=2
    for i in range(10):
        add(f"rand_{i}", "negative_random")
    for i in range(5):
        add(f"lowcx_{i}", "negative_lowcomplexity")

    # write
    with open(a.out_fasta, "w") as fa, open(a.out_truth, "w") as tr:
        tr.write("fid\tklass\ttrue_mono\ttrue_hor\tdepth\tperiodic\tn_mono\tflank\tbounds\n")
        for f in fixtures:
            fa.write(f">{f['fid']}\n{f['seq']}\n")
            tr.write("\t".join(str(x) for x in [
                f['fid'], f['klass'], f['P'], f['hor'], f['depth'], f['periodic'],
                f['n_mono'], f['flank'], ",".join(map(str, f['bounds']))]) + "\n")
    print(f"wrote {len(fixtures)} fixtures to {a.out_fasta}")
    from collections import Counter
    for k, n in sorted(Counter(f['klass'] for f in fixtures).items()):
        print(f"  {k}: {n}")

if __name__ == "__main__":
    main()
