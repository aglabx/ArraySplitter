#!/usr/bin/env python3
"""Score ArraySplitter on the T5 factorial fixtures.
Inputs: factorial_truth.tsv, ArraySplitter summary.tsv + monomers.tsv (same multi-FASTA run).
Outputs: per-class + per-axis metrics (period recovery, HOR-detection, boundary P/R/F1,
FP-rate on negatives, nesting-depth signal)."""
import sys, bisect, statistics
from collections import defaultdict, Counter

def load_truth(path):
    t = {}
    with open(path) as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            t[c[0]] = dict(klass=c[1], mono=int(c[2]), hor=int(c[3]), depth=int(c[4]),
                           periodic=int(c[5]), n_mono=int(c[6]), flank=int(c[7]),
                           bounds=[int(x) for x in c[8].split(",") if x])
    return t

def load_summary(path):
    s = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            fid = c[col["array_id"]]
            def g(n, d=0):
                try: return float(c[col[n]])
                except: return d
            s[fid] = dict(alen=g("array_length"), hor=int(g("hor_period")),
                          mono=int(g("mono_period")), pc=c[col["period_classes"]] if "period_classes" in col else "")
    return s

def load_monomer_cuts(path):
    cuts = defaultdict(list)  # fid -> list of monomer lengths (in order)
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(hdr)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= col["type"]: continue
            if c[col["type"]] == "base_monomer":
                cuts[c[col["array_id"]]].append(int(c[col["length"]]))
    out = {}
    for fid, lens in cuts.items():
        b, pos = [], 0
        for L in lens: pos += L; b.append(pos)
        out[fid] = b[:-1]  # internal cuts
    return out

def _near(b, arr):
    i = bisect.bisect_left(arr, b); best = 10**9
    if i < len(arr): best = min(best, abs(arr[i]-b))
    if i > 0: best = min(best, abs(arr[i-1]-b))
    return best

def boundary_prf(reported, true, tol):
    if not true and not reported: return 1.0, 1.0, 1.0
    if not true: return 0.0, float('nan'), 0.0
    if not reported: return float('nan'), 0.0, 0.0
    rs, ts = sorted(reported), sorted(true)
    prec = sum(_near(b, ts) <= tol for b in rs) / len(rs)
    rec = sum(_near(b, rs) <= tol for b in ts) / len(ts)
    f1 = 2*prec*rec/(prec+rec) if (prec+rec) > 0 else 0.0
    return prec, rec, f1

def boundary_prf_periodaware(reported, true, period, tol):
    """Period-aware boundary F1: the monomer 'first base' is an arbitrary rotation, so
    a constant phase offset must be folded out (paper §3.4). Search the global offset
    (from pairwise cut differences mod period) that maximises F1, then report it."""
    if not true and not reported: return 1.0, 1.0, 1.0
    if not true: return 0.0, float('nan'), 0.0
    if not reported: return float('nan'), 0.0, 0.0
    ts = sorted(true)
    cand = {0}
    for r in reported[:40]:
        for t in ts[:40]:
            cand.add((r - t) % period)
    best = (0.0, 0.0, -1.0)
    for off in cand:
        shifted = [r - off for r in reported]
        p, rc, f1 = boundary_prf(shifted, ts, tol)
        if f1 > best[2]:
            best = (p, rc, f1)
    return best

def main():
    truth = load_truth(sys.argv[1])
    summ = load_summary(sys.argv[2])
    cuts = load_monomer_cuts(sys.argv[3])

    per_class = defaultdict(lambda: dict(n=0, mono_ok=0, hor_ok=0, hordet_tp=0, hordet_fp=0,
                                         hordet_tn=0, hordet_fn=0, f1=[], fp_neg=0))

    for fid, tr in truth.items():
        s = summ.get(fid)
        if s is None: continue
        pc = per_class[tr["klass"]]; pc["n"] += 1

        if tr["periodic"]:
            # mono period recovery
            tol_m = max(2, 0.05*tr["mono"])
            mono_ok = abs(s["mono"] - tr["mono"]) <= tol_m
            pc["mono_ok"] += mono_ok
            # HOR period recovery (only HOR fixtures)
            if tr["hor"] > 0:
                tol_h = max(5, 0.05*tr["hor"])
                pc["hor_ok"] += abs(s["hor"] - tr["hor"]) <= tol_h
            # HOR-detection confusion
            pred_hor = s["hor"] >= 1.5*max(s["mono"], 1)
            true_hor = tr["hor"] > 0
            if pred_hor and true_hor: pc["hordet_tp"] += 1
            elif pred_hor and not true_hor: pc["hordet_fp"] += 1
            elif not pred_hor and not true_hor: pc["hordet_tn"] += 1
            else: pc["hordet_fn"] += 1
            # boundary F1 (restrict to array region if flanked)
            rep = cuts.get(fid, [])
            fl = tr["flank"]; alen = int(s["alen"])
            if fl > 0:
                rep = [x for x in rep if fl <= x <= alen-fl]
            tol_b = max(5, 0.1*tr["mono"])
            p, r, f1 = boundary_prf_periodaware(rep, tr["bounds"], max(tr["mono"], 1), tol_b)
            if f1 == f1: pc["f1"].append(f1)
        else:
            # negative: FP = claimed a real short period (not the array-length degenerate)
            alen = int(s["alen"]) or 1
            fp = 0 < s["mono"] < 0.5*alen
            pc["fp_neg"] += fp

    # ---- report ----
    print("=== T5 factorial benchmark — ArraySplitter accuracy ===\n")
    print(f"{'class':<24} {'n':>4} {'mono%':>6} {'HOR%':>6} {'HORdet(P/R)':>12} {'bndF1':>6}")
    tot_n, tot_mono, tot_f1 = 0, 0, []
    for k in sorted(per_class):
        pc = per_class[k]
        n = pc["n"]
        # mono recovery denom = periodic fixtures in class
        periodic_n = sum(1 for fid, tr in truth.items() if tr["klass"]==k and tr["periodic"])
        hor_n = sum(1 for fid, tr in truth.items() if tr["klass"]==k and tr["hor"]>0)
        mono_pct = 100*pc["mono_ok"]/periodic_n if periodic_n else float('nan')
        hor_pct = 100*pc["hor_ok"]/hor_n if hor_n else float('nan')
        tp,fp,tn,fn = pc["hordet_tp"],pc["hordet_fp"],pc["hordet_tn"],pc["hordet_fn"]
        hp = tp/(tp+fp) if (tp+fp) else float('nan')
        hr = tp/(tp+fn) if (tp+fn) else float('nan')
        f1 = statistics.mean(pc["f1"]) if pc["f1"] else float('nan')
        extra = f"  FP-rate={100*pc['fp_neg']/n:.0f}%" if k.startswith("negative") else ""
        print(f"{k:<24} {n:>4} {mono_pct:>5.0f}% {hor_pct:>5.0f}% {hp:>5.2f}/{hr:<5.2f} {f1:>5.2f}{extra}")
        if periodic_n:
            tot_mono += pc["mono_ok"]; tot_n += periodic_n
        tot_f1 += pc["f1"]

    print()
    print(f"OVERALL (periodic fixtures): mono-period recovery "
          f"{100*tot_mono/tot_n:.1f}% (n={tot_n}); "
          f"mean boundary F1 {statistics.mean(tot_f1):.3f} (n={len(tot_f1)})")
    # FP on negatives overall
    neg = [fid for fid,tr in truth.items() if not tr['periodic']]
    fp = sum(1 for fid in neg if fid in summ and 0 < summ[fid]['mono'] < 0.5*(int(summ[fid]['alen']) or 1))
    print(f"FP-rate on {len(neg)} negative controls: {100*fp/len(neg):.0f}%")

if __name__ == "__main__":
    main()
