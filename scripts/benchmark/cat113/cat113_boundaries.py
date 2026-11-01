#!/usr/bin/env python3
"""ArraySplitter monomer boundaries on the cat 113-bp fields against nhmmer monomer calls (round 5, item 3).

    python3 -I cat113_boundaries.py RUN_DIR INPUTS REFERENCE_DIR > boundaries_vs_nhmmer.txt

RUN_DIR: run_cat113.sh workdir (asplit_new/, asplit_v182/); INPUTS: pat.fa, mat.fa; REFERENCE_DIR:
scripts/benchmark/cat113/reference (calls_pat.bed, calls_mat.bed, 0-based BED in assembly coordinates).
Lanes: v1.9.0 and v1.8.2 at D 32000. Scored fields: status ok and a 108-116-bp monomer period.
Per lane and haplotype:
  placement   level-1 units from the exact top-level tiling in .lengths; deeper units and base monomers by
              their sequence inside the parent window; unplaced monomers are named and excluded.
  raw         ArraySplitter boundaries inside the span of nhmmer calls; share within 5 / 10 bp of the nearest nhmmer
              boundary; null = the same boundaries shifted by every offset 1..P-1, averaged.
  offset      per field, the modal signed offset to the nearest nhmmer boundary is removed (one number per field);
              share within 5 / 10 bp after it; null = share of span positions within the same tolerance of a shifted
              nhmmer boundary (a random position).
  matching    nhmmer call matched when an ArraySplitter monomer has both ends within 10 bp of it (raw and after the
              field offset); recall = matched calls / calls in scored fields; precision = matched monomers / monomers
              inside the span.
Tolerances: ASM-cat113-boundary-tolerance. A field id missing from the FASTA, or a call outside every field, stops.
"""
import bisect
import csv
import os
import re
import sys
from collections import Counter, defaultdict

csv.field_size_limit(1 << 30)
RUN, INP, REF = sys.argv[1:4]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TOL = (5, 10)
MATCH = 10


def fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def field_coords(aid):
    m = re.match(r"[^|]+\|([^:]+):(\d+)-(\d+)\|", aid)
    if not m:
        sys.exit(f"cat113_boundaries: cannot parse coordinates of {aid}")
    return m.group(1), int(m.group(2)) - 1, int(m.group(3))  # acc, 0-based start, end


def read_lengths(path):
    out, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = line[1:].split()[0]
        elif cur is not None:
            out[cur] = [int(x) for x in line.split()]
    return out


def place(seq, orient, top_lengths, hors_rows, mono_rows):
    """Base monomers as [s, e) in forward field coordinates.

    Level-1 units are positioned from the exact top-level tiling in .lengths (oriented array): the units of type
    `monomer`, sorted by idx, must form one contiguous run of it, flanks fill the ends. Deeper units and base
    monomers are found by their sequence inside their parent's window, in file order. Returns (intervals, unplaced).
    """
    L = len(seq)
    o = seq if orient == "fwd" else seq.translate(COMP)[::-1]
    if sum(top_lengths) != L:
        sys.exit(f"cat113_boundaries: .lengths sums to {sum(top_lengths)}, field is {L}")
    units = sorted((int(r["idx"]), int(r["length"])) for r in hors_rows if r["level"] == "1" and r["type"] == "monomer")
    ul = [l for _, l in units]
    k = next((k for k in range(len(top_lengths) - len(ul) + 1) if top_lengths[k:k + len(ul)] == ul), None)
    if k is None:
        return [], len(mono_rows)
    pos = sum(top_lengths[:k])
    win = {}
    for idx, l in units:
        win[(1, idx)] = [pos, pos + l, pos]  # start, end, cursor
        pos += l
    for r in hors_rows:  # deeper HOR levels, in file order
        lev = r["level"]
        if lev in ("1", "-", "NA", "") or r["type"] != "monomer":
            continue
        par = win.get((int(lev) - 1, int(r["parent_idx"])))
        if par is None:
            continue
        p = o.find(r["sequence"], par[2], par[1])
        if p >= 0:
            win[(int(lev), int(r["idx"]))] = [p, p + len(r["sequence"]), p]
            par[2] = p + len(r["sequence"])
    out, miss = [], 0
    for r in mono_rows:
        par = win.get((int(r["parent_level"]), int(r["parent_idx"])))
        p = o.find(r["sequence"], par[2], par[1]) if par else -1
        if p < 0:
            miss += 1
            continue
        par[2] = p + len(r["sequence"])
        out.append((p, p + len(r["sequence"])) if orient == "fwd" else (L - p - len(r["sequence"]), L - p))
    return sorted(out), miss


def nearest(sorted_b, x):
    i = bisect.bisect_left(sorted_b, x)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(sorted_b):
            d = x - sorted_b[j]
            if best is None or abs(d) < abs(best):
                best = d
    return best


def run_lane(lane, hap, seqs, calls_by_field):
    summ = {r["array_id"]: r for r in csv.DictReader(open(os.path.join(RUN, lane, f"{hap}_d32k.summary.tsv")), delimiter="\t")}
    monos = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(RUN, lane, f"{hap}_d32k.monomers.tsv")), delimiter="\t"):
        if r["type"] == "base_monomer" and r["parent_idx"] not in ("-", "NA", ""):
            monos[r["array_id"]].append(r)
    hrows = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(RUN, lane, f"{hap}_d32k.hors.tsv")), delimiter="\t"):
        hrows[r["array_id"]].append(r)
    tops = read_lengths(os.path.join(RUN, lane, f"{hap}_d32k.lengths"))
    agg = Counter()
    raw_hits = {t: 0 for t in TOL}; raw_null = {t: 0.0 for t in TOL}
    off_hits = {t: 0 for t in TOL}; off_null = {t: 0.0 for t in TOL}
    n_b = 0; n_fields = 0; unplaced = []
    matched_calls = matched_calls_off = n_calls = 0
    matched_mono = matched_mono_off = n_mono_span = 0
    for aid, s in summ.items():
        if s["status"] != "ok":
            continue
        P = int(float(s["mono_period"]))
        if not 108 <= P <= 116:
            continue
        if aid not in seqs:
            sys.exit(f"cat113_boundaries: {aid} not in FASTA")
        calls = calls_by_field.get(aid, [])
        if not calls:
            agg["scored_fields_without_calls"] += 1
            continue
        iv, miss = place(seqs[aid], s["orientation"], tops[aid], hrows[aid], monos[aid])
        if miss:
            unplaced.append(f"{aid.split('|')[0]}:{miss}")
        agg["placed"] += len(iv); agg["unplaced"] += miss
        n_fields += 1
        cb = sorted({x for c in calls for x in c})
        lo, hi = calls[0][0], max(c[1] for c in calls)
        ab = sorted({x for a in iv for x in a})
        inb = [b for b in ab if lo <= b <= hi]
        if not inb:
            continue
        n_b += len(inb)
        d = [nearest(cb, b) for b in inb]
        for t in TOL:
            raw_hits[t] += sum(1 for x in d if abs(x) <= t)
            tot = 0
            for k in range(1, P):
                tot += sum(1 for b in inb if lo <= b + k <= hi and abs(nearest(cb, b + k)) <= t)
            raw_null[t] += tot / (P - 1)
        clipped = [x for x in d if abs(x) <= P // 2]
        off = Counter(clipped).most_common(1)[0][0] if clipped else 0
        d2 = [nearest(cb, b - off) for b in inb]
        shifted = sorted(c + off for c in cb)
        span = hi - lo + 1
        for t in TOL:
            off_hits[t] += sum(1 for x in d2 if abs(x) <= t)
            covered = 0  # positions of the span within t of a shifted nhmmer boundary
            last = lo - 1
            for c in shifted:
                a0, a1 = max(lo, c - t, last + 1), min(hi, c + t)
                if a1 >= a0:
                    covered += a1 - a0 + 1
                    last = a1
            off_null[t] += len(inb) * covered / span
        starts = sorted(iv)
        st = [a for a, _ in starts]
        for (cs, ce) in calls:
            n_calls += 1
            for shift, which in ((0, "raw"), (off, "off")):
                i = bisect.bisect_left(st, cs + shift - MATCH)
                ok = False
                while i < len(st) and st[i] <= cs + shift + MATCH:
                    if abs(starts[i][1] - (ce + shift)) <= MATCH:
                        ok = True
                        break
                    i += 1
                if ok:
                    if which == "raw":
                        matched_calls += 1
                    else:
                        matched_calls_off += 1
        cstarts = sorted(calls)
        cst = [a for a, _ in cstarts]
        for (a, b) in iv:
            if a < lo or b > hi:
                continue
            n_mono_span += 1
            for shift, which in ((0, "raw"), (off, "off")):
                i = bisect.bisect_left(cst, a - shift - MATCH)
                ok = False
                while i < len(cst) and cst[i] <= a - shift + MATCH:
                    if abs(cstarts[i][1] - (b - shift)) <= MATCH:
                        ok = True
                        break
                    i += 1
                if ok:
                    if which == "raw":
                        matched_mono += 1
                    else:
                        matched_mono_off += 1
    print(f"\n## {lane} {hap} D32000: scored fields {n_fields} (+{agg['scored_fields_without_calls']} without calls), "
          f"monomers placed {agg['placed']}, unplaced {agg['unplaced']} {' '.join(unplaced)}")
    if n_b == 0:
        print("no boundaries inside nhmmer spans")
        return
    for t in TOL:
        print(f"raw    within {t:>2} bp: {100 * raw_hits[t] / n_b:5.1f}%   null (random phase) {100 * raw_null[t] / n_b:5.1f}%")
    for t in TOL:
        print(f"offset within {t:>2} bp: {100 * off_hits[t] / n_b:5.1f}%   null (random position) {100 * off_null[t] / n_b:5.1f}%")
    print(f"boundaries scored {n_b}; nhmmer calls {n_calls}: matched raw {matched_calls} ({100 * matched_calls / n_calls:.1f}%), "
          f"after offset {matched_calls_off} ({100 * matched_calls_off / n_calls:.1f}%); ArraySplitter monomers in span {n_mono_span}: "
          f"matched raw {matched_mono} ({100 * matched_mono / max(1, n_mono_span):.1f}%), after offset {matched_mono_off} "
          f"({100 * matched_mono_off / max(1, n_mono_span):.1f}%)")


for hap in ("pat", "mat"):
    seqs = fasta(os.path.join(INP, f"{hap}.fa"))
    fields = [(field_coords(a), a) for a in seqs]
    calls_by_field = defaultdict(list)
    for line in open(os.path.join(REF, f"calls_{hap}.bed")):
        acc, s, e = line.split("\t")[:3]
        s, e = int(s), int(e)
        hit = [a for (fa, fs, fe), a in fields if fa == acc and fs <= s and e <= fe]
        if not hit:
            continue  # calls outside the kept fields (fields need >= MIN_CALLS)
        fs = field_coords(hit[0])[1]
        calls_by_field[hit[0]].append((s - fs, e - fs))
    for a in calls_by_field:
        calls_by_field[a].sort()
    n_in = sum(len(v) for v in calls_by_field.values())
    print(f"# {hap}: nhmmer calls inside fields {n_in}")
    for lane in ("asplit_new", "asplit_v182"):
        run_lane(lane, hap, seqs, calls_by_field)
