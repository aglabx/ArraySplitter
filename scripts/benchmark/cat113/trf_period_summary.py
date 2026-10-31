#!/usr/bin/env python3
"""Per array: which periods TRF finds and how much of the array each covers (trf_on_arrays.sh output).

    python3 -I trf_period_summary.py TRF_ARRAYS.tsv ARRAYS.bed BINS > SUMMARY.tsv

BINS: comma-separated period ranges, e.g. 100-126,440-530. Rows per array: the array length; the union of all TRF
repeats; for each bin the union of the repeats whose period lies in it and their best consensus (highest score);
then the top 5 periods by union coverage. TRF reports overlapping repeats at several periods for one region, so
coverages of different periods add up to more than the array. A label in the TSV but not in the BED, or a row whose
numbers do not parse, stops the run. A BED label with no TRF row is printed as no_trf_rows, never as 0 % (TRF may
have failed on it, or found nothing; trf_on_arrays.sh lists failures in FAILED.tsv).
"""
import sys
from collections import defaultdict


def union(iv):
    tot, cur_s, cur_e = 0, None, None
    for s, e in sorted(iv):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                tot += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s + 1
    return tot


def main():
    tsv, bed, bins_arg = sys.argv[1:4]
    bins = []
    for b in bins_arg.split(","):
        lo, hi = b.split("-")
        bins.append((int(lo), int(hi)))
    length, order = {}, []
    with open(bed) as f:
        for line in f:
            if not line.strip():
                continue
            _, s, e, label = line.rstrip("\n").split("\t")[:4]
            length[label] = int(e) - int(s)
            order.append(label)
    rows = defaultdict(list)
    with open(tsv) as f:
        head = next(f).rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["label"] not in length:
                sys.exit(f"trf_period_summary: label {r['label']} not in {bed}")
            rows[r["label"]].append((int(r["start"]), int(r["end"]), int(r["period"]), float(r["score"]),
                                     float(r["pct_match"]), r["consensus"]))
    print("label\tlength\ttrf_rows\tunion_all_bp\tunion_all_pct\t" +
          "\t".join(f"bin_{lo}_{hi}_bp\tbin_{lo}_{hi}_pct\tbin_{lo}_{hi}_best" for lo, hi in bins) +
          "\ttop_periods_by_union")
    for label in order:
        rr = rows.get(label, [])
        L = length[label]
        if not rr:
            print("\t".join([label, str(L), "no_trf_rows"]))
            continue
        out = [label, str(L), str(len(rr))]
        u = union([(s, e) for s, e, *_ in rr])
        out += [str(u), f"{100 * u / L:.1f}"]
        for lo, hi in bins:
            sel = [x for x in rr if lo <= x[2] <= hi]
            ub = union([(s, e) for s, e, *_ in sel])
            best = max(sel, key=lambda x: x[3]) if sel else None
            out += [str(ub), f"{100 * ub / L:.1f}",
                    f"p{best[2]} {best[4]:.0f}% {best[5]}" if best else "-"]
        byp = defaultdict(list)
        for s, e, p, *_ in rr:
            byp[p].append((s, e))
        top = sorted(((union(v), p) for p, v in byp.items()), reverse=True)[:5]
        out.append("; ".join(f"{p}:{100 * c / L:.1f}%" for c, p in top))
        print("\t".join(out))


if __name__ == "__main__":
    main()
