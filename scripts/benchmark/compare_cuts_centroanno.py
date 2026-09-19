#!/usr/bin/env python3
"""Compare ArraySplitter cuts vs centroAnno cuts on a panel.

Mirror of state/scripts/smoke_asplit_vs_fastan/compare_cuts_zfinch.py but for
centroAnno output instead of FasTAN. centroAnno's per-array per-monomer
decomposition is direct (no GDB / ONEview / M+P line gymnastics) — each
`<array>_decomposedResult.csv` is one row per detected monomer:

    array_id, template, start, end, identity, monomer_length

So the "cuts" for an array are just the `start` column, and the period is
the mode of `monomer_length`.

Outputs three files per panel:
- <out>.summary.json   — top-line metrics (matches the §3.1 table cells)
- <out>.per_array.tsv  — per-array period delta + per-array Jaccard
- stdout               — human-readable summary

Python stdlib only.
"""
import argparse, bisect, csv, collections, glob, json, os, statistics, sys

# ArraySplitter hors.tsv / monomers.tsv carry per-monomer consensus strings that
# can run to multiple megabytes on big HOR arrays — Python's default
# csv.field_size_limit (131072) chokes on them.
csv.field_size_limit(sys.maxsize)


# ---------- ArraySplitter ----------

def parse_asplit(summary_tsv, hors_tsv):
    """Read ArraySplitter summary + hors → {aid: {period, cuts}}.

    Mirrors compare_cuts_zfinch.py logic: period is `hor_period or mono_period`
    from summary.tsv; cuts are reconstructed as cumulative monomer lengths
    inside each array, with the left-flank offset accounted for.
    """
    asplit = {}
    with open(summary_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            aid = row['array_id']
            try:
                hp = int(row.get('hor_period') or 0)
                mp = int(row.get('mono_period') or 0)
            except ValueError:
                continue
            p = hp if hp else mp
            if p <= 0:
                continue
            asplit[aid] = {'period': p, 'cuts': []}

    with open(hors_tsv) as f:
        cur_aid = None
        flank = 0
        cuts = []
        def flush(aid):
            if aid in asplit:
                asplit[aid]['cuts'] = cuts[:]
        for row in csv.DictReader(f, delimiter='\t'):
            aid = row['array_id']
            if aid != cur_aid:
                if cur_aid is not None:
                    flush(cur_aid)
                cur_aid = aid
                flank = 0
                cuts = []
            t = row.get('type', '')
            try:
                L = int(row['length'])
            except (KeyError, ValueError):
                continue
            if t == 'flank' and row.get('source') == 'left_flank':
                flank = L
            elif t == 'monomer':
                cuts.append(flank + L if not cuts else cuts[-1] + L)
        if cur_aid is not None:
            flush(cur_aid)
    return asplit


# ---------- centroAnno ----------

def parse_centroanno(centro_dir):
    """Walk *_decomposedResult.csv → {aid: {period, cuts}}.

    Per-row format (no header): aid, template, start, end, identity, monomer_length.
    `period` is the mode of monomer_length; `cuts` is the sorted `start` column.
    array_id is taken from the filename rather than the first CSV column, because
    centroAnno can write a "region-qualified" array_id inside (rare on this panel
    but happens for arrays split into >1 detected region) — the filename is the
    canonical input fasta record id.
    """
    out = {}
    files = sorted(glob.glob(os.path.join(centro_dir, '*_decomposedResult.csv')))
    for path in files:
        fname = os.path.basename(path)
        aid = fname[:-len('_decomposedResult.csv')]
        lengths = []
        cuts = []
        with open(path) as f:
            for line in f:
                parts = line.rstrip('\n').split(',')
                if len(parts) < 6:
                    continue
                try:
                    start = int(parts[2])
                    L = int(parts[5])
                except ValueError:
                    continue
                cuts.append(start)
                lengths.append(L)
        if not lengths:
            continue
        # mode of monomer length (centroAnno's notion of "monomer")
        period = collections.Counter(lengths).most_common(1)[0][0]
        cuts.sort()
        out[aid] = {'period': period, 'cuts': cuts, 'n_rows': len(lengths)}
    return out


# ---------- comparison ----------

def fold_delta(d, P):
    return ((d + P // 2) % P) - P // 2 if P > 0 else d


def nearest_phase_delta_fast(c, sorted_folded, P):
    """Min |fold_delta(c - o, P)| over sorted_folded[].

    sorted_folded is the asplit cuts collapsed by `x % P`, sorted ASC. The
    nearest folded delta is between (c % P) and the closest entry in
    sorted_folded, accounting for the circular wrap at P.
    """
    cm = c % P
    idx = bisect.bisect_left(sorted_folded, cm)
    candidates = []
    if idx < len(sorted_folded):
        candidates.append(sorted_folded[idx])
    if idx > 0:
        candidates.append(sorted_folded[idx - 1])
    if sorted_folded:
        candidates.append(sorted_folded[0] + P)  # wrap forward
        candidates.append(sorted_folded[-1] - P)  # wrap backward
    best = None
    for o in candidates:
        d = cm - o
        # Already lives in [-P, +P]; collapse to [-P/2, +P/2] once.
        if d > P // 2:
            d -= P
        elif d < -P // 2:
            d += P
        if best is None or abs(d) < abs(best):
            best = d
    return best if best is not None else 0


def compare(asplit, centro, label):
    common = sorted(set(asplit) & set(centro))
    print(f'[{label}]', file=sys.stderr)
    print(f'  ArraySplitter arrays: {len(asplit)}', file=sys.stderr)
    print(f'  centroAnno    arrays: {len(centro)}', file=sys.stderr)
    print(f'  Common arrays:        {len(common)}', file=sys.stderr)

    period_match_strict = 0
    period_match_loose = 0
    period_deltas = []
    n_cuts_deltas = []
    all_phase_deltas = []
    per_array = []
    no_cuts_pair = 0

    for aid in common:
        a = asplit[aid]; c = centro[aid]
        pa, pc = a['period'], c['period']
        pd = pc - pa
        period_deltas.append(pd)
        if abs(pd) <= max(5, 0.05 * pa):
            period_match_strict += 1
        if abs(pd) <= max(20, 0.10 * pa):
            period_match_loose += 1

        na, nc = len(a['cuts']), len(c['cuts'])
        n_cuts_deltas.append(nc - na)

        # Period-aware Jaccard requires both cut sets non-empty and comparable period
        if a['cuts'] and c['cuts'] and abs(pd) <= max(20, 0.10 * pa):
            P = pa
            sorted_folded = sorted(x % P for x in a['cuts'])
            phase = [nearest_phase_delta_fast(c2, sorted_folded, P) for c2 in c['cuts']]
            all_phase_deltas.extend(phase)
            med = statistics.median(phase) if phase else 0
            per_array.append({
                'array_id': aid,
                'asplit_period': pa, 'centro_period': pc,
                'asplit_cuts': na, 'centro_cuts': nc,
                'phase_median_bp': round(med, 1),
                'pct_within_0bp':  round(100 * sum(1 for x in phase if abs(x) <= 0)  / len(phase), 1),
                'pct_within_10bp': round(100 * sum(1 for x in phase if abs(x) <= 10) / len(phase), 1),
                'pct_within_100bp':round(100 * sum(1 for x in phase if abs(x) <= 100)/ len(phase), 1),
            })
        else:
            no_cuts_pair += 1
            per_array.append({
                'array_id': aid,
                'asplit_period': pa, 'centro_period': pc,
                'asplit_cuts': na, 'centro_cuts': nc,
                'phase_median_bp': None,
                'pct_within_0bp': None,
                'pct_within_10bp': None,
                'pct_within_100bp': None,
            })

    n = len(period_deltas) or 1
    out = {
        'panel': label,
        'n_asplit': len(asplit),
        'n_centro': len(centro),
        'n_common': len(common),
        'period_match_strict_pct': round(100 * period_match_strict / n, 1),
        'period_match_loose_pct':  round(100 * period_match_loose / n, 1),
        'period_match_strict_n':   period_match_strict,
        'period_match_loose_n':    period_match_loose,
        'n_cuts_delta_median':     statistics.median(n_cuts_deltas) if n_cuts_deltas else None,
        'n_phase_deltas':          len(all_phase_deltas),
        'arrays_in_phase_pool':    sum(1 for r in per_array if r['phase_median_bp'] is not None),
        'arrays_skipped_no_overlap': no_cuts_pair,
    }
    if all_phase_deltas:
        for tol in (0, 10, 100):
            inside = sum(1 for d in all_phase_deltas if abs(d) <= tol)
            out[f'cut_jaccard_pct_{tol}bp'] = round(100 * inside / len(all_phase_deltas), 1)
    return out, per_array, all_phase_deltas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--panel', required=True, help='label, e.g. zfinch or hg002')
    ap.add_argument('--centro-dir', required=True)
    ap.add_argument('--asplit-summary', required=True)
    ap.add_argument('--asplit-hors', required=True)
    ap.add_argument('--out', required=True, help='output prefix')
    args = ap.parse_args()

    asplit = parse_asplit(args.asplit_summary, args.asplit_hors)
    centro = parse_centroanno(args.centro_dir)
    summary, per_array, _phase_deltas = compare(asplit, centro, args.panel)

    # Persist
    with open(args.out + '.summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    with open(args.out + '.per_array.tsv', 'w') as f:
        if per_array:
            w = csv.DictWriter(f, fieldnames=list(per_array[0].keys()), delimiter='\t')
            w.writeheader()
            for r in per_array:
                w.writerow(r)

    # Stdout block — paper-grade
    print()
    print(f'=== {args.panel} : ArraySplitter vs centroAnno ===')
    print(f'  ArraySplitter arrays with period:  {summary["n_asplit"]}')
    print(f'  centroAnno    arrays with period:  {summary["n_centro"]}')
    print(f'  Common arrays:                     {summary["n_common"]}')
    print(f'  Period agreement (|delta| <= max(5,5%)):  '
          f'{summary["period_match_strict_n"]}/{summary["n_common"]} '
          f'({summary["period_match_strict_pct"]}%)')
    print(f'  Period agreement (|delta| <= max(20,10%)):'
          f'{summary["period_match_loose_n"]}/{summary["n_common"]} '
          f'({summary["period_match_loose_pct"]}%)')
    print(f'  Period-aware cut pool: {summary["n_phase_deltas"]} cuts '
          f'across {summary["arrays_in_phase_pool"]} arrays')
    if 'cut_jaccard_pct_0bp' in summary:
        print(f'  Cut Jaccard +/-0 bp:  {summary["cut_jaccard_pct_0bp"]}%')
        print(f'  Cut Jaccard +/-10 bp: {summary["cut_jaccard_pct_10bp"]}%')
        print(f'  Cut Jaccard +/-100 bp:{summary["cut_jaccard_pct_100bp"]}%')
    print()
    print(f'Wrote: {args.out}.summary.json + {args.out}.per_array.tsv')


if __name__ == '__main__':
    main()
