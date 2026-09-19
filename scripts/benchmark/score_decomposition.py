#!/usr/bin/env python3
"""Score ArraySplitter, FasTAN, and centroAnno on synthetic tandem-repeat fixtures
with known ground truth.

Inputs (passed via flags):
  --truth          sim_truth.tsv (one row per fixture, columns described in
                   generate_fixtures.py)
  --asplit-summary asplit_sim.summary.tsv
  --asplit-hors    asplit_sim.hors.tsv
  --centro-dir     dir holding centroAnno per-fixture *_decomposedResult.csv
  --fastan-ano     fastan_sim.1ano (decoded via ONEview, optional)
  --out            output prefix

Per-fixture metrics:
  - period_delta_mono = |reported_mono_period − true_mono_period|
  - period_delta_hor  = |reported_hor_period  − true_hor_period|
  - recursion_depth_delta = reported_max_level − true_max_level
  - cut_recall, cut_precision at ±5 bp (computed only on flat / single-level
    fixtures where the true cut grid is unambiguously known)

Aggregate per (tool × difficulty) cell into a summary JSON.
"""
import argparse, csv, collections, glob, json, os, statistics, subprocess, sys

csv.field_size_limit(sys.maxsize)


# ---------- truth ----------

def read_truth(path):
    truth = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            truth[row['fixture_id']] = {
                'archetype':        row['archetype'],
                'array_length_bp':  int(row['array_length_bp']),
                'mono_period':      int(row['true_mono_period']),
                'hor_period':       int(row['true_hor_period']),
                'sub_period':       int(row['true_sub_period']),
                'n_monomers':       int(row['n_monomers']),
                'divergence_pct':   float(row['divergence_pct']),
                'difficulty':       row['difficulty'],
                'flank_bp':         int(row['flank_bp']),
            }
    return truth


# ---------- ArraySplitter ----------

def read_asplit(summary_path, hors_path):
    out = {}
    with open(summary_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            aid = row['array_id']
            try:
                hp = int(row.get('hor_period') or 0)
                mp = int(row.get('mono_period') or 0)
            except ValueError:
                continue
            # Levels — count distinct periods in `period_classes` if present
            classes = row.get('period_classes', '')
            level_count = 1
            if classes:
                level_count = max(1, len([c for c in classes.split(',') if c.strip()]))
            out[aid] = {
                'mono_period': mp,
                'hor_period':  hp,
                'level_count': level_count,
                'cuts':        [],
            }

    # Reconstruct cut positions from hors.tsv (cumsum of monomer lengths after left flank)
    with open(hors_path) as f:
        cur_aid = None
        flank = 0
        cuts = []
        def flush(aid):
            if aid in out:
                out[aid]['cuts'] = cuts[:]
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
    return out


# ---------- centroAnno ----------

def read_centro(centro_dir):
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
        period = collections.Counter(lengths).most_common(1)[0][0]
        cuts.sort()
        out[aid] = {'mono_period': period, 'hor_period': period, 'level_count': 1, 'cuts': cuts}
    return out


# ---------- FasTAN ----------

def read_fastan(ano_path, oneview_path='/home/akomissarov/asplit_bench/FASTGA/ONEview'):
    """Parse FasTAN .1ano (decoded via ONEview) into {aid: {mono_period, cuts}}.

    FasTAN ANO format used here:
      S <len> <name>      — scaffold record
      M <idx> <beg> <end> — repeat region (idx = 0-based scaffold index)
      L <count> <unit_bp> — detected unit length(s); we take the dominant one
      P <n> <c1 c2 ...>   — cut positions RELATIVE to M's <beg>
    """
    if not os.path.exists(ano_path):
        return {}
    ano = subprocess.check_output([oneview_path, ano_path], text=True)

    out = {}
    scaf_names = []
    cur_aid = None
    cur_beg = None
    cur_unit = None
    raw = collections.defaultdict(lambda: {'units': [], 'cuts': []})

    def _flush():
        if cur_aid is None or cur_unit is None:
            return
        raw[cur_aid]['units'].append(cur_unit)
        # cuts in cur stayed local — we attached them after each P line directly to raw

    for line in ano.splitlines():
        if line.startswith('S '):
            parts = line.split(maxsplit=2)
            if len(parts) >= 3:
                scaf_names.append(parts[2])
        elif line.startswith('M '):
            _flush()
            parts = line.split()
            idx = int(parts[1])
            cur_aid = scaf_names[idx] if idx < len(scaf_names) else f'unknown_{idx}'
            cur_beg = int(parts[2])
            cur_unit = None
        elif line.startswith('L ') and cur_aid is not None:
            parts = line.split()
            cur_unit = int(parts[2])
        elif line.startswith('P ') and cur_aid is not None and cur_beg is not None:
            parts = line.split()
            n = int(parts[1])
            cuts = [cur_beg + int(x) for x in parts[2:2 + n]]
            raw[cur_aid]['cuts'].extend(cuts)
    _flush()

    for aid, e in raw.items():
        if not e['units']:
            continue
        # Dominant unit length (mode)
        unit = collections.Counter(e['units']).most_common(1)[0][0]
        out[aid] = {'mono_period': unit, 'hor_period': unit, 'level_count': 1,
                    'cuts': sorted(e['cuts'])}
    return out


# ---------- score ----------

def score_period(reported, truth_val):
    """Return |delta|, but mark missing reports as None and 'failure-test' truth as N/A."""
    if reported in (0, None):
        return None
    if truth_val in (-1, 0):
        return None
    return abs(reported - truth_val)


def expected_cuts(truth):
    """Given the truth dict, generate the expected cut positions (start-of-each-monomer,
    not counting the trailing flank). Only meaningful for flat / single-level fixtures
    where the cut grid is unambiguous (i.e., n_monomers > 0)."""
    if truth['n_monomers'] <= 0 or truth['mono_period'] <= 0:
        return None
    flank = truth['flank_bp']
    period = truth['mono_period']
    return [flank + i * period for i in range(truth['n_monomers'])]


def cut_recall_precision(true_cuts, reported_cuts, tol=5):
    if not true_cuts or not reported_cuts:
        return None, None
    reported_set = sorted(reported_cuts)
    # nearest-neighbour matching with tolerance
    matched_true = 0
    for tc in true_cuts:
        # binary-search nearest reported cut
        import bisect
        idx = bisect.bisect_left(reported_set, tc)
        for k in (idx - 1, idx):
            if 0 <= k < len(reported_set) and abs(reported_set[k] - tc) <= tol:
                matched_true += 1
                break
    matched_reported = 0
    for rc in reported_cuts:
        idx = sorted(true_cuts).__contains__  # we can use plain set with rounding
        # easier: linear search bounded by tol since true_cuts is small here
        for tc in true_cuts:
            if abs(tc - rc) <= tol:
                matched_reported += 1
                break
    recall = matched_true / len(true_cuts)
    precision = matched_reported / len(reported_cuts)
    return recall, precision


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True)
    ap.add_argument('--asplit-summary', required=True)
    ap.add_argument('--asplit-hors', required=True)
    ap.add_argument('--centro-dir')
    ap.add_argument('--fastan-ano')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    truth = read_truth(args.truth)
    asplit = read_asplit(args.asplit_summary, args.asplit_hors)
    centro = read_centro(args.centro_dir) if args.centro_dir else {}
    fastan = read_fastan(args.fastan_ano) if args.fastan_ano else {}

    per_fix = []
    for fid, t in sorted(truth.items()):
        row = {'fixture_id': fid, 'archetype': t['archetype'], 'difficulty': t['difficulty'],
               'true_mono': t['mono_period'], 'true_hor': t['hor_period'],
               'true_sub': t['sub_period'], 'true_n_monomers': t['n_monomers']}

        true_cuts = expected_cuts(t)

        for tool, store in [('asplit', asplit), ('centro', centro), ('fastan', fastan)]:
            r = store.get(fid)
            if not r:
                row[f'{tool}_mono']   = None
                row[f'{tool}_hor']    = None
                row[f'{tool}_levels'] = 0
                row[f'{tool}_delta_mono'] = None
                row[f'{tool}_delta_hor']  = None
                row[f'{tool}_recall']    = None
                row[f'{tool}_precision'] = None
                continue
            row[f'{tool}_mono']   = r['mono_period']
            row[f'{tool}_hor']    = r['hor_period']
            row[f'{tool}_levels'] = r.get('level_count', 1)
            row[f'{tool}_delta_mono'] = score_period(r['mono_period'], t['mono_period'])
            row[f'{tool}_delta_hor']  = score_period(r['hor_period'],  t['hor_period'])
            recall, precision = cut_recall_precision(true_cuts, r.get('cuts') or [])
            row[f'{tool}_recall']    = recall
            row[f'{tool}_precision'] = precision
        per_fix.append(row)

    # Aggregate per difficulty
    by_diff = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in per_fix:
        d = r['difficulty']
        for tool in ('asplit', 'centro', 'fastan'):
            for metric in ('delta_mono', 'delta_hor', 'recall', 'precision'):
                v = r[f'{tool}_{metric}']
                if v is not None:
                    by_diff[(d, tool)][metric].append(v)

    summary = {'n_fixtures': len(per_fix), 'difficulties': {}, 'per_fixture': per_fix}
    for (d, tool), vals in by_diff.items():
        block = summary['difficulties'].setdefault(d, {}).setdefault(tool, {})
        for metric, lst in vals.items():
            block[f'{metric}_n']      = len(lst)
            block[f'{metric}_mean']   = round(statistics.mean(lst), 3) if lst else None
            block[f'{metric}_median'] = round(statistics.median(lst), 3) if lst else None
            block[f'{metric}_max']    = round(max(lst), 3) if lst else None

    with open(args.out + '.summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    keys = ['fixture_id', 'archetype', 'difficulty', 'true_mono', 'true_hor', 'true_sub',
            'true_n_monomers']
    for tool in ('asplit', 'centro', 'fastan'):
        keys += [f'{tool}_mono', f'{tool}_hor', f'{tool}_levels',
                 f'{tool}_delta_mono', f'{tool}_delta_hor',
                 f'{tool}_recall', f'{tool}_precision']
    with open(args.out + '.per_fixture.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader()
        for r in per_fix:
            w.writerow(r)

    # Human-readable stdout
    print()
    print('=== Per-fixture (compact) ===')
    print('fixture                       diff         true_m  true_h   as_m   as_h   ce_m   ft_m')
    for r in per_fix:
        am = r['asplit_mono'] or '-'
        ah = r['asplit_hor']  or '-'
        cm = r['centro_mono'] or '-'
        fm = r['fastan_mono'] or '-'
        print(f"{r['fixture_id']:30s}{r['difficulty']:13s}"
              f"{r['true_mono']:>7}{r['true_hor']:>8}"
              f"{str(am):>7}{str(ah):>7}{str(cm):>7}{str(fm):>7}")
    print()
    print('=== Aggregate by difficulty (median |delta_mono| in bp) ===')
    for d, by in sorted(summary['difficulties'].items()):
        as_med = by.get('asplit', {}).get('delta_mono_median', '-')
        ce_med = by.get('centro', {}).get('delta_mono_median', '-')
        ft_med = by.get('fastan', {}).get('delta_mono_median', '-')
        print(f'  {d:15s}  asplit={as_med}  centro={ce_med}  fastan={ft_med}')

    print()
    print(f'Wrote {args.out}.summary.json + {args.out}.per_fixture.tsv')


if __name__ == '__main__':
    main()
