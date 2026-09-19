#!/usr/bin/env python3
"""Compare ArraySplitter --method autocorr vs --method classic on the same panel.

Two outputs:
- per-array TSV: array_id, autocorr_hor, autocorr_mono, classic_hor, classic_mono,
  hor_agreement (within ±5/5%), mono_agreement
- summary JSON: counts of arrays where the two methods agree / disagree on period.

stdlib only.
"""
import argparse, csv, json, statistics, sys

csv.field_size_limit(sys.maxsize)


def read_summary(path):
    out = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            try:
                hp = int(row.get('hor_period') or 0)
                mp = int(row.get('mono_period') or 0)
                nm = int(row.get('hor_n_monomers') or 0)
            except ValueError:
                continue
            out[row['array_id']] = {'hor': hp, 'mono': mp, 'n_mono': nm,
                                     'autocorr': float(row.get('hor_autocorr') or 0)}
    return out


def period_agree(a, b, tol_bp=5, tol_pct=0.05):
    if a == 0 or b == 0:
        return False
    return abs(a - b) <= max(tol_bp, tol_pct * max(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--autocorr-summary', required=True)
    ap.add_argument('--classic-summary', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    auto = read_summary(args.autocorr_summary)
    classic = read_summary(args.classic_summary)
    common = sorted(set(auto) & set(classic))

    rows = []
    hor_agree_n = 0
    mono_agree_n = 0
    both_agree_n = 0
    auto_has_hor = 0
    classic_has_hor = 0
    for aid in common:
        a = auto[aid]; c = classic[aid]
        ha = period_agree(a['hor'], c['hor'])
        ma = period_agree(a['mono'], c['mono'])
        if ha: hor_agree_n += 1
        if ma: mono_agree_n += 1
        if ha and ma: both_agree_n += 1
        if a['hor'] > a['mono'] * 1.5: auto_has_hor += 1
        if c['hor'] > c['mono'] * 1.5: classic_has_hor += 1
        rows.append({'array_id': aid, 'auto_hor': a['hor'], 'auto_mono': a['mono'],
                     'classic_hor': c['hor'], 'classic_mono': c['mono'],
                     'hor_agree': ha, 'mono_agree': ma,
                     'auto_has_hor': a['hor'] > a['mono'] * 1.5,
                     'classic_has_hor': c['hor'] > c['mono'] * 1.5})

    summary = {
        'n_auto_arrays':    len(auto),
        'n_classic_arrays': len(classic),
        'n_common':         len(common),
        'auto_only':        len(set(auto) - set(classic)),
        'classic_only':     len(set(classic) - set(auto)),
        'hor_period_agree_n':  hor_agree_n,
        'hor_period_agree_pct': round(100 * hor_agree_n / max(len(common), 1), 1),
        'mono_period_agree_n': mono_agree_n,
        'mono_period_agree_pct': round(100 * mono_agree_n / max(len(common), 1), 1),
        'both_agree_n':       both_agree_n,
        'autocorr_emits_hor_n': auto_has_hor,
        'classic_emits_hor_n':  classic_has_hor,
    }

    with open(args.out + '.summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    with open(args.out + '.per_array.tsv', 'w') as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader()
            for r in rows:
                w.writerow(r)

    print(json.dumps(summary, indent=2))
    print(f'\nWrote {args.out}.summary.json + {args.out}.per_array.tsv')


if __name__ == '__main__':
    main()
