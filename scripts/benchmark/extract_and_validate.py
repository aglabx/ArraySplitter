#!/usr/bin/env python3
"""Extract 191 bp HOR consensus sequences from ArraySplitter zfinch output,
MAFFT-align them, build a per-chromosome confidence table, and run a
shuffled-baseline control to test whether the recursive 10-13 bp sub-unit
signal is reproducible vs. a composition-matched random sequence.

Inputs:
  --asplit-summary zf_5000.summary.tsv  (full path on server)
  --asplit-hors    zf_5000.hors.tsv     (full path on server)
  --asplit-bin     /path/to/arraysplitter (with --max-period support)
  --chr-pattern    regex (default chr(16|34|35)_(mat|pat))
  --hor-period-min int (default 180)
  --hor-period-max int (default 210)
  --shuffles       int (default 10)
  --out-prefix     str

Outputs:
  <out>.consensuses.fasta    — extracted 191 bp consensus seqs
  <out>.per_array.tsv        — per-array confidence row
  <out>.shuffled.tsv         — shuffled-baseline result per (array × shuffle)
  <out>.shuffled.summary.tsv — aggregate: fraction of shuffles with 10-15 bp sub-HOR signal
"""
import argparse, csv, os, random, re, statistics, subprocess, sys, tempfile

csv.field_size_limit(sys.maxsize)


def extract_consensuses(summary_tsv, chr_pattern, hor_min, hor_max):
    """Return list of (array_id, hor_period, hor_consensus, hor_autocorr, hor_mean_ed_tmpl, hor_cv)."""
    pat = re.compile(chr_pattern)
    out = []
    with open(summary_tsv) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            aid = row['array_id']
            if not pat.search(aid):
                continue
            try:
                hp = int(row.get('hor_period') or 0)
            except ValueError:
                continue
            if not (hor_min <= hp <= hor_max):
                continue
            cons = row.get('hor_consensus', '').strip()
            if not cons:
                continue
            out.append({
                'array_id':      aid,
                'hor_period':    hp,
                'hor_consensus': cons,
                'hor_autocorr':  float(row.get('hor_autocorr') or 0),
                'hor_n_monomers':int(row.get('hor_n_monomers') or 0),
                'hor_mean_ed_tmpl': float(row.get('hor_mean_ed_tmpl') or 0),
                'hor_cv':        float(row.get('hor_cv') or 0),
            })
    return out


def write_fasta(records, path):
    with open(path, 'w') as f:
        for r in records:
            f.write(f">{r['array_id']}\n{r['hor_consensus']}\n")


def shuffle_seq(rng, seq):
    chars = list(seq)
    rng.shuffle(chars)
    return ''.join(chars)


def make_tandem_fasta(seq, n_copies, out_path, header):
    """Tile the consensus n_copies times and wrap as a FASTA — emulates what
    ArraySplitter sees for a real array of HOR consensus repeats."""
    with open(out_path, 'w') as f:
        f.write(f">{header}\n")
        long_seq = seq * n_copies
        for i in range(0, len(long_seq), 70):
            f.write(long_seq[i:i + 70] + '\n')


def run_arraysplitter(asplit_bin, fasta_path, out_prefix, max_period=100000, threads=2):
    """Run ArraySplitter and return the dict from summary.tsv (single record)."""
    subprocess.run([asplit_bin, '-i', fasta_path, '-o', out_prefix,
                    '--max-period', str(max_period), '-t', str(threads),
                    '--method', 'autocorr'],
                   check=True, capture_output=True)
    summary_path = out_prefix + '.summary.tsv'
    with open(summary_path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            return row
    return None


def parse_period_classes(s):
    """Parse "13:29, 21:1023" → {13: 29, 21: 1023}."""
    out = {}
    for part in (s or '').split(','):
        part = part.strip()
        if not part or ':' not in part:
            continue
        k, v = part.split(':', 1)
        try:
            out[int(k.strip())] = int(v.strip())
        except ValueError:
            continue
    return out


def has_sub_period_in_range(period_classes_str, lo, hi):
    """True iff any period in 10-15 bp (inclusive) appears in the recursion descent."""
    cls = parse_period_classes(period_classes_str)
    return any(lo <= p <= hi for p in cls.keys())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--asplit-summary', required=True)
    ap.add_argument('--asplit-bin', required=True)
    ap.add_argument('--chr-pattern', default=r'chr(16|34|35)_(mat|pat)')
    ap.add_argument('--hor-period-min', type=int, default=180)
    ap.add_argument('--hor-period-max', type=int, default=210)
    ap.add_argument('--shuffles', type=int, default=10)
    ap.add_argument('--sub-period-lo', type=int, default=10)
    ap.add_argument('--sub-period-hi', type=int, default=15)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    consensuses = extract_consensuses(args.asplit_summary, args.chr_pattern,
                                       args.hor_period_min, args.hor_period_max)
    print(f'Extracted {len(consensuses)} consensus sequences matching {args.chr_pattern}', file=sys.stderr)

    # Write FASTA of consensuses
    fa_path = args.out_prefix + '.consensuses.fasta'
    write_fasta(consensuses, fa_path)
    print(f'Wrote {fa_path}', file=sys.stderr)

    # Per-array TSV (just the confidence-relevant fields)
    per_array_path = args.out_prefix + '.per_array.tsv'
    with open(per_array_path, 'w') as f:
        f.write('\t'.join(['array_id', 'hor_period', 'hor_autocorr', 'hor_n_monomers',
                            'hor_mean_ed_tmpl', 'hor_cv', 'confidence_score']) + '\n')
        for r in consensuses:
            # Crude confidence: autocorr × (1 − ed_per_bp), where ed_per_bp ≈ mean_ed_tmpl / hor_period
            ed_per_bp = r['hor_mean_ed_tmpl'] / max(r['hor_period'], 1)
            confidence = round(r['hor_autocorr'] * max(0.0, 1.0 - ed_per_bp), 3)
            f.write('\t'.join([r['array_id'], str(r['hor_period']), f"{r['hor_autocorr']:.4f}",
                                str(r['hor_n_monomers']), f"{r['hor_mean_ed_tmpl']:.2f}",
                                f"{r['hor_cv']:.4f}", str(confidence)]) + '\n')
    print(f'Wrote {per_array_path}', file=sys.stderr)

    # Shuffled-baseline control
    # For each consensus, run ArraySplitter on (a) the genuine consensus tiled 200×
    # and (b) 10 shuffled copies tiled 200×. Report whether a 10–15 bp sub-period
    # appears in the period_classes column.
    n_copies = 200  # gives ~38 kb tandem array, well above 10 kb
    shuf_path = args.out_prefix + '.shuffled.tsv'
    with open(shuf_path, 'w') as f:
        f.write('\t'.join(['array_id', 'replicate', 'is_shuffled',
                           'asplit_mono_period', 'asplit_hor_period',
                           'period_classes', 'has_sub_period_10_15']) + '\n')

        with tempfile.TemporaryDirectory() as td:
            for r in consensuses:
                aid = r['array_id']
                cons = r['hor_consensus']
                # genuine run
                fa = os.path.join(td, f'{aid}_genuine.fasta')
                make_tandem_fasta(cons, n_copies, fa, f'{aid}_genuine')
                op = os.path.join(td, f'{aid}_genuine_out')
                try:
                    row = run_arraysplitter(args.asplit_bin, fa, op)
                except subprocess.CalledProcessError as e:
                    row = None
                if row:
                    pcls = row.get('period_classes', '')
                    has = has_sub_period_in_range(pcls, args.sub_period_lo, args.sub_period_hi)
                    f.write('\t'.join([aid, '0', 'False',
                                       row.get('mono_period', ''),
                                       row.get('hor_period', ''),
                                       pcls, str(has)]) + '\n')
                # shuffled runs
                for k in range(args.shuffles):
                    s_seq = shuffle_seq(rng, cons)
                    fa = os.path.join(td, f'{aid}_s{k}.fasta')
                    make_tandem_fasta(s_seq, n_copies, fa, f'{aid}_s{k}')
                    op = os.path.join(td, f'{aid}_s{k}_out')
                    try:
                        row = run_arraysplitter(args.asplit_bin, fa, op)
                    except subprocess.CalledProcessError:
                        row = None
                    if row:
                        pcls = row.get('period_classes', '')
                        has = has_sub_period_in_range(pcls, args.sub_period_lo, args.sub_period_hi)
                        f.write('\t'.join([aid, str(k + 1), 'True',
                                           row.get('mono_period', ''),
                                           row.get('hor_period', ''),
                                           pcls, str(has)]) + '\n')
    print(f'Wrote {shuf_path}', file=sys.stderr)

    # Aggregate summary
    genuine_hits = 0
    shuffle_hits = 0
    shuffle_total = 0
    genuine_total = 0
    with open(shuf_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            is_shuf = row['is_shuffled'] == 'True'
            has = row['has_sub_period_10_15'] == 'True'
            if is_shuf:
                shuffle_total += 1
                if has:
                    shuffle_hits += 1
            else:
                genuine_total += 1
                if has:
                    genuine_hits += 1

    sum_path = args.out_prefix + '.shuffled.summary.tsv'
    with open(sum_path, 'w') as f:
        f.write('\t'.join(['category', 'hits', 'total', 'rate']) + '\n')
        f.write('\t'.join(['genuine_consensus', str(genuine_hits), str(genuine_total),
                           f'{genuine_hits / max(genuine_total, 1):.3f}']) + '\n')
        f.write('\t'.join(['composition_shuffled', str(shuffle_hits), str(shuffle_total),
                           f'{shuffle_hits / max(shuffle_total, 1):.3f}']) + '\n')

    print(f'Wrote {sum_path}', file=sys.stderr)
    print(f'\nGenuine consensus: {genuine_hits}/{genuine_total} have a 10-15 bp sub-period in recursive descent.')
    print(f'Composition-shuffled controls: {shuffle_hits}/{shuffle_total} have a 10-15 bp sub-period.')


if __name__ == '__main__':
    main()
