import csv, sys, collections
csv.field_size_limit(10**9)
prefix = sys.argv[1]
hors = {}
types = {}
for r in csv.DictReader(open(prefix + '.hors.tsv'), delimiter='\t'):
    if r['type'] in ('monomer', 'flank', 'sub_hor'):
        hors[(r['array_id'], int(r['level']), int(r['idx']))] = r['sequence']
        types[(r['array_id'], int(r['level']), int(r['idx']))] = r['type']
stat = collections.defaultdict(collections.Counter)
for r in csv.DictReader(open(prefix + '.monomers.tsv'), delimiter='\t'):
    if r['type'] not in ('base_monomer', 'monomer'): continue
    aid = r['array_id']; lvl = int(r['parent_level']); pidx = int(r['parent_idx']); seq = r['sequence']
    key = (aid, lvl, pidx)
    if key not in hors:
        stat[aid]['parent_row_missing'] += 1; continue
    ok = seq in hors[key]
    if ok and types[key] != 'flank':
        stat[aid]['join_ok'] += 1
    else:
        # would the join work with +1 (skipping a left flank)?
        alt = (aid, lvl, pidx + 1)
        if lvl == 1 and alt in hors and seq in hors[alt]:
            stat[aid]['join_WRONG_but_idx+1_is_right' + ('(points at flank row)' if types[key]=='flank' else '')] += 1
        else:
            stat[aid]['join_WRONG' + ('(points at flank row)' if types[key]=='flank' else '')] += 1
tot = collections.Counter()
for aid in sorted(stat):
    print(f"{aid[:42]:42s}", dict(stat[aid]))
    tot.update({k.split('(')[0]: v for k, v in stat[aid].items()})
print('TOTAL', dict(tot))
