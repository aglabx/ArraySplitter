import csv, collections, sys, bisect
csv.field_size_limit(10**9)

def load(prefix):
    """Return dict aid -> dict(summary=row, top=[(start,end,type)], base=[(start,end)], n_subhor=int, max_level=int)
    Coordinates are in the INPUT orientation (flipped back if the tool reverse-complemented)."""
    summ = {r['array_id']: r for r in csv.DictReader(open(prefix + '.summary.tsv'), delimiter='\t')}
    top = collections.defaultdict(list); sub = collections.Counter(); maxlvl = collections.Counter()
    for r in csv.DictReader(open(prefix + '.hors.tsv'), delimiter='\t'):
        if r['type'] in ('monomer', 'flank') and r['level'] == '1':
            top[r['array_id']].append((int(r['idx']), r['type'], int(r['length'])))
        elif r['type'] == 'sub_hor':
            sub[r['array_id']] += 1
            maxlvl[r['array_id']] = max(maxlvl[r['array_id']], int(r['level']))
    base = collections.defaultdict(list)
    for r in csv.DictReader(open(prefix + '.monomers.tsv'), delimiter='\t'):
        if r['type'] in ('base_monomer', 'monomer'):
            base[r['array_id']].append((int(r['idx']), int(r['length'])))
    out = {}
    for aid, s in summ.items():
        L = int(s['array_length']); rev = s['orientation'] == 'rev'
        t = sorted(top[aid]); b = sorted(base[aid])
        pos = 0; bi = 0; tops = []; bases = []
        for idx, typ, ln in t:
            tops.append((pos, pos + ln, typ))
            if typ == 'monomer':
                c = 0
                while c < ln and bi < len(b):
                    l = b[bi][1]; bases.append((pos + c, pos + c + l)); c += l; bi += 1
            pos += ln
        if rev:
            tops = sorted((L - e, L - s_, ty) for s_, e, ty in tops)
            bases = sorted((L - e, L - s_) for s_, e in bases)
        out[aid] = dict(summary=s, top=tops, base=bases, n_subhor=sub[aid], max_level=max(1, maxlvl[aid]), L=L)
    return out

def boundary_score(pred_cuts, true_bounds, P, tol=2):
    """Rotation-tolerant boundary score. For each predicted cut, delta to nearest true boundary;
    modal delta r = rotation. TP = predicted cuts with |delta - r| <= tol (circular in P)."""
    tb = sorted(true_bounds)
    lo, hi = tb[0], tb[-1]
    pc = sorted(set(c for c in pred_cuts if lo - P <= c <= hi + P))
    if not pc: return dict(r=None, prec=0.0, rec=0.0, f1=0.0, n_pred=0)
    deltas = []
    for c in pc:
        i = bisect.bisect_left(tb, c)
        cand = [tb[j] for j in (i - 1, i) if 0 <= j < len(tb)]
        t = min(cand, key=lambda x: abs(x - c))
        deltas.append(c - t)
    cnt = collections.Counter(d % P for d in deltas)
    # best rotation with tolerance window
    best_r, best = 0, -1
    for r in cnt:
        tot = sum(cnt[(r + k) % P] for k in range(-tol, tol + 1))
        if tot > best: best, best_r = tot, r
    def ok(d): 
        x = (d - best_r) % P
        return min(x, P - x) <= tol
    tp = sum(1 for d in deltas if ok(d))
    # recall: true boundaries (shifted by r) that have a predicted cut within tol
    pcs = pc
    rec_hits = 0; n_true = 0
    r_signed = best_r if best_r <= P // 2 else best_r - P
    for t in tb:
        target = t + r_signed
        if target < lo or target > hi: continue
        n_true += 1
        i = bisect.bisect_left(pcs, target - tol)
        if i < len(pcs) and pcs[i] <= target + tol: rec_hits += 1
    prec = tp / len(pc); rec = rec_hits / max(1, n_true)
    f1 = 2 * prec * rec / max(1e-9, prec + rec)
    return dict(r=r_signed, prec=prec, rec=rec, f1=f1, n_pred=len(pc))
