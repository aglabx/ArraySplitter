"""Null control of the paper's 'period-aware cut agreement' on the REAL zebra finch panel.
Comparator = uniformly random cut positions (same number of cuts as ArraySplitter made), scored with the
paper's own functions from scripts/benchmark/compare_cuts_centroanno.py."""
import sys, random, importlib.util, csv, collections
PREFIX, REPO = sys.argv[1], sys.argv[2]   # output prefix of an arraysplitter run, repo path
spec = importlib.util.spec_from_file_location('cmp', REPO + '/scripts/benchmark/compare_cuts_centroanno.py')
cmp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cmp)
asp = cmp.parse_asplit(PREFIX + '.summary.tsv', PREFIX + '.hors.tsv')
L = {r['array_id']: int(r['array_length']) for r in csv.DictReader(open(PREFIX + '.summary.tsv'), delimiter='\t')}
rng = random.Random(1)
pool = collections.defaultdict(list); per_class = collections.defaultdict(lambda: collections.defaultdict(list)); ncuts = collections.Counter()
for aid, a in asp.items():
    P, cuts = a['period'], a['cuts']
    if not cuts or P <= 0: continue
    sf = sorted(x % P for x in cuts)
    rnd = [rng.randrange(0, L[aid]) for _ in cuts]                       # nonsense tool: random cuts
    shift = rng.randrange(1, P) if P > 1 else 0
    rot = [c + shift for c in cuts]                                      # perfect tool, different rotation
    cls = '5-7 bp' if P <= 7 else ('8-200 bp' if P <= 200 else ('201-1000 bp' if P <= 1000 else '>1000 bp'))
    ncuts[cls] += len(cuts)
    for label, other in (('random', rnd), ('rotated-perfect', rot)):
        d = [cmp.nearest_phase_delta_fast(c, sf, P) for c in other]
        pool[label] += d; per_class[cls][label] += d
def pct(d): return tuple(round(100 * sum(1 for x in d if abs(x) <= t) / len(d), 1) for t in (0, 10, 100))
tot = sum(ncuts.values())
print(f"arrays scored: {len(asp)}, ArraySplitter cuts pooled: {tot}")
print("share of pooled cuts by top-level period class:", {k: f"{100*v/tot:.1f}%" for k, v in ncuts.items()})
print("\n'agreement' with ArraySplitter at +-0 / +-10 / +-100 bp under the paper's period-aware folding:")
print(f"   RANDOM cut positions (pooled)              : {pct(pool['random'])}")
print(f"   ArraySplitter's own cuts, rotated (pooled) : {pct(pool['rotated-perfect'])}")
print("   paper, FasTAN vs ArraySplitter             : (85.6, 95.2, 99.7)")
for cls in ('5-7 bp', '8-200 bp', '201-1000 bp', '>1000 bp'):
    if per_class[cls]['random']:
        print(f"   {cls:12s} random {str(pct(per_class[cls]['random'])):22s} rotated-perfect {pct(per_class[cls]['rotated-perfect'])}")
