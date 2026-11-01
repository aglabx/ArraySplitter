#!/usr/bin/env python3
"""Do the children of every split unit tile it? (stronger than check_join.py)

    python3 -I check_tiling.py PREFIX [PREFIX ...]

PREFIX = output prefix of an arraysplitter run (PREFIX.hors.tsv, PREFIX.monomers.tsv).
check_join.py only asks whether a child's sequence occurs somewhere in its parent, which in a repeat does not
show that it was assigned to the right copy or position. Here, for every parent unit (hors.tsv rows of type `monomer` or `sub_hor`) that has children (hors.tsv rows one level down, and monomers.tsv rows pointing to it),
the parent counts as tiled when its children, in some order, concatenate exactly to it (walked from position 0,
each step must match one unused child; sub-HORs and base monomers of one parent live in different files, so output
order is not used). Prints per prefix: parents, tiled, not tiled (named, first 10), by total child length.
Exit 1 if a file is missing; the counts are the result, not a pass/fail gate.
"""
import csv
import sys
sys.setrecursionlimit(100000)
from collections import defaultdict

csv.field_size_limit(1 << 30)


def tiles(parent, kids, budget=200000):
    """True / False if the children, in some order, do / do not concatenate exactly to the parent; None if the
    backtracking search used up its step budget (reported as undecided, never as tiled)."""
    left = defaultdict(int)
    for s in kids:
        left[s] += 1
    if sum(len(s) * c for s, c in left.items()) != len(parent):
        return False
    lengths = sorted({len(s) for s in kids}, reverse=True)
    steps = [0]

    def go(cur):
        if cur == len(parent):
            return True
        steps[0] += 1
        if steps[0] > budget:
            raise TimeoutError
        for L in lengths:
            s = parent[cur:cur + L]
            if len(s) == L and left.get(s, 0) > 0:
                left[s] -= 1
                if go(cur + L):
                    return True
                left[s] += 1
        return False

    try:
        return go(0)
    except TimeoutError:
        return None


def check(prefix):
    parents = {}
    children = defaultdict(list)  # (array, level, idx) -> [seq]
    with open(prefix + ".hors.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lev = r["level"]
            if not lev.isdigit():
                continue
            if r["type"] in ("monomer", "sub_hor"):
                parents[(r["array_id"], int(lev), r["idx"])] = r["sequence"]
            if int(lev) >= 2 and r["parent_idx"] not in ("-", "NA", "", "-1") and r["type"] in ("sub_hor", "monomer", "flank"):
                children[(r["array_id"], int(lev) - 1, r["parent_idx"])].append(r["sequence"])
    with open(prefix + ".monomers.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["parent_idx"] in ("-", "NA", "") or not r["parent_level"].isdigit():
                continue
            children[(r["array_id"], int(r["parent_level"]), r["parent_idx"])].append(r["sequence"])
    n = tiled = 0
    bad = defaultdict(list)
    for key, kids in children.items():
        if key not in parents:
            bad["parent row missing"].append(key)
            continue
        ps = parents[key]
        n += 1
        res = tiles(ps, kids)
        if res:
            tiled += 1
            continue
        if res is None:
            bad["undecided (search budget)"].append(key)
            continue
        tot = sum(len(s) for s in kids)
        reason = ("children shorter than the parent" if tot < len(ps) else
                  "children longer than the parent" if tot > len(ps) else "same length, no exact tiling")
        bad[reason].append(key)
    print(f"{prefix}: parents with children {n}, tiled exactly {tiled} ({100 * tiled / max(1, n):.2f}%)")
    for reason, keys in sorted(bad.items(), key=lambda x: -len(x[1])):
        ex = ", ".join(f"{a.split('|')[0]}:L{l}:{i}" for a, l, i in keys[:10])
        print(f"  {reason}: {len(keys)}  e.g. {ex}")


for p in sys.argv[1:]:
    check(p)
