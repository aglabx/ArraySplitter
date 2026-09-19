#!/usr/bin/env python3
"""Minimal tandem-array simulator with exact truth (monomer boundaries)."""
import random

def rand_seq(n, rng, gc=0.5):
    w = [(1-gc)/2, gc/2, gc/2, (1-gc)/2]
    return ''.join(rng.choices('ACGT', weights=w, k=n))

def mutate(seq, rng, sub=0.0, indel=0.0):
    out = []
    for ch in seq:
        r = rng.random()
        if r < indel/2:            # deletion
            continue
        if r < indel:              # insertion after
            out.append(ch); out.append(rng.choice('ACGT')); continue
        if rng.random() < sub:
            out.append(rng.choice([c for c in 'ACGT' if c != ch]))
        else:
            out.append(ch)
    return ''.join(out)

def flat_array(mono_len, n_copies, rng, sub=0.02, indel=0.0, gc=0.5, flank=100):
    """Star-phylogeny flat array: every copy = independently mutated consensus."""
    cons = rand_seq(mono_len, rng, gc)
    parts, bounds, pos = [], [], flank
    for _ in range(n_copies):
        m = mutate(cons, rng, sub, indel)
        bounds.append(pos); pos += len(m); parts.append(m)
    bounds.append(pos)
    seq = rand_seq(flank, rng) + ''.join(parts) + rand_seq(flank, rng)
    return seq, bounds

def hor_array(mono_len, k, n_hor, rng, between=0.2, within=0.01, indel=0.0, gc=0.5, flank=100):
    """HOR: k monomers each diverged `between` from an ancestral monomer; HOR copies mutated by `within`."""
    anc = rand_seq(mono_len, rng, gc)
    hor_monos = [mutate(anc, rng, between, 0.0) for _ in range(k)]
    parts, mono_bounds, hor_bounds, pos = [], [], [], flank
    for _ in range(n_hor):
        hor_bounds.append(pos)
        for m in hor_monos:
            mm = mutate(m, rng, within, indel)
            mono_bounds.append(pos); pos += len(mm); parts.append(mm)
    mono_bounds.append(pos); hor_bounds.append(pos)
    seq = rand_seq(flank, rng) + ''.join(parts) + rand_seq(flank, rng)
    return seq, mono_bounds, hor_bounds

def write_fasta(path, recs):
    with open(path, 'w') as f:
        for name, seq in recs:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
