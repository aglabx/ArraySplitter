#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse

import os
from collections import Counter
from tqdm import tqdm

import editdistance as ed
from intervaltree import IntervalTree

from ArraySplitter.core_functions.io.fasta_reader import \
    sc_iter_fasta_file
from ArraySplitter.core_functions.io.satellome_reader import \
    sc_iter_satellome_file
from ArraySplitter.core_functions.io.trf_reader import sc_iter_trf_file
from ArraySplitter.core_functions.tools.fs_tree import \
    build_fs_tree_from_sequence


def get_top1_nucleotide(array):
    ### Step 1. Find the most frequent nucleotide (TODO: check all nucleotides and find with the best final score
    c = Counter()
    for n in "ACTG":
        c[n] = array.count(n)
        # print(n, array.count(n))
    return c.most_common(1)[0][0]


def get_fs_tree(array, top1_nucleotide, cutoff):
    ### Step 2. Build fs_tree (TODO:  optimize it for long sequences)
    names_ = [i for i in range(len(array)) if array[i] == top1_nucleotide]
    positions_ = names_[::]
    # print(f"Starting positions: {len(positions_)}")
    if not cutoff:
        if len(array) > 100_000:
            cutoff = 300
        else:
            cutoff = 10
    fs_tree, abrupted_nodes = build_fs_tree_from_sequence(
        array, top1_nucleotide, names_, positions_, cutoff
    )
    return fs_tree


def compute_hints(array, fs_tree, depth):
    ### Step 3. Find a list of hints (hint is the sequenece for array cutoff)
    cid1, seq1, names1, positions1, parents1, children1, cov1, loop1, extend1 = fs_tree[
        0
    ]

    cid2length = {}
    cid2length[0] = (1, len(names1))
    queue = [0]
    while queue:
        start = queue.pop(0)
        cid, seq, names, positions, parents, children, cov, loop, extend = fs_tree[
            start
        ]
        for child in children:
            cid2length[child] = (cid2length[cid][0] + 1, len(fs_tree[child][2]))
            queue.append(child)

    hints = []
    # print(cid2length)
    for length in range(depth):
        found = False
        c = Counter()
        tree = IntervalTree()
        for cid, (L, N) in cid2length.items():
            if L == length:
                found = True
                for start in fs_tree[cid][2]:
                    tree.addi(start, start + L)
                start = fs_tree[cid][2][0]
                end = fs_tree[cid][3][0]
                seq = array[start : end + 1]
                c[seq] = len(fs_tree[cid][2])

        tree.merge_overlaps()
        coverage = 0
        pred_sizes = Counter()
        pred_seqs = Counter()
        for start, end, _ in tree:
            coverage += abs(end - start)
            pred_sizes[abs(end - start)] += 1
            pred_seqs[array[start:end]] += 1
        rcoverage = 0
        for seq, cov in c.items():
            rcoverage += len(seq) * cov
        pcoverage = 0
        for seq, cov in pred_seqs.most_common(10):
            pcoverage += len(seq) * cov
        if found:
            # print(f"-{length}--{100.*pcoverage/len(array)}--{100.*coverage/len(array)}--{100.*rcoverage//len(array)}--")
            hints.append(c.most_common(500)[0])
    return hints


def compute_cuts(array, hints):
    ### Step 4. Find optimal cutoff
    best_cut_seq = None
    best_cut_score = 0.0
    best_period = None
    for sp, _ in hints:
        if len(sp) < 2:
            continue
        t = 0
        c = Counter()
        for x in array.split(sp):
            if len(x + sp) == 44:
                t += 1
            c[len(x + sp)] += 1
        a, b = c.most_common(1)[0]
        v = sum(c.values())
        score = b / v
        if score > best_cut_score:
            best_cut_score = score
            best_cut_seq = sp
            best_period = c.most_common(1)[0][0]
            # print(100*b//v, sp, t, c.most_common(3))
    return best_cut_seq, best_cut_score, best_period


### Step 5a. Try to cut long monomers to expected
def refine_repeat_even(repeat, best_period):
    if len(repeat) % best_period == 0:
        start = 0
        for _ in range(len(repeat) // best_period):
            yield repeat[start : start + best_period]
            start += best_period
    else:
        yield repeat


def decompose_array_iter1(array, best_cut_seq, best_period, verbose=True):
    t = 0
    repeats2count = Counter()
    decomposition = []
    cuts = array.split(best_cut_seq)
    for ii, x in enumerate(cuts):
        if ii + 1 != len(cuts):
            repeat = x + best_cut_seq
        else:
            repeat = x
        if not repeat:
            continue
        if len(repeat) == best_period:
            t += 1
            if verbose:
                print(len(repeat), repeat)
            repeats2count[repeat] += 1
            decomposition.append(repeat)
        else:
            # print(len(repeat), best_period)
            for repeat in refine_repeat_even(repeat, best_period):
                if len(repeat) == best_period:
                    t += 1
                if verbose:
                    print(len(repeat), repeat)
                repeats2count[repeat] += 1
                decomposition.append(repeat)
    if verbose:
        print(t)
    return decomposition, repeats2count


### Step 5b. Try to cut long monomers to expected
def refine_repeat_odd(repeat, best_period, most_common_monomer, verbose=False):
    if len(repeat) / best_period > 1.3:
        n = len(most_common_monomer)
        optimal_cut = 0
        best_ed = n
        for i in range(len(repeat) - n + 1):
            rep_b = repeat[i : i + n]
            dist = ed.eval(most_common_monomer, rep_b)
            if dist < best_ed:
                best_ed = dist
                optimal_cut = i
                if verbose:
                    print(
                        "Optimal cut",
                        best_ed,
                        optimal_cut,
                        len(repeat[:optimal_cut]),
                        len(repeat[optimal_cut:]),
                    )
        if best_ed < n / 2:
            if optimal_cut == 0:
                optimal_cut += n
            a = repeat[:optimal_cut]
            b = repeat[optimal_cut:]
            if min(len(a), len(b)) < 0:  # n/3:
                yield repeat
            else:
                if a:
                    yield a
                if b:
                    yield b
        else:
            yield repeat
    else:
        yield repeat


def decompose_array_iter2(decomposition, best_period, repeats2count_ref, verbose=True):
    t = 0
    repeats2count = Counter()
    refined_decomposition = []
    most_common_monomer = None
    for monomer, tf in repeats2count_ref.most_common(1000):
        if len(monomer) == best_period:
            most_common_monomer = monomer
            break
    assert most_common_monomer
    for repeat in decomposition:
        if verbose:
            print("Repeat under consideration", len(repeat), repeat)
        for repeat in refine_repeat_odd(
            repeat, best_period, most_common_monomer, verbose=verbose
        ):
            if len(repeat) == best_period:
                t += 1
            if verbose:
                print("Added:", len(repeat), repeat)
            repeats2count[repeat] += 1
            refined_decomposition.append(repeat)
    if verbose:
        print(t)
    return (
        refined_decomposition,
        repeats2count,
        len(refined_decomposition) != len(decomposition),
    )


def print_monomers(decomposition, repeats2count, best_period):
    start2tf = Counter()
    for monomer in decomposition:
        start2tf[monomer[:5]] += 1
    print(start2tf)

    most_common_monomer = None
    for monomer, tf in repeats2count.most_common(1000):
        if len(monomer) == best_period:
            most_common_monomer = monomer
            break
    assert most_common_monomer
    for repeat in decomposition:
        print(
            len(repeat),
            start2tf[repeat[:5]],
            repeat,
            ed.eval(repeat, most_common_monomer),
        )


def print_pause_clean(decomposition, repeats2count, best_period):
    print_monomers(decomposition, repeats2count, best_period)
    input("?")


#   clear_output(wait=True)


def decompose_array(array, depth=500, cutoff=20, verbose=True):
    ### Step 1. Find the most frequent nucleotide (TODO: check all nucleotides and find the one with the best final score)
    top1_nucleotide = get_top1_nucleotide(array)
    # print("top1_nucleotide:", top1_nucleotide)
    # top1_nucleotide = "A"
    ### Step 2. Build fs_tree (TODO: optimize it for long sequences)
    fs_tree = get_fs_tree(array, top1_nucleotide, cutoff=cutoff)
    ### Step 3. Find a list of hints (hint is the sequence for array cutting)
    ### Here I defined it as a sequence with maximal coverage in the original array for a given length
    hints = compute_hints(array, fs_tree, depth)
    ### Step 4. Find the optimal cut sequence and the best period
    ### Defined as the maximal fraction of the cut sequence to the total cut sequence
    best_cut_seq, best_cut_score, best_period = compute_cuts(array, hints)

    ### Step 5. Cut the array
    ### The first iteration finds monomer frequencies
    decomposition, repeats2count = decompose_array_iter1(
        array, best_cut_seq, best_period, verbose=False
    )
    
    assert "".join(decomposition) == array
    ### The second iteration tries to cut longer monomers to the expected length
    changed = True
    while changed:
        decomposition, repeats2count, changed = decompose_array_iter2(
            decomposition, best_period, repeats2count, verbose=False
        )
        assert "".join(decomposition) == array

    ### TODO: The third iteration tries to glue short dangling monomers to the nearest monomer

    return decomposition, repeats2count, best_cut_seq, best_cut_score, best_period


def get_array_generator(input_file, format):
    '''Get array generator by format.'''
    if format == "fasta":
        return sc_iter_fasta_file(input_file)
    if format == "trf":
        return sc_iter_trf_file(input_file)
    if format == "satellome":
        return sc_iter_satellome_file(input_file)
    
    print(f"Unknown format: {format}")
    exit(1)
    

def main(input_file, output_prefix, format, threads):
    """Main function."""

    sequences = get_array_generator(input_file, format)
    total = 0
    for _ in sequences:
        total += 1
    sequences = get_array_generator(input_file, format)

    print(f"Start processing")
    print(f"Output prefix: {output_prefix}")

    depth = 500
    cutoff = 20
    verbose = False

    output_file = f"{output_prefix}.decomposed.fasta"

    with open(output_file, "w") as fw:
        for header, array in tqdm(sequences, total=total):
            # print(len(array), end=" ")
            (
                decomposition,
                repeats2count,
                best_cut_seq,
                best_cut_score,
                best_period,
            ) = decompose_array(array, depth=depth, cutoff=cutoff, verbose=verbose)

            # print("best period:", best_period, "len:", len(decomposition))
            # print_pause_clean(decomposition, repeats2count, best_period)

            fw.write(f">{header} {best_period}\n")
            fw.write(" ".join(decomposition) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="De novo decomposition of satellite DNA arrays into monomers"
    )
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument(
        "--format",
        help="Input format: fasta, trf [fasta]",
        required=False,
        default="fasta",
    )
    parser.add_argument("-o", "--output", help="Output prefix", required=True)
    parser.add_argument(
        "-t", "--threads", help="Number of threads", required=False, default=4
    )
    args = parser.parse_args()

    input_file = args.input
    output_prefix = args.output
    format = args.format
    threads = int(args.threads)

    if not os.path.isfile(input_file):
        print(f"File {input_file} not found")
        exit(1)

    main(input_file, output_prefix, format, threads)
