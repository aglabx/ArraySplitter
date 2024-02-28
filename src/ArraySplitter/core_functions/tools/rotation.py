#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 20.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from .sequences import get_revcomp
from collections import Counter
from tqdm import tqdm

def best_kmer_for_start(arrays, k=5):
    k = 5
    kmers = Counter()
    for array in arrays:
        repeat = "".join(array)
        prev_kmer = None
        for i in range(len(repeat)-k+1):
            kmer = repeat[i:i+k]
            if kmer == prev_kmer:
                continue
            prev_kmer = kmer
            kmers[kmer] += 1

    best_kmer = None
    max_f1 = 0
    for kmer, _ in tqdm(kmers.most_common(100)):
        copy_number = []
        for array in arrays:
            for repeat in array:
                copy_number.append(repeat.count(kmer))

        n = len(copy_number)
        r = len([x for x in copy_number if x > 0])/n
        p = len([x for x in copy_number if x == 1])/n
        f1 = 2 * (r * p) / (r + p)
        if f1 < 0.5:
            continue
        if f1 > max_f1:
            best_kmer = (r, p, f1, kmer)
            max_f1 = f1
    return max_f1, best_kmer

def rotate_arrays(arrays_, starting_kmer=None):

    if not starting_kmer:
        arrays = []
        for array in arrays_:
            feature1 = array.count("C") < array.count("G")
            feature2 = array.count("A") < array.count("T")
            if feature1 and feature2 or (not feature1 and not feature2):
                if feature1 and feature2:
                    arrays.append(array.split())
                else:
                    arrays.append(get_revcomp(array).split())
            else:
                if feature1 and not feature2:
                    arrays.append(get_revcomp(array).split())
                else:
                    arrays.append(array.split())
        max_f1, (r, p, f1, best_kmer) = best_kmer_for_start(arrays)
    else:
        best_kmer = starting_kmer
        best_rev_kmer = get_revcomp(best_kmer)
        arrays = []
        for array in arrays_:
            forward_feature = array.count(best_rev_kmer) < array.count(best_kmer)
            if forward_feature:
                arrays.append(array.split())
            else:
                arrays.append(get_revcomp(array).split())
        
    new_arrays = []
    for array in arrays:
        new_array = []
        prev_repeat = ""
        queue = array[::]
        while queue:
            repeat = prev_repeat + queue.pop(0)
            pos = repeat.find(best_kmer, 1)
            if pos == -1:
                prev_repeat = repeat
                continue
            new_array.append(repeat[:pos])
            prev_repeat = repeat[pos:]
        if prev_repeat:
            while True:
                pos = prev_repeat.find(best_kmer, 1)
                if pos == -1:
                    new_array.append(prev_repeat)
                    break
                new_array.append(prev_repeat[:pos])
                prev_repeat = prev_repeat[pos:]
            
        assert "".join(array) == "".join(new_array)
        
        new_arrays.append(new_array)

    return new_arrays
