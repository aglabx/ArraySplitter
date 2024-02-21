#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import deque

import heapq

class WeightedValueHeap:
    def __init__(self):
        self.heap = []
    
    def insert(self, weight, value):
        heapq.heappush(self.heap, (weight, value))
    
    def pop(self):
        if self.heap:
            return heapq.heappop(self.heap)
        return None
    
    def peek(self):
        if self.heap:
            return self.heap[0]
        return None

def update(
    fs_tree,
    queue,
    nucl,
    fs_x,
    fs_xp,
    current_cid,
    cid,
    cutoff,
):
    if len(fs_x) > cutoff:
        fs_tree[cid] = (
            cid,
            [nucl],
            fs_x,
            fs_xp,
            None,
            [],
        )
        fs_tree[current_cid][5].append(cid)
        queue.append(cid)
        cid += 1
    return cid


def build_fs_tree_from_sequence(array, starting_seq_, names_, positions_, cutoff):
    fs_tree = {}
    queue = []
    cid = 0
    seq = starting_seq_
    names = names_[::]
    positions = positions_[::]
    children = []
    fs = (cid, [seq], names, positions, None, children)
    queue.append(cid)
    fs_tree[cid] = fs
    cid += 1
    cutoff = cutoff

    while queue:
        qcid = queue.pop()

        fs = fs_tree[qcid]

        fs_a, fs_c, fs_g, fs_t = [], [], [], []
        fs_ap, fs_cp, fs_gp, fs_tp = [], [], [], []

        for ii, pos in enumerate(fs[3]):
            current_cid = fs[0]
            seq = fs[1]
            name = fs[2][ii]

            if pos + 1 == len(array):
                continue

            nucl = array[pos + 1]
            if nucl == "A":
                fs_a.append(name)
                fs_ap.append(pos + 1)
            elif nucl == "C":
                fs_c.append(name)
                fs_cp.append(pos + 1)
            elif nucl == "G":
                fs_g.append(name)
                fs_gp.append(pos + 1)
            elif nucl == "T":
                fs_t.append(name)
                fs_tp.append(pos + 1)
            
        cid = update(
            fs_tree,
            queue,
            "A",
            fs_a,
            fs_ap,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            "C",
            fs_c,
            fs_cp,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            "G",
            fs_g,
            fs_gp,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            "T",
            fs_t,
            fs_tp,
            current_cid,
            cid,
            cutoff,
        )

    return fs_tree


def iter_fs_tree_from_sequence(array, starting_seq_, names_, positions_, cutoff):
    heap = []
    names = names_[::]
    positions = positions_[::]
    fs = (1, names, positions)
    heapq.heappush(heap, (1, fs))
    cutoff = cutoff
    while heap:
        L, fs = heapq.heappop(heap)

        fs_a, fs_c, fs_g, fs_t = [], [], [], []
        fs_ap, fs_cp, fs_gp, fs_tp = [], [], [], []

        for ii, pos in enumerate(fs[2]):
            
            name = fs[1][ii]

            if pos + 1 == len(array):
                continue

            nucl = array[pos + 1]
            if nucl == "A":
                fs_a.append(name)
                fs_ap.append(pos + 1)
            elif nucl == "C":
                fs_c.append(name)
                fs_cp.append(pos + 1)
            elif nucl == "G":
                fs_g.append(name)
                fs_gp.append(pos + 1)
            elif nucl == "T":
                fs_t.append(name)
                fs_tp.append(pos + 1)
            
        if len(fs_a) > cutoff:
            fs = (
                L+1,
                fs_a,
                fs_ap,
            )
            heapq.heappush(heap, (L+1, fs))
            yield fs

        if len(fs_c) > cutoff:
            fs = (
                L+1,
                fs_c,
                fs_cp,
            )
            heapq.heappush(heap, (L+1, fs))
            yield fs

        if len(fs_g) > cutoff:
            fs = (
                L+1,
                fs_g,
                fs_gp,
            )
            heapq.heappush(heap, (L+1, fs))
            yield fs

        if len(fs_t) > cutoff:
            fs = (
                L+1,
                fs_t,
                fs_tp,
            )
            heapq.heappush(heap, (L+1, fs))
            yield fs
