#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import deque

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
    
    queue = deque()
    cid = 0
    seq = starting_seq_
    names = names_[::]
    positions = positions_[::]
    fs = (cid, [seq], names, positions)
    queue.appendleft(fs)
    cid += 1
    cutoff = cutoff
    while queue:
        fs = queue.popleft()

        fs_a, fs_c, fs_g, fs_t = [], [], [], []
        fs_ap, fs_cp, fs_gp, fs_tp = [], [], [], []

        for ii, pos in enumerate(fs[3]):
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
            
        if len(fs_a) > cutoff:
            fs = (
                cid,
                ["A"],
                fs_a,
                fs_ap,
            )
            queue.appendleft(fs)
            cid += 1
            yield fs

        if len(fs_c) > cutoff:
            fs = (
                cid,
                ["C"],
                fs_c,
                fs_cp,
            )
            queue.appendleft(fs)
            cid += 1
            yield fs

        if len(fs_g) > cutoff:
            fs = (
                cid,
                ["G"],
                fs_g,
                fs_gp,
            )
            queue.appendleft(fs)
            cid += 1
            yield fs

        if len(fs_t) > cutoff:
            fs = (
                cid,
                ["T"],
                fs_t,
                fs_tp,
            )
            queue.appendleft(fs)
            cid += 1
            yield fs
        