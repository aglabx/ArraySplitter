#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 02.02.2024
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


def update(
    fs_tree,
    queue,
    abrupted_nodes,
    nucl,
    fs_x,
    fs_xp,
    loop,
    extend,
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
            [current_cid],
            [],
            [len(fs_xp)],
            loop,
            extend,
        )
        fs_tree[current_cid][5].append(cid)
        queue.append(cid)
        cid += 1
    else:
        for ii, fs in enumerate(fs_x):
            abrupted_nodes[fs] = (fs_xp[ii], nucl)
    return cid


def build_fs_tree_from_sequence(array, starting_seq_, names_, positions_, cutoff):
    fs_tree = {}
    abrupted_nodes = {}
    queue = []
    cid = 0
    loop = None
    seq = starting_seq_
    names = names_[::]
    positions = positions_[::]
    parents = []
    children = []
    extend = True
    coverage = [len(positions)]
    fs = (cid, [seq], names, positions, parents, children, coverage, loop, extend)
    queue.append(cid)
    fs_tree[cid] = fs
    cid += 1
    cutoff = cutoff

    while queue:
        qcid = queue.pop()

        fs = fs_tree[qcid]

        assert qcid == fs[0]

        if not fs[-1]:
            continue

        fs_a, fs_c, fs_g, fs_t, fs_n, fs_start_n = [], [], [], [], [], []
        fs_ap, fs_cp, fs_gp, fs_tp, fs_np, fs_start_p = [], [], [], [], [], []

        for ii, pos in enumerate(fs[3]):
            current_cid = fs[0]
            seq = fs[1]
            name = fs[2][ii]
            loop = None
            extend = True

            if pos + 1 == len(array):
                fs_n.append(name)
                fs_np.append(pos + 1)
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
            else:
                fs_n.append(name)
                fs_np.append(pos + 1)

        if len(fs_start_n) > cutoff:
            fs_tree[0][4].append(current_cid)
            fs_tree[current_cid][5].append(0)
        if fs_start_n:
            for ii, fs in enumerate(fs_start_n):
                abrupted_nodes[fs] = (fs_start_p[ii], "L")

        cid = update(
            fs_tree,
            queue,
            abrupted_nodes,
            "A",
            fs_a,
            fs_ap,
            loop,
            extend,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            abrupted_nodes,
            "C",
            fs_c,
            fs_cp,
            loop,
            extend,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            abrupted_nodes,
            "G",
            fs_g,
            fs_gp,
            loop,
            extend,
            current_cid,
            cid,
            cutoff,
        )
        cid = update(
            fs_tree,
            queue,
            abrupted_nodes,
            "T",
            fs_t,
            fs_tp,
            loop,
            extend,
            current_cid,
            cid,
            cutoff,
        )

        for ii, fs in enumerate(fs_n):
            abrupted_nodes[fs] = (fs_np[ii], "N")

    return fs_tree, abrupted_nodes
