"""C1, C2: brute-force enumeration of the graphs of Q_n (S2 eq-count (5.11), eq-theta-rule (5.13),
cubic interaction vertices), independent of enumerate_q3.py.

Vertex 0 = O, vertices 1..2n+1 = interaction vertices.  For every p, every labelled forest F on
{1..2n+1} with (2n+1-p) edges (all edge subsets of K_{2n+1}, acyclicity by union-find), every
spoke multiset (multiplicities s_v >= 0, sum s_v = n+p, deg_F(v)+s_v <= 3): keep the graphs with
2|I_S| <= 3|V_S| - 4 for every vertex subset S with I_S nonempty (induced, with multiplicity) and
connected.  Iso classes with O fixed by networkx VF2 (MultiGraph, node colour O); automorphism
groups by brute force over all permutations of 1..2n+1.  Orbit check: #labelled = sum (2n+1)!/|Aut|.
"""
import itertools
import math
import sys
import time

import networkx as nx
import numpy as np


def forests(k, ne):
    pairs = list(itertools.combinations(range(1, k + 1), 2))
    for sub in itertools.combinations(pairs, ne):
        par = list(range(k + 1))

        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        ok = True
        for a, b in sub:
            ra, rb = f(a), f(b)
            if ra == rb:
                ok = False
                break
            par[ra] = rb
        if ok:
            yield sub


def spoke_multisets(cap, total):
    k = len(cap)

    def rec(i, left):
        if i == k:
            if left == 0:
                yield ()
            return
        for s in range(min(cap[i], left) + 1):
            for rest in rec(i + 1, left - s):
                yield (s,) + rest
    yield from rec(0, total)


def enumerate_n(n):
    k = 2 * n + 1
    nvtx = k + 1
    subsets = [S for S in range(1, 1 << nvtx)]
    sizes = np.array([bin(S).count('1') for S in subsets])
    # pair indicator over subsets
    pv = {}
    for a in range(nvtx):
        for b in range(a, nvtx):
            pv[(a, b)] = np.array([((S >> a) & 1) & ((S >> b) & 1) for S in subsets], dtype=np.int16)
    res = {}
    for p in range(1, k + 1):
        ne = k - p
        cand = []
        nfor = 0
        for F in forests(k, ne):
            nfor += 1
            deg = [0] * (k + 1)
            for a, b in F:
                deg[a] += 1
                deg[b] += 1
            if max(deg[1:]) > 3:
                continue
            cap = [3 - deg[v] for v in range(1, k + 1)]
            for s in spoke_multisets(cap, n + p):
                edges = list(F) + [(0, v + 1) for v in range(k) for _ in range(s[v])]
                cand.append(edges)
        good = []
        B = 20000
        for i0 in range(0, len(cand), B):
            batch = cand[i0:i0 + B]
            E = np.zeros((len(batch), len(subsets)), dtype=np.int16)
            for gi, edges in enumerate(batch):
                acc = np.zeros(len(subsets), dtype=np.int16)
                for a, b in edges:
                    acc = acc + pv[(min(a, b), max(a, b))]
                E[gi] = acc
            viol = ((E > 0) & (2 * E > 3 * sizes[None, :] - 4)).any(axis=1)
            for gi in np.nonzero(~viol)[0]:
                edges = batch[gi]
                G = nx.MultiGraph()
                G.add_nodes_from(range(nvtx))
                G.add_edges_from(edges)
                if nx.is_connected(G):
                    good.append(edges)
        res[p] = (nfor, len(cand), good)
        print('n=%d p=%d: forests %d, candidates (deg<=3) %d, passing count + connected %d'
              % (n, p, nfor, len(cand), len(good)), flush=True)
    return res


def mg(edges, nvtx):
    G = nx.MultiGraph()
    for v in range(nvtx):
        G.add_node(v, root=(v == 0))
    G.add_edges_from(edges)
    return G


def classes(goods, nvtx):
    reps = []
    counts = []
    nm = lambda a, b: a['root'] == b['root']
    for edges in goods:
        G = mg(edges, nvtx)
        for i, (H, _) in enumerate(reps):
            if nx.is_isomorphic(G, H, node_match=nm):
                counts[i] += 1
                break
        else:
            reps.append((G, edges))
            counts.append(1)
    return reps, counts


def aut_bruteforce(edges, nvtx):
    key = sorted(tuple(sorted(e)) for e in edges)
    auts = []
    for perm in itertools.permutations(range(1, nvtx)):
        mp = (0,) + perm
        k2 = sorted(tuple(sorted((mp[a], mp[b]))) for a, b in edges)
        if k2 == key:
            auts.append(mp)
    return auts


def skeleton_chains(edges, nvtx):
    G = mg(edges, nvtx)
    deg = dict(G.degree())
    branch = [v for v in range(nvtx) if deg[v] != 2]
    lengths = []
    used = set()
    for b in branch:
        for (u, w, kk) in G.edges(b, keys=True):
            e0 = tuple(sorted((u, w))) + (kk,)
            if e0 in used:
                continue
            L = 1
            used.add(e0)
            prev, cur = b, w
            while deg[cur] == 2:
                nxt = [(x, kk2) for (_, x, kk2) in G.edges(cur, keys=True) if tuple(sorted((cur, x))) + (kk2,) not in used]
                x, kk2 = nxt[0]
                used.add(tuple(sorted((cur, x))) + (kk2,))
                prev, cur = cur, x
                L += 1
            lengths.append(L)
    return sorted(lengths)


EXPLORER = {
    'G1': [('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), ('O', 5), (5, 6), (6, 7), (7, 'O')],
    'G2': [('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
    'G3': [('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
    'G4': [('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
    'G5': [('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
}
EXPL_AUT = {'G1': 2, 'G2': 2, 'G3': 1, 'G4': 8, 'G5': 2}
EXPL_OUT = {'G1': [1, 2, 4, 6, 7], 'G2': [1, 2, 4, 6, 7], 'G3': [1, 2, 5, 6, 7], 'G4': [1, 2, 4, 6, 7], 'G5': [1, 2, 5, 6, 7]}
EXPL_AUTGEN = {'G1': [{1: 7, 7: 1, 2: 6, 6: 2, 3: 5, 5: 3}], 'G2': [{6: 7, 7: 6}], 'G3': [],
               'G4': [{1: 2, 2: 1}, {6: 7, 7: 6}, {1: 6, 6: 1, 2: 7, 7: 2, 3: 5, 5: 3}], 'G5': [{1: 2, 2: 1}]}


def conv(E):
    return [(0 if a == 'O' else a, 0 if b == 'O' else b) for a, b in E]


if __name__ == '__main__':
    t0 = time.time()
    for n in (1, 2):
        res = enumerate_n(n)
        allg = [g for p in res for g in res[p][2]]
        reps, counts = classes(allg, 2 * n + 2)
        for (G, e), c in zip(reps, counts):
            au = aut_bruteforce(e, 2 * n + 2)
            print('  n=%d class: edges %s  degO %d  labelled %d  |Aut| %d  (2n+1)!/|Aut| %d  chains %s'
                  % (n, sorted(tuple(sorted(x)) for x in e), G.degree(0), c, len(au),
                     math.factorial(2 * n + 1) // len(au), skeleton_chains(e, 2 * n + 2)))
    n = 3
    res = enumerate_n(n)
    allg = [g for p in res for g in res[p][2]]
    print('n=3 passing graphs by p:', {p: len(res[p][2]) for p in res})
    reps, counts = classes(allg, 8)
    print('n=3: %d iso classes (O fixed)' % len(reps))
    nm = lambda a, b: a['root'] == b['root']
    tot = 0
    for (G, e), c in zip(reps, counts):
        au = aut_bruteforce(e, 8)
        tot += c
        deg = dict(G.degree())
        T = G.subgraph(range(1, 8))
        match = [k for k, E in EXPLORER.items() if nx.is_isomorphic(G, mg(conv(E), 8), node_match=nm)]
        print('  class ~ %s: degO %d, G-O tree %s, tree degrees %s, outputs %d, labelled %d, |Aut| %d, 7!/|Aut| %d, chains %s'
              % (match, deg[0], nx.is_tree(nx.Graph(T)) and T.number_of_edges() == 6,
                 sorted([deg[v] for v in range(1, 8)], reverse=True), sum(3 - deg[v] for v in range(1, 8)),
                 c, len(au), 5040 // len(au), skeleton_chains(e, 8)))
    print('sum labelled = %d' % tot)
    print('explorer labelled graphs, direct data:')
    for k, E in EXPLORER.items():
        e = conv(E)
        au = aut_bruteforce(e, 8)
        deg = {v: sum((a == v) + (b == v) for a, b in e) for v in range(8)}
        outs = [v for v in range(1, 8) if deg[v] == 2]
        gens_ok = all(tuple([0] + [g.get(v, v) for v in range(1, 8)]) in au for g in EXPL_AUTGEN[k])
        # group generated by stated generators
        grp = {tuple(range(8))}
        frontier = list(grp)
        gens = [tuple([0] + [g.get(v, v) for v in range(1, 8)]) for g in EXPL_AUTGEN[k]]
        while frontier:
            x = frontier.pop()
            for g in gens:
                y = tuple(g[x[i]] for i in range(8))
                if y not in grp:
                    grp.add(y)
                    frontier.append(y)
        print('  %s: |Aut| brute %d (stated %d), stated generators are automorphisms %s, generate group of order %d, '
              'outputs (deg 2) %s (stated %s), degO %d'
              % (k, len(au), EXPL_AUT[k], gens_ok, len(grp), outs, EXPL_OUT[k], deg[0]))
    print('time %.1fs' % (time.time() - t0))
