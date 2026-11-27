"""Rooted Laman graphs of Q_n on C^H with cubic vertices.

Counting [drv]: |I| = Hn+1, |V| = (H-1)n+2, (H-1)|I_S| <= H|V_S| - H - 1.  With Gamma \\ O having c components,
b loops and s spokes: summing the subgraph bounds over the components of Gamma \\ O gives c = 1, s = n + 1 - b.
  tree      (b = 0): s = n+1 spokes, V_b = (H-1)n+1 bulk vertices, condition per connected subtree with s' marks:
            |subtree| >= (H-1)(s'-1)+1; H = 1: one vertex with n+1 parallel spokes (banana).
  unicyclic (b = 1): s = n spokes, all subgraph bounds checked by brute force.
Every unmarked leaf violates the bound (removing it), so leaves are marked; outputs at v = 3 - deg(v).
Labelled form: bulk vertices 0..V_b-1, root V_b; bulk edges first (oriented away from a BFS root), spokes (v, root)
in vertex order.  |Aut| = automorphisms of Gamma fixing O.
"""
import itertools
from functools import lru_cache

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def laman_ok(nv, edges, H):
    """brute-force: (H-1)|I_S| <= H|V_S| - H - 1 for every vertex subset with |V_S| >= 2 (induced)."""
    for k in range(2, nv + 1):
        for S in itertools.combinations(range(nv), k):
            Ss = set(S)
            m = sum(1 for a, b in edges if a in Ss and b in Ss)
            if m == 0:
                continue
            if (H - 1) * m > H * k - H - 1:
                return False
    return True


def _aut_order(G, marks):
    GM = GraphMatcher(G, G, node_match=lambda a, b: a['m'] == b['m'])
    return sum(1 for _ in GM.isomorphisms_iter())


def _label(G, marks):
    """BFS labelling from a marked vertex of minimal degree; tree/bulk edges oriented parent -> child."""
    start = sorted(marks, key=lambda v: (G.degree(v), v))[0]
    order = list(nx.bfs_tree(G, start).nodes())
    lab = {v: i for i, v in enumerate(order)}
    bulk = []
    seen = set()
    for u, w in nx.bfs_edges(G, start):
        bulk.append((lab[u], lab[w]))
        seen.add(frozenset((u, w)))
    extra = [(lab[u], lab[w]) for u, w in G.edges() if frozenset((u, w)) not in seen]
    extra = [tuple(sorted(e)) for e in extra]
    nb = G.number_of_nodes()
    spokes = [(lab[v], nb) for v in sorted(marks, key=lambda x: lab[x])]
    return nb, bulk + sorted(extra) + spokes


def _marked_classes(G, nmarks, maxdeg=3):
    """orbit representatives of mark sets: leaves marked, marked vertices of degree <= maxdeg-1."""
    leaves = [v for v in G if G.degree(v) <= 1]
    cand = [v for v in G if G.degree(v) == 2 and v not in leaves] if maxdeg == 3 else []
    if any(G.degree(v) > maxdeg - 1 for v in leaves):
        return []
    extra = nmarks - len(leaves)
    if extra < 0 or extra > len(cand):
        return []
    auts = list(GraphMatcher(G, G).isomorphisms_iter())
    remaining = {frozenset(leaves) | frozenset(c) for c in itertools.combinations(cand, extra)}
    reps = []
    while remaining:
        ms = min(remaining, key=lambda s: sorted(s))
        orbit = {frozenset(a[v] for v in ms) for a in auts}
        remaining -= orbit
        reps.append(ms)
    return reps


def rooted_trees(H, n):
    """list of dicts: nb, edges (labelled), aut, kind='tree'."""
    if n == 0:
        return [dict(nb=1, edges=[(0, 1)], aut=1, kind='Q0')]
    if H == 1:
        if n + 1 > 3:
            return []
        return [dict(nb=1, edges=[(0, 1)] * (n + 1), aut=1, kind='banana')]  # aut: see note in kernels
    V = (H - 1) * n + 1
    trees = [nx.empty_graph(1)] if V == 1 else list(nx.nonisomorphic_trees(V))
    out = []
    for T in trees:
        if max(dict(T.degree()).values(), default=0) > 3:
            continue
        for ms in _marked_classes(T, n + 1):
            ok = True
            # subtree condition
            for k in range(2, V):
                if not ok:
                    break
                for S in itertools.combinations(T.nodes(), k):
                    if not nx.is_connected(T.subgraph(S)):
                        continue
                    s = len(ms & set(S))
                    if s and k < (H - 1) * (s - 1) + 1:
                        ok = False
                        break
            if not ok:
                continue
            G = T.copy()
            nx.set_node_attributes(G, {v: (v in ms) for v in G}, 'm')
            nb, edges = _label(G, ms)
            out.append(dict(nb=nb, edges=edges, aut=_aut_order(G, ms), kind='tree'))
    return out


def rooted_unicyclic(H, n):
    """Gamma \\ O connected unicyclic, n spokes."""
    if H < 3 or n < 1:
        return []
    V = (H - 1) * n + 1
    seen = []
    for T in nx.nonisomorphic_trees(V):
        for u, w in itertools.combinations(T.nodes(), 2):
            if T.has_edge(u, w):
                continue
            G = T.copy()
            G.add_edge(u, w)
            if max(dict(G.degree()).values()) > 3:
                continue
            cyc = nx.cycle_basis(G)[0]
            if len(cyc) < H + 1:
                continue
            if any(nx.is_isomorphic(G, X) for X in seen):
                continue
            seen.append(G)
    out = []
    for G0 in seen:
        for ms in _marked_classes(G0, n):
            G = G0.copy()
            nx.set_node_attributes(G, {v: (v in ms) for v in G}, 'm')
            nb, edges = _label(G, ms)
            if not laman_ok(nb + 1, edges, H):
                continue
            out.append(dict(nb=nb, edges=edges, aut=_aut_order(G, ms), kind='unicyclic'))
    return out


def outputs(nb, edges):
    deg = [0] * nb
    for a, b in edges:
        if a < nb:
            deg[a] += 1
        if b < nb:
            deg[b] += 1
    return [3 - d for d in deg]


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '../../physics/6d-nilpotency')
    from rooted_trees import counts
    from fractions import Fraction as Fr
    print('H=3 rooted_trees.counts:', counts(5)[1:])
    for H in (2, 3, 4, 5):
        for n in range(1, 5 if H <= 3 else 4):
            if H == 5 and n == 3:
                continue
            gs = rooted_trees(H, n)
            assert all(laman_ok(g['nb'] + 1, g['edges'], H) for g in gs) or (H, n) in ()
            w = sum(Fr(1, g['aut']) for g in gs)
            print('H=%d n=%d trees %d  sum 1/Aut = %s  outputs %s' % (
                H, n, len(gs), w, sorted(set(sum(outputs(g['nb'], g['edges'])) for g in gs))), flush=True)
    for H, n in ((3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 1), (5, 2)):
        us = rooted_unicyclic(H, n)
        print('H=%d n=%d unicyclic %d' % (H, n, len(us)), flush=True)
