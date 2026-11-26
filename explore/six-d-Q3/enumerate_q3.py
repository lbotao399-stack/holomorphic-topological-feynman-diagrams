"""Graphs of Q_3 of holomorphic Chern-Simons theory on C^{3|2} (S2 Section 5.4).

Rule (S2 eq-theta-rule, (5.13)):  deg_Gamma O = n + p,  Gamma - O a forest of p trees,
counting (S2 eq-count, (5.11)):  |I| = 3n+1, |V| = 2n+2,  2|I_S| <= 3|V_S| - 4  (I_S nonempty),
cubic interaction: every vertex other than O has degree <= 3 in Gamma.

Enumeration (independent of htk.catalog): all forests F on 2n+1 = 7 vertices with p components,
all multisets of spokes O - v (deg O = n + p) with deg_F(v) + spokes(v) <= 3, then the subgraph
bound on all vertex subsets, then isomorphism classes (O fixed).  Cross-check against
htk.catalog.tight_graphs(3, 3) filtered by cone_over_forest.
"""
import itertools
import sys
from collections import Counter

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')

n = 3
NV = 2 * n + 1          # vertices of Gamma - O
O = 'O'


def sparse_ok(G):
    """2|I_S| <= 3|V_S| - 4 for every vertex subset with at least one edge (induced subgraphs
    maximize |I_S| at fixed V_S).  G may be a MultiGraph."""
    nodes = list(G.nodes())
    for k in range(2, len(nodes) + 1):
        for S in itertools.combinations(nodes, k):
            m = G.subgraph(S).number_of_edges()
            if m and 2 * m > 3 * k - 4:
                return False
    return True


def forests(nv, p):
    """Non-isomorphic forests on nv vertices with p components (via partitions into trees)."""
    out = []
    def parts(m, k, maxp):
        if k == 0:
            if m == 0:
                yield []
            return
        for a in range(min(m, maxp), 0, -1):
            for rest in parts(m - a, k - 1, a):
                yield [a] + rest
    for part in parts(nv, p, nv):
        trees_by_size = [list(nx.nonisomorphic_trees(a)) if a > 1 else [nx.empty_graph(1)] for a in part]
        for combo in itertools.product(*trees_by_size):
            F = nx.Graph()
            off = 0
            for T in combo:
                F.add_edges_from((u + off, v + off) for u, v in T.edges())
                F.add_nodes_from(v + off for v in T.nodes())
                off += T.number_of_nodes()
            out.append(F)
    # dedupe
    uniq = []
    for F in out:
        if not any(nx.is_isomorphic(F, G) for G in uniq):
            uniq.append(F)
    return uniq


def classes(graphs):
    reps = []
    for G in graphs:
        found = False
        for R in reps:
            GM = GraphMatcher(R, G, node_match=lambda a, b: a.get('O', False) == b.get('O', False))
            if GM.is_isomorphic():
                found = True
                break
        if not found:
            reps.append(G)
    return reps


results = {}
for p in range(1, NV + 1):
    degO = n + p
    cands = []
    for F in forests(NV, p):
        if max((d for _, d in F.degree()), default=0) > 3:
            continue
        verts = list(F.nodes())
        # spokes: multiset of degO vertices (multi-edges allowed at this stage, removed by the bound)
        for spokes in itertools.combinations_with_replacement(verts, degO):
            cnt = Counter(spokes)
            if any(F.degree(v) + cnt[v] > 3 for v in verts):
                continue
            G = nx.MultiGraph()
            G.add_nodes_from(verts)
            G.add_edges_from(F.edges())
            G.add_node(O, O=True)
            for v in spokes:
                G.add_edge(O, v)
            assert G.number_of_edges() == 3 * n + 1 and G.number_of_nodes() == 2 * n + 2
            if not sparse_ok(G):
                continue
            cands.append(nx.Graph(G))
    for G in cands:
        nx.set_node_attributes(G, {v: (v == O) for v in G.nodes()}, 'O')
    reps = classes(cands)
    results[p] = reps
    print('p = %d, deg O = %d: candidates passing the bound %d, isomorphism classes %d'
          % (p, degO, len(cands), len(reps)))

reps = results[1]


def describe(G):
    H = G.copy()
    H.remove_node(O)
    degs = sorted((d for _, d in G.degree()), reverse=True)
    deg2 = sorted(v for v in H.nodes() if G.degree(v) == 2)
    aut = sum(1 for _ in GraphMatcher(G, G, node_match=lambda a, b: a['O'] == b['O']).isomorphisms_iter())
    return degs, aut


print()
for k, G in enumerate(reps):
    H = G.copy(); H.remove_node(O)
    degs, aut = describe(G)
    tree_deg = sorted(dict(H.degree()).values(), reverse=True)
    spokes = sorted(v for v in G.neighbors(O))
    print('class %d: degrees %s, tree degree sequence %s, tree diameter %d, |Aut| = %d'
          % (k + 1, degs, tree_deg, nx.diameter(H), aut))
    print('    tree edges', sorted(tuple(sorted(e)) for e in H.edges()), ' spokes at', spokes)

# cross-check with htk.catalog
from htk import catalog
tg = catalog.tight_graphs(3, 3)
sel = []
for g in tg:
    Gx = nx.Graph()
    Gx.add_edges_from(g.edges)
    if nx.number_of_edges(Gx) != g.ne:
        continue
    roots = catalog.cone_over_forest(Gx, 3)
    for (Ov, p) in roots:
        H = nx.relabel_nodes(Gx, {Ov: O})
        nx.set_node_attributes(H, {v: (v == O) for v in H.nodes()}, 'O')
        sel.append(H)
selc = classes(sel)
print('\nhtk.catalog: tight graphs (D=3, n=3): %d; rooted at an admissible O: %d; classes: %d'
      % (len(tg), len(sel), len(selc)))
match = all(any(GraphMatcher(R, G, node_match=lambda a, b: a['O'] == b['O']).is_isomorphic() for R in selc) for G in reps)
print('same classes as this enumeration:', match and len(selc) == len(reps))
