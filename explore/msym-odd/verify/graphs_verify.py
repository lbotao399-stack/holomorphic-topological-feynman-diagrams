"""Independent enumeration of the two-loop graphs of Q_2 for D = H + T = 2, 3, 4 (claim 'graphs-Q2').
Two-loop connected graphs with min degree >= 2 are subdivisions of theta, dumbbell and figure-eight
cores (a degree-1 vertex v gives S = Gamma minus v with (D-1)|I_S| > D|V_S| - D - 1, shown below
by the check on every edge subset).  Conditions: |I| = 2D + 1, |V| = 2D, and for every nonempty edge
subset S: (D-1)|I_S| <= D|V_S| - D - 1.  O admissible: Gamma \\ O a forest and deg O = n + p
(S2 eq-theta-rule, n = 2, p = #components of Gamma \\ O).  |Aut| = automorphisms fixing O."""
import itertools
import networkx as nx
from networkx.algorithms.isomorphism import MultiGraphMatcher


def subdivide(core_edges, lens):
    """core edges (u, v) with lengths; returns edge list on integer vertices."""
    nodes = sorted({x for e in core_edges for x in e})
    lab = {u: i for i, u in enumerate(nodes)}
    nxt = len(nodes)
    out = []
    for (u, v), L in zip(core_edges, lens):
        prev = lab[u]
        for j in range(L):
            if j == L - 1:
                cur = lab[v]
            else:
                cur = nxt; nxt += 1
            out.append((prev, cur))
            prev = cur
    return nxt, out


def ok_counting(nv, edges, D):
    if len(edges) != 2 * D + 1 or nv != 2 * D:
        return False
    E = len(edges)
    for r in range(1, E + 1):
        for S in itertools.combinations(range(E), r):
            VS = {x for e in S for x in edges[e]}
            if (D - 1) * len(S) > D * len(VS) - D - 1:
                return False
    return True


def canon(nv, edges):
    g = nx.MultiGraph(); g.add_nodes_from(range(nv)); g.add_edges_from(edges)
    return g


for D in (2, 3, 4):
    E = 2 * D + 1
    found = []
    cores = [('theta', [('a', 'b')] * 3, 3), ('dumbbell', [('a', 'a'), ('a', 'b'), ('b', 'b')], 3),
             ('eight', [('a', 'a'), ('a', 'a')], 2)]
    for name, ce, k in cores:
        for lens in itertools.product(range(1, E + 1), repeat=k):
            if sum(lens) != E:
                continue
            nv, ed = subdivide(ce, lens)
            g = canon(nv, ed)
            if any(nx.is_isomorphic(g, h) for _, _, h, _ in found):
                continue
            # also count graphs with a pendant path attached (degree-1 vertices) to confirm exclusion
            found.append((name, lens, g, ok_counting(nv, ed, D)))
    # pendant-vertex graphs: attach one leaf to any vertex of any core subdivision with 2D edges
    pend_ok = 0
    for name, ce, k in cores:
        for lens in itertools.product(range(1, E), repeat=k):
            if sum(lens) != E - 1:
                continue
            nv, ed = subdivide(ce, lens)
            for v in range(nv):
                if ok_counting(nv + 1, ed + [(v, nv)], D):
                    pend_ok += 1
    tight = [(nm, ln, g) for nm, ln, g, ok in found if ok]
    print("D=%d: candidate two-loop cores x subdivisions with |I| = %d: %d; satisfying the counting: %s; with a pendant vertex: %d"
          % (D, E, len(found), [(nm, ln) for nm, ln, _ in tight], pend_ok))
    for nm, ln, g in tight:
        res = []
        for O in g.nodes():
            h = g.copy(); h.remove_node(O)
            p = nx.number_connected_components(h)
            forest = h.number_of_edges() == h.number_of_nodes() - p
            if forest and g.degree(O) == 2 + p:
                # automorphisms fixing O
                nx.set_node_attributes(g, {u: (u == O) for u in g.nodes()}, 'O')
                aut = sum(1 for _ in MultiGraphMatcher(g, g, node_match=lambda a, b: a['O'] == b['O']).isomorphisms_iter())
                nverts = g.number_of_nodes() - 1
                outputs = 3 * nverts - 2 * (g.number_of_edges() - g.degree(O)) - g.degree(O)
                res.append((O, g.degree(O), p, aut, outputs))
        degs = sorted(d for _, d in g.degree())
        print("   %s%s: degrees %s; admissible O (vertex, deg O, p, |Aut fixing O|, outputs): %s" % (nm, ln, degs, res))
