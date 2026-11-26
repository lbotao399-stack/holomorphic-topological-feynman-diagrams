"""Graphs with a nonzero holomorphic kernel on C^d, and the reach of the reduction moves.

Counting (LLC main.tex Lemmas "Laman count", "Laman bound"; 6d note eq. (5.11)):
    (d-1)|I| = d|V| - d - 1,              i.e. |I| = d n + 1, |V| = (d-1) n + 2,
    (d-1)|I_S| <= d|V_S| - d - 1          for every induced S with an edge.

Classes reported for each (d, n):
    total   all connected graphs satisfying the two conditions
    SP      series-parallel (no K4 minor): reducible by series eliminations and parallel merges
    CH      reducible to an edge by replacing k >= 2 parallel chains with sum a_i = d(k-1)+1 by an edge
    YD      d = 2 only: reducible to an edge by LLC moves (a) series-parallel and
            (b) star-triangle with two existing edges, staying inside Laman graphs
    cubic   maximal degree <= 3

d = 2 totals are checked against the known numbers of Laman graphs 1, 1, 1, 3, 13, 70, 608.
"""
import itertools
import sys
import warnings

import networkx as nx

warnings.filterwarnings("ignore", category=UserWarning, module="networkx")
from networkx.generators.atlas import graph_atlas_g

LAMAN_COUNTS = {3: 1, 4: 1, 5: 3, 6: 13, 7: 70, 8: 608}


def sparse_ok(G, d):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * len(nodes)
    for a, b in G.edges():
        adj[idx[a]] |= 1 << idx[b]
        adj[idx[b]] |= 1 << idx[a]
    for S in range(3, 1 << len(nodes)):
        k = bin(S).count("1")
        if k < 2:
            continue
        e2, T, i = 0, S, 0
        while T:
            if T & 1:
                e2 += bin(adj[i] & S).count("1")
            T >>= 1
            i += 1
        e = e2 // 2
        if e and (d - 1) * e > d * k - d - 1:
            return False
    return True


def tight_ok(G, d):
    V, E = G.number_of_nodes(), G.number_of_edges()
    return (d - 1) * E == d * V - d - 1 and nx.is_connected(G) and sparse_ok(G, d)


def dedupe(graphs):
    buckets = {}
    for G in graphs:
        key = (tuple(sorted(dict(G.degree()).values())), nx.weisfeiler_lehman_graph_hash(G, iterations=4))
        lst = buckets.setdefault(key, [])
        if not any(nx.is_isomorphic(G, H) for H in lst):
            lst.append(G)
    return [G for lst in buckets.values() for G in lst]


ATLAS = graph_atlas_g()


def atlas_sparse(V, E, d):
    return [G for G in ATLAS if G.number_of_nodes() == V and G.number_of_edges() == E and sparse_ok(G, d)]


def enumerate_tight(d, n):
    V, E = (d - 1) * n + 2, d * n + 1
    if V <= 7:
        return [G for G in atlas_sparse(V, E, d) if nx.is_connected(G)]
    if V != 8:
        raise ValueError("only |V| <= 8 is implemented")
    # every tight graph on 8 vertices has a vertex of degree 2 (d >= 2) or 3 (d = 2); remove it
    out = []
    for k in ((2, 3) if d == 2 else (2,)):
        for H in atlas_sparse(7, E - k, d):
            for nb in itertools.combinations(H.nodes(), k):
                G = H.copy()
                G.add_node(7)
                G.add_edges_from((7, v) for v in nb)
                if tight_ok(G, d):
                    out.append(G)
    return dedupe(out)


def is_sp(G):
    """No K4 minor: delete vertices of degree <= 2, joining the neighbours of a degree-2 vertex."""
    H = G.copy()
    while H.number_of_nodes() > 2:
        v = next((u for u in H.nodes() if H.degree(u) <= 2), None)
        if v is None:
            return False
        nb = list(H.neighbors(v))
        H.remove_node(v)
        if len(nb) == 2:
            H.add_edge(*nb)
    return True


def yd_reducible(G):
    """d = 2: LLC move (a) (degree 2, neighbours adjacent) and move (b) (degree 3, two of the
    three neighbour pairs adjacent), each step required to give a Laman graph."""
    if G.number_of_nodes() == 2:
        return True
    for v in list(G.nodes()):
        nb = list(G.neighbors(v))
        if len(nb) == 2 and G.has_edge(*nb):
            H = G.copy()
            H.remove_node(v)
            if tight_ok(H, 2) and yd_reducible(H):
                return True
        if len(nb) == 3:
            pairs = list(itertools.combinations(nb, 2))
            present = [p for p in pairs if G.has_edge(*p)]
            if len(present) == 2:
                H = G.copy()
                H.remove_node(v)
                H.add_edge(*next(p for p in pairs if not G.has_edge(*p)))
                if tight_ok(H, 2) and yd_reducible(H):
                    return True
    return False


def chains_between(G, p, q):
    """Chains from p to q whose interior vertices all have degree two: list of interiors."""
    out = []
    for w in G.neighbors(p):
        prev, cur, interior = p, w, []
        while True:
            if cur == q:
                out.append(tuple(interior))
                break
            if cur == p or G.degree(cur) != 2:
                break
            interior.append(cur)
            prev, cur = cur, next(u for u in G.neighbors(cur) if u != prev)
    return out


def chain_reducible(G, d):
    """Reducible to one edge by moves replacing k >= 2 parallel chains of lengths a_i with
    sum a_i = d (k - 1) + 1 by a single edge (E1 of PLAN.md)."""
    if G.number_of_nodes() == 2:
        return True
    nodes = list(G.nodes())
    for p, q in itertools.combinations(nodes, 2):
        ch = chains_between(G, p, q)
        for k in range(2, len(ch) + 1):
            for sub in itertools.combinations(ch, k):
                if sum(len(x) + 1 for x in sub) != d * (k - 1) + 1:
                    continue
                H = G.copy()
                H.remove_nodes_from(v for x in sub for v in x)
                H.add_edge(p, q)
                if chain_reducible(H, d):
                    return True
    return False


def skeleton(G):
    """Suppress degree-2 vertices; returns the multigraph on vertices of degree >= 3 as a
    sorted list of (multiplicity-labelled) chain lengths between branch vertices."""
    branch = [v for v in G.nodes() if G.degree(v) != 2]
    if not branch:
        return "cycle", [G.number_of_nodes()]
    chains = []
    seen = set()
    for p in branch:
        for w in G.neighbors(p):
            prev, cur, length = p, w, 1
            interior = []
            while G.degree(cur) == 2:
                interior.append(cur)
                nxt = next(u for u in G.neighbors(cur) if u != prev)
                prev, cur, length = cur, nxt, length + 1
            key = frozenset(interior) if interior else frozenset((p, cur))
            if key in seen:
                continue
            seen.add(key)
            chains.append((min(p, cur), max(p, cur), length))
    S = nx.MultiGraph()
    S.add_nodes_from(branch)
    for a, b, L in chains:
        S.add_edge(a, b, length=L)
    return S, sorted(L for _, _, L in chains)


def describe(G):
    S, lengths = skeleton(G)
    if S == "cycle":
        return f"cycle C{lengths[0]}"
    und = nx.Graph(S)
    name = f"skeleton |V|={S.number_of_nodes()} |E|={S.number_of_edges()}"
    if S.number_of_nodes() == 4 and und.number_of_edges() == 6 and S.number_of_edges() == 6:
        name = "K4"
    elif S.number_of_nodes() == 6 and S.number_of_edges() == 9 and nx.is_isomorphic(und, nx.complete_bipartite_graph(3, 3)):
        name = "K33"
    elif S.number_of_nodes() == 6 and S.number_of_edges() == 9 and nx.is_isomorphic(und, nx.circular_ladder_graph(3)):
        name = "prism"
    degs = sorted((d for _, d in G.degree()), reverse=True)
    return f"{name}, chain lengths {lengths}, degrees {degs}"


def report(d, n, show=6):
    Gs = enumerate_tight(d, n)
    V, E = (d - 1) * n + 2, d * n + 1
    sp = [G for G in Gs if is_sp(G)]
    non_sp = [G for G in Gs if not is_sp(G)]
    cubic = [G for G in Gs if max(dict(G.degree()).values()) <= 3]
    cubic_non_sp = [G for G in non_sp if max(dict(G.degree()).values()) <= 3]
    ch = [G for G in Gs if chain_reducible(G, d)]
    line = (f"d={d} n={n} |V|={V} |I|={E}: total={len(Gs)} SP={len(sp)} nonSP={len(non_sp)} CH={len(ch)}"
            f" cubic={len(cubic)} cubic_nonSP={len(cubic_non_sp)}")
    assert all(is_sp(G) for G in ch)
    if d == 2:
        yd = [G for G in Gs if yd_reducible(G)]
        line += f" YD={len(yd)}"
        if V in LAMAN_COUNTS:
            line += f" [Laman count check: {len(Gs) == LAMAN_COUNTS[V]}]"
        rest = [G for G in Gs if not yd_reducible(G)]
    else:
        rest = non_sp
    print(line, flush=True)
    for G in rest[:show]:
        print("    not reduced:", describe(G), "| edges", sorted(tuple(sorted(e)) for e in G.edges()))
    if len(rest) > show:
        print(f"    ... {len(rest) - show} more")
    return Gs


def two_loop_thetas(d):
    """n = 2 for any d: sum(deg - 2) = 2 leaves theta graphs with path lengths a <= b <= c,
    a + b + c = 2d + 1 and every cycle of length >= d + 1, i.e. c <= d.  (A vertex of degree 4
    would need two cycles of total length 2d + 1 < 2(d + 1).)"""
    return [(a, b, 2 * d + 1 - a - b) for a in range(1, d + 1) for b in range(a, d + 1)
            if b <= 2 * d + 1 - a - b <= d]


def cone_over_forest(G, n, max_deg=3):
    """Selection rule of holomorphic Chern-Simons on C^{3|2} (6d note eq. (5.13)): a vertex O with
    Gamma - O a forest of p trees and deg O = n + p; all other vertices of degree <= max_deg."""
    out = []
    for O in G.nodes():
        H = G.copy()
        H.remove_node(O)
        if not nx.is_forest(H):
            continue
        p = nx.number_connected_components(H)
        if G.degree(O) == n + p and all(G.degree(v) <= max_deg for v in H.nodes()):
            out.append((O, p))
    return out


def report_cone(n):
    Gs = enumerate_tight(3, n)
    hits = []
    for G in Gs:
        picks = cone_over_forest(G, n)
        if picks:
            # inequivalent choices of O: group by (degree of O, sorted tree sizes)
            kinds = set()
            for O, p in picks:
                H = G.copy()
                H.remove_node(O)
                kinds.add((G.degree(O), p, tuple(sorted(len(c) for c in nx.connected_components(H)))))
            hits.append((G, kinds))
    print(f"d=3 n={n}: graphs of Q_n in hCS on C^(3|2) (cone over a forest, deg O = n+p): {len(hits)}"
          f" of {len(Gs)}; all SP: {all(is_sp(G) for G, _ in hits)};"
          f" all CH: {all(chain_reducible(G, 3) for G, _ in hits)}", flush=True)
    for G, kinds in hits:
        print("    ", describe(G), "| (deg O, p, tree sizes):", sorted(kinds))


if __name__ == "__main__":
    jobs = [(2, n) for n in range(1, 7)] + [(3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 1), (6, 1)]
    if len(sys.argv) > 1:
        jobs = [tuple(int(x) for x in a.split(",")) for a in sys.argv[1:]]
    for d, n in jobs:
        report(d, n)
    if len(sys.argv) == 1:
        for d in range(2, 7):
            print(f"d={d} n=2 theta graphs (a,b,c): {two_loop_thetas(d)}")
        for n in (1, 2, 3):
            report_cone(n)
