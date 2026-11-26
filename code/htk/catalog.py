"""Layer e): graph utilities, taken from S6 (code/enumerate_graphs.py) with the dimension d
read as D = H + T.

Counting:  (D-1)|I| = D|V| - D - 1  (|I| = D n + 1, |V| = (D-1) n + 2),
           (D-1)|I_S| <= D|V_S| - D - 1 for every induced S with an edge.

Functions on networkx graphs (as in S6): sparse_ok, tight_ok, enumerate_tight, is_sp,
yd_reducible, chains_between, chain_reducible, skeleton, describe, two_loop_thetas,
cone_over_forest.  Functions returning htk.Graph objects: tight_graphs, theta_graphs.
"""
import itertools
import warnings

import networkx as nx

from .graph import Graph, from_networkx, theta

__all__ = ["sparse_ok", "tight_ok", "enumerate_tight", "is_sp", "yd_reducible", "chains_between",
           "chain_reducible", "skeleton", "describe", "two_loop_thetas", "cone_over_forest",
           "tight_graphs", "theta_graphs", "is_series_parallel", "is_chain_reducible", "subgraph_bound_ok",
           "LAMAN_COUNTS"]

LAMAN_COUNTS = {3: 1, 4: 1, 5: 3, 6: 13, 7: 70, 8: 608}
_ATLAS = None


def _atlas():
    global _ATLAS
    if _ATLAS is None:
        warnings.filterwarnings("ignore", category=UserWarning, module="networkx")
        from networkx.generators.atlas import graph_atlas_g
        _ATLAS = graph_atlas_g()
    return _ATLAS


def sparse_ok(G, d):
    """(d-1)|I_S| <= d|V_S| - d - 1 for every induced subgraph S with an edge."""
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


def _dedupe(graphs):
    buckets = {}
    for G in graphs:
        key = (tuple(sorted(dict(G.degree()).values())), nx.weisfeiler_lehman_graph_hash(G, iterations=4))
        lst = buckets.setdefault(key, [])
        if not any(nx.is_isomorphic(G, H) for H in lst):
            lst.append(G)
    return [G for lst in buckets.values() for G in lst]


def _atlas_sparse(V, E, d):
    return [G for G in _atlas() if G.number_of_nodes() == V and G.number_of_edges() == E and sparse_ok(G, d)]


def enumerate_tight(d, n):
    """All connected simple graphs with |V| = (d-1) n + 2, |I| = d n + 1 satisfying the subgraph
    bound, up to isomorphism (networkx graphs).  Implemented for |V| <= 8 (S6)."""
    V, E = (d - 1) * n + 2, d * n + 1
    if V <= 7:
        return [G for G in _atlas_sparse(V, E, d) if nx.is_connected(G)]
    if V != 8:
        raise ValueError("only |V| <= 8 is implemented")
    out = []
    for k in ((2, 3) if d == 2 else (2,)):
        for H in _atlas_sparse(7, E - k, d):
            for nb in itertools.combinations(H.nodes(), k):
                G = H.copy()
                G.add_node(7)
                G.add_edges_from((7, v) for v in nb)
                if tight_ok(G, d):
                    out.append(G)
    return _dedupe(out)


def is_sp(G):
    """No K4 minor: delete vertices of degree <= 2, joining the neighbours of a degree-2 vertex."""
    H = nx.Graph(G)
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
    """d = 2: S1 move (a) (degree 2, neighbours adjacent) and move (b) (degree 3, two of the three
    neighbour pairs adjacent), each step required to give a Laman graph."""
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
    sum a_i = d (k - 1) + 1 by a single edge."""
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
    """Suppress degree-2 vertices: (multigraph on the vertices of degree != 2, sorted chain lengths)."""
    branch = [v for v in G.nodes() if G.degree(v) != 2]
    if not branch:
        return "cycle", [G.number_of_nodes()]
    chains, seen = [], set()
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
        return "cycle C%d" % lengths[0]
    und = nx.Graph(S)
    name = "skeleton |V|=%d |E|=%d" % (S.number_of_nodes(), S.number_of_edges())
    if S.number_of_nodes() == 4 and und.number_of_edges() == 6 and S.number_of_edges() == 6:
        name = "K4"
    elif S.number_of_nodes() == 6 and S.number_of_edges() == 9 and nx.is_isomorphic(und, nx.complete_bipartite_graph(3, 3)):
        name = "K33"
    elif S.number_of_nodes() == 6 and S.number_of_edges() == 9 and nx.is_isomorphic(und, nx.circular_ladder_graph(3)):
        name = "prism"
    degs = sorted((d for _, d in G.degree()), reverse=True)
    return "%s, chain lengths %s, degrees %s" % (name, lengths, degs)


def two_loop_thetas(d):
    """n = 2: theta graphs with path lengths a <= b <= c, a + b + c = 2d + 1, c <= d."""
    return [(a, b, 2 * d + 1 - a - b) for a in range(1, d + 1) for b in range(a, d + 1)
            if b <= 2 * d + 1 - a - b <= d]


def cone_over_forest(G, n, max_deg=3):
    """Selection rule of holomorphic Chern-Simons on C^{3|2} (S2 eq-theta-rule): a vertex O with
    G - O a forest of p trees and deg O = n + p; all other vertices of degree <= max_deg.
    Returns the list of (O, p)."""
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


# ---------------------------------------------------------------------- htk.Graph interface

def tight_graphs(D, n, root=None):
    """The graphs of enumerate_tight(D, n) as htk.Graph objects (from_networkx conventions)."""
    return [from_networkx(G, root=root, name="D%d_n%d_%d" % (D, n, k))[0]
            for k, G in enumerate(enumerate_tight(D, n))]


def theta_graphs(D):
    """[(pqr, Graph, interior)] for the two-loop theta graphs at dimension D."""
    return [(pqr,) + theta(*pqr) for pqr in two_loop_thetas(D)]


def is_series_parallel(G):
    return is_sp(G.to_networkx() if isinstance(G, Graph) else G)


def is_chain_reducible(G, D):
    H = nx.Graph(G.to_networkx()) if isinstance(G, Graph) else G
    return chain_reducible(H, D)


def subgraph_bound_ok(G, D):
    H = nx.Graph(G.to_networkx()) if isinstance(G, Graph) else G
    return sparse_ok(H, D)
