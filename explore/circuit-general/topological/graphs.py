"""graphs.py -- own enumeration of admissible (D-Laman, D = H+T) graphs and named graphs.

Admissible: connected, |I| = D n + 1, |V| = (D-1) n + 2, and for every vertex set S, |S| >= 2,
    f(S) = D |V_S| - D - 1 - (D-1) |I_S| >= 0      (Gaiotto-Kulp-Wu, Higher.tex l.633-655, n -> D).
D = 1: the banana with n+1 parallel edges (two vertices).
D >= 2: no multiple edges, no bridges, no cut vertices; Gamma = subdivision of a 2-connected loopless
multigraph skeleton with minimum degree >= 3 (n >= 2) or of a single cycle (n = 1).
Reduction of the local check: adding a vertex with d edges into S changes f by D - (D-1) d; a chain of
length l with both ends in S changes f by l - D when completed; partial chains raise f.  Hence
min f = min over skeleton sets Q of  D|Q| - D - 1 + sum_{chains in Q} (l == 1 ? 1-D : min(0, l-D)).
"""
from itertools import combinations, combinations_with_replacement
import networkx as nx
from tb import Graph


def _multigraph_ok(s, edges):
    deg = [0] * s
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    if min(deg) < 3:
        return False
    Gs = nx.Graph()
    Gs.add_nodes_from(range(s))
    Gs.add_edges_from(edges)
    if not nx.is_connected(Gs):
        return False
    if s >= 3 and not nx.is_biconnected(Gs):
        return False
    return True


def skeletons(n):
    """2-connected loopless multigraphs, min degree >= 3, cyclomatic number n >= 2, up to isomorphism"""
    out = []
    for s in range(2, 2 * (n - 1) + 1):
        m = s + n - 1
        pairs = list(combinations(range(s), 2))
        found = []
        for ms in combinations_with_replacement(pairs, m):
            if not _multigraph_ok(s, ms):
                continue
            MG = nx.MultiGraph()
            MG.add_nodes_from(range(s))
            MG.add_edges_from(ms)
            if any(nx.is_isomorphic(MG, H) for H in found):
                continue
            found.append(MG)
            out.append((s, list(ms)))
    return out


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def _local_ok(D, s, sk_edges, lens):
    """min over S of f(S): skeleton part Q (|Q| >= 2); a chain of length l with both ends in Q changes f by
    l - D when all its interior vertices are added (forced for l = 1, optional for l >= 2); partial
    chains only raise f (each added vertex of induced degree 1 adds +1)."""
    for q in range(2, s + 1):
        for Q in combinations(range(s), q):
            Qs = set(Q)
            f = D * q - D - 1
            for (a, b), l in zip(sk_edges, lens):
                if a in Qs and b in Qs:
                    f += (1 - D) if l == 1 else min(0, l - D)
            if f < 0:
                return False
    return True


def subdivide(s, sk_edges, lens):
    """simple graph from skeleton with chain lengths; returns (nv, edge list)"""
    nv = s
    E = []
    for (a, b), l in zip(sk_edges, lens):
        prev = a
        for j in range(l):
            if j == l - 1:
                nxt = b
            else:
                nxt = nv
                nv += 1
            E.append((prev, nxt))
            prev = nxt
    return nv, E


def admissible(D, n):
    """list of Graph objects (vertices relabelled so that the edge list is sorted, root = last label)"""
    N = D * n + 1
    if D == 1:
        return [Graph(2, [(0, 1)] * (n + 1), name="banana%d" % (n + 1))]
    if n == 1:
        return [Graph(N, [(i, (i + 1) % N) for i in range(N)], name="C%d" % N)]
    cands = []
    for s, sk in skeletons(n):
        m = len(sk)
        for lens in compositions(N, m):
            if not _local_ok(D, s, sk, lens):
                continue
            nv, E = subdivide(s, sk, lens)
            if len(set(tuple(sorted(e)) for e in E)) < len(E):
                continue   # multiple edge
            cands.append((nv, E, s, sk, lens))
    # dedupe
    buckets = {}
    out = []
    for nv, E, s, sk, lens in cands:
        Gx = nx.Graph()
        Gx.add_nodes_from(range(nv))
        Gx.add_edges_from(E)
        h = nx.weisfeiler_lehman_graph_hash(Gx)
        lst = buckets.setdefault(h, [])
        if any(nx.is_isomorphic(Gx, H) for H in lst):
            continue
        lst.append(Gx)
        # relabel: sorted nodes, edges from smaller to larger label, sorted (as htk catalog.from_networkx)
        edges = sorted(tuple(sorted(e)) for e in Gx.edges())
        out.append(Graph(nv, edges, name="D%dn%d_sk%d_%s" % (D, n, s, "".join(map(str, lens)))))
    return out


# ---------------------------------------------------------------- named graphs

def banana(k):
    return Graph(2, [(0, 1)] * k, name="banana%d" % k)


def theta(p, q, r):
    """two trivalent vertices N=0 and N'=nv-1 (root), three paths of p, q, r edges oriented N -> N'"""
    nv = 2 + (p - 1) + (q - 1) + (r - 1)
    root = nv - 1
    edges = []
    v = 1
    for L in (p, q, r):
        prev = 0
        for j in range(L):
            if j == L - 1:
                nxt = root
            else:
                nxt = v
                v += 1
            edges.append((prev, nxt))
            prev = nxt
    return Graph(nv, edges, name="theta%d%d%d" % (p, q, r))


def bitriangle_gkw():
    """GKW Higher.tex l.1475-1500 (fig:bitriangle), edges t12, t23, t24, t31, t43 (the order of their
    dVol_t, l.1524-1532); GKW vertices v1..v4 -> 0..3, root = v4."""
    return Graph(4, [(0, 1), (1, 2), (1, 3), (2, 0), (3, 2)], name="bitriangle")
