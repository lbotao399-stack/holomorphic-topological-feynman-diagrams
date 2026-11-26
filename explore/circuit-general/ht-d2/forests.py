"""Cubic product forests (marked graphs of A*B, LLC Theorem b-00249) and spiders (LLC Definition
b-00380) as htk graphs.

Tree vertex = (spokes, children); every vertex has |spokes| + |children| = 2, spokes in {A, B},
at most one spoke to each insertion.  Root = unique vertex of full degree two (carries the source).
Graph: interaction vertices 0..n-1 numbered depth-first (trees in canonical order), A = n,
B = n+1 = htk root; edges in lexicographic order, oriented from the smaller label; (A,B) last.
"""
import itertools
import os
import sys
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "code"))

from htk import Graph  # noqa: E402

LEAF = (("A", "B"), ())


def _key(node):
    spokes, ch = node
    return (len(spokes), spokes, tuple(_key(c) for c in ch))


@lru_cache(maxsize=None)
def trees(m):
    """Canonical rooted cubic product trees with m vertices."""
    if m == 1:
        return (LEAF,)
    out = []
    for X in ("A", "B"):
        for c in trees(m - 1):
            out.append(((X,), (c,)))
    for m1 in range(1, (m - 1) // 2 + 1):
        m2 = m - 1 - m1
        T1, T2 = trees(m1), trees(m2)
        for i, c1 in enumerate(T1):
            for j, c2 in enumerate(T2):
                if m1 == m2 and j < i:
                    continue
                pair = tuple(sorted((c1, c2), key=_key))
                out.append(((), pair))
    # dedupe
    seen, res = set(), []
    for t in out:
        k = _key(t)
        if k not in seen:
            seen.add(k)
            res.append(t)
    return tuple(res)


def size(node):
    return 1 + sum(size(c) for c in node[1])


def forests(n):
    """Multisets of trees with total size n (canonical order)."""
    out = []

    def parts(rem, maxkey, acc):
        if rem == 0:
            out.append(tuple(acc))
            return
        for m in range(1, rem + 1):
            for t in trees(m):
                k = (m, _key(t))
                if maxkey is not None and k > maxkey:
                    continue
                parts(rem - m, k, acc + [t])

    parts(n, None, [])
    # canonical: sort descending by key handled by maxkey ordering
    return out


def label(node):
    """LLC notation: leaf = (A,B); (X, child); (child1, child2)."""
    spokes, ch = node
    if not ch:
        return "(A,B)"
    if len(ch) == 1:
        return "(%s,%s)" % (spokes[0], label(ch[0]))
    return "(%s,%s)" % (label(ch[0]), label(ch[1]))


def forest_graph(forest):
    """htk Graph of the marked product graph; returns (Graph, roots) with roots = list of the
    interaction vertices carrying the surviving sources (one per tree)."""
    n = sum(size(t) for t in forest)
    A, B = n, n + 1
    edges, roots = [], []
    counter = [0]

    def walk(node, parent):
        v = counter[0]
        counter[0] += 1
        if parent is not None:
            edges.append((parent, v))
        for s in node[0]:
            edges.append((v, A if s == "A" else B))
        for c in node[1]:
            walk(c, v)
        return v

    for t in forest:
        roots.append(walk(t, None))
    edges.append((A, B))
    edges = sorted((min(a, b), max(a, b)) for a, b in edges)
    return Graph(n + 2, edges, " ".join(label(t) for t in forest)), roots


def spider_graph(arms):
    """Spider of LLC Definition b-00380 with arm lengths arms = (l_1, ..., l_a): arm i has vertices
    v_{i,1..l_i} towards the hub A; every arm vertex and A joined to B.  Vertices: arms in order
    (v_{i,1} first), then A = n, B = n+1 = root.  Edge order of LLC: arm by arm, (v_{i,j}, v_{i,j+1})
    before (v_{i,j}, B), (A, B) last.  Returns (Graph, outer) with outer = [v_{i,1}]."""
    n = sum(arms)
    A, B = n, n + 1
    edges, outer, v = [], [], 0
    for l in arms:
        ids = list(range(v, v + l))
        outer.append(ids[0])
        for j, x in enumerate(ids):
            nxt = ids[j + 1] if j + 1 < l else A
            edges.append((x, nxt))
            edges.append((x, B))
        v += l
    edges.append((A, B))
    return Graph(n + 2, edges, "spider%s" % (tuple(arms),)), outer


if __name__ == "__main__":
    for n in range(1, 7):
        print(n, len(forests(n)), [len(trees(m)) for m in range(1, n + 1)])


def tadpole_graph(c, L):
    """Single-source Q_n graph at (H,T) = (1,1) (REPORT.md Section 3): cycle z = 0, 1, ..., c-1
    (cycle edges (i, i+1 mod c)), tail z - c - c+1 - ... - c+L-1 (q = last tail vertex, q = z if
    L = 0), insertion O = n+1 = root joined to every interaction vertex except z; n = c + L - 1.
    Source at q, -source at O.  Edges in lexicographic order.  Returns (Graph, q)."""
    n = c + L - 1
    O = n + 1
    edges = [(i, (i + 1) % c) for i in range(c)]
    prev = 0
    for j in range(L):
        v = c + j
        edges.append((prev, v))
        prev = v
    q = prev
    for v in range(1, n + 1):
        edges.append((v, O))
    edges = sorted((min(a, b), max(a, b)) for a, b in edges)
    return Graph(n + 2, edges, "tadpole(c=%d,L=%d)" % (c, L)), q
