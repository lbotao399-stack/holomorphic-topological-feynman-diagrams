"""Exact particular currents from one spanning tree, without tree enumeration.

For a rooted tree, the current from a vertex towards its parent is the sum of
the sources in that vertex's subtree.  Each subtree is visited once.  This
constructs B^t y0 = lambda; the loop solve supplies the electrical correction.
The sources may occupy every vertex and have any number of components.
"""
from fractions import Fraction

from .brute import sources

__all__ = ["tree_flow"]


def tree_flow(G, lam, tree=None):
    """Return exact particular currents, as an ``ne x H`` list of Fractions.

    ``lam`` accepts the same rational source formats as ``htk.brute.sources``.
    The root is ``G.root``.  With no ``tree``, use the graph's deterministic
    Kruskal tree.  A supplied tree is a sequence of edge indices; its validity
    is checked by traversal.  Parallel and reversed edges keep their original
    indices and orientations.  Zero sources give the zero flow.

    Once a spanning tree is chosen, the work is O(ne + nv H), and no spanning
    tree polynomial, matrix inverse, or determinant is constructed.
    """
    if G.nv < 1:
        raise ValueError("a rooted graph must have at least one vertex")
    lam_nr = sources(G, lam)
    H = len(lam_nr[0]) if lam_nr else 0
    if any(len(row) != H for row in lam_nr):
        raise ValueError("source vectors must have a common dimension")
    tree = G.spanning_tree() if tree is None else list(tree)
    if (len(tree) != G.nv - 1 or len(set(tree)) != len(tree)
            or any(not isinstance(e, int) or not 0 <= e < G.ne for e in tree)):
        raise ValueError("tree must contain nv-1 distinct valid edge indices")
    adjacency = [[] for _ in range(G.nv)]
    for e in tree:
        u, v = G.edges[e]
        adjacency[u].append((v, e))
        adjacency[v].append((u, e))
    order = [G.root]
    parent = {G.root: (None, None)}
    for v in order:
        for u, e in adjacency[v]:
            if u not in parent:
                parent[u] = (v, e)
                order.append(u)
    if len(order) != G.nv:
        raise ValueError("edge set is not a spanning tree")
    total = [list(row) for row in lam_nr] + [[Fraction(0)] * H]
    current = [[Fraction(0)] * H for _ in G.edges]
    for v in reversed(order[1:]):
        u, e = parent[v]
        sign = 1 if G.edges[e][0] == v else -1
        current[e] = [sign * value for value in total[v]]
        for i in range(H):
            total[u][i] += total[v][i]
    return current
