"""Graph encoding of htk (see CONVENTIONS.md, Section 1).

Graph(nv, edges): vertices 0..nv-1, root = nv-1 (grounded vertex, position of the insertion),
edges = ordered list of oriented pairs (tail, head) = (e(0), e(1)).
"""
from .linalg import det, perm_sign

__all__ = [
    "Graph", "polygon", "theta", "banana", "fan", "box", "graph_A", "graph_B", "ladder",
    "glued_polygons", "from_networkx",
]


class Graph:
    def __init__(self, nv, edges, name=None):
        self.nv = int(nv)
        self.edges = [(int(a), int(b)) for a, b in edges]
        self.name = name
        for a, b in self.edges:
            if not (0 <= a < self.nv and 0 <= b < self.nv):
                raise ValueError("edge (%d,%d) outside 0..%d" % (a, b, self.nv - 1))
            if a == b:
                raise ValueError("self-loop at %d" % a)

    # ------------------------------------------------------------------ basic data
    @property
    def ne(self):
        return len(self.edges)

    @property
    def root(self):
        return self.nv - 1

    @property
    def loops(self):
        return self.ne - self.nv + 1

    def __repr__(self):
        nm = "" if self.name is None else self.name + ": "
        return "Graph(%snv=%d, edges=%s)" % (nm, self.nv, self.edges)

    def degree(self, v):
        return sum((a == v) + (b == v) for a, b in self.edges)

    def incidence(self):
        """Grounded incidence matrix B (ne x (nv-1)):  B[e][v] = +1 if v = e(0), -1 if v = e(1)."""
        g = self.nv - 1
        B = [[0] * g for _ in self.edges]
        for e, (a, b) in enumerate(self.edges):
            if a < g:
                B[e][a] += 1
            if b < g:
                B[e][b] -= 1
        return B

    def is_connected(self):
        return len(self.spanning_tree(strict=False)) == self.nv - 1

    def tight(self, D):
        """Counting condition (D-1)|I| = D|V| - D - 1  (|I| = D n + 1, |V| = (D-1) n + 2)."""
        return (D - 1) * self.ne == D * self.nv - D - 1

    # ------------------------------------------------------------------ trees and cycles
    def spanning_tree(self, strict=True):
        """Kruskal tree in edge order: edge e is taken iff it joins two components of the
        edges taken before it.  Returns the sorted list of tree edge indices."""
        parent = list(range(self.nv))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        tree = []
        for e, (a, b) in enumerate(self.edges):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
                tree.append(e)
        if strict and len(tree) != self.nv - 1:
            raise ValueError("graph is not connected")
        return tree

    def cotree(self, tree=None):
        tree = self.spanning_tree() if tree is None else list(tree)
        s = set(tree)
        return [e for e in range(self.ne) if e not in s]

    def _check_tree(self, tree):
        tree = sorted(tree)
        if len(tree) != self.nv - 1:
            raise ValueError("a spanning tree has nv-1 edges")
        B = self.incidence()
        if det([B[e] for e in tree]) == 0:
            raise ValueError("edge set %s is not a spanning tree" % tree)
        return tree

    def fundamental_cycles(self, tree=None):
        """Fundamental cycle matrix C_T (n x ne) of the spanning tree T: row r belongs to the
        r-th cotree edge c_r (increasing order), C[r][c_r] = +1, and the cycle runs along c_r and
        back through T; C[r][e] = +1 (-1) if e is traversed along (against) its orientation."""
        tree = self.spanning_tree() if tree is None else self._check_tree(tree)
        adj = {v: [] for v in range(self.nv)}
        for e in tree:
            a, b = self.edges[e]
            adj[a].append((b, e, +1))
            adj[b].append((a, e, -1))
        C = []
        for c in self.cotree(tree):
            a, b = self.edges[c]
            # tree path from b back to a
            prev = {b: None}
            stack = [b]
            while stack:
                v = stack.pop()
                if v == a:
                    break
                for (u, e, s) in adj[v]:
                    if u not in prev:
                        prev[u] = (v, e, s)
                        stack.append(u)
            row = [0] * self.ne
            row[c] = 1
            v = a
            while prev[v] is not None:
                u, e, s = prev[v]
                row[e] += s          # e was traversed from u to v with sign s
                v = u
            C.append(row)
        return C

    def sigma(self, tree=None):
        """sigma_T = eps_T det(B_T) = +-1.  eps_T: sign of the permutation that moves the tree
        edges (increasing order) in front of the cotree edges (increasing order); B_T: rows of the
        grounded incidence matrix on the tree edges.  det[B^t ; C_T] = sigma_T * (number of
        spanning trees)."""
        tree = self.spanning_tree() if tree is None else self._check_tree(tree)
        order = tree + self.cotree(tree)
        inv = [0] * self.ne
        for pos, e in enumerate(order):
            inv[e] = pos
        B = self.incidence()
        d = det([B[e] for e in tree])
        assert d in (1, -1)
        return perm_sign(inv) * int(d)

    def cycle_matrix(self, tree=None, positive=True):
        """Cycle matrix used in the definition of the reduced kernel K.
        positive=True: the fundamental cycle matrix of T with its first row multiplied by
        sigma_T, so that det[B^t ; C] > 0 (positively oriented basis of the cycle lattice).
        positive=False: the plain fundamental cycle matrix C_T."""
        C = self.fundamental_cycles(tree)
        if positive and C and self.sigma(tree) < 0:
            C[0] = [-x for x in C[0]]
        return C

    # ------------------------------------------------------------------ relabellings
    def reverse_edges(self, which):
        which = set(which)
        return Graph(self.nv, [(b, a) if e in which else (a, b) for e, (a, b) in enumerate(self.edges)],
                     self.name)

    def permute_edges(self, order):
        """New graph whose j-th edge is the old edge order[j]."""
        if sorted(order) != list(range(self.ne)):
            raise ValueError("order must be a permutation of the edges")
        return Graph(self.nv, [self.edges[e] for e in order], self.name)

    def relabel(self, new_label):
        """New graph in which old vertex v is called new_label[v]."""
        if sorted(new_label) != list(range(self.nv)):
            raise ValueError("new_label must be a permutation of the vertices")
        return Graph(self.nv, [(new_label[a], new_label[b]) for a, b in self.edges], self.name)

    def to_networkx(self):
        import networkx as nx
        G = nx.MultiGraph()
        G.add_nodes_from(range(self.nv))
        for e, (a, b) in enumerate(self.edges):
            G.add_edge(a, b, key=e)
        return G


def from_networkx(G, root=None, name=None):
    """Graph from a networkx (multi)graph: vertices relabelled 0..nv-1 in sorted order with `root`
    (default: the largest label) last; edges in sorted order, oriented from the smaller new label
    to the larger one.  Returns (Graph, dict old label -> new label)."""
    nodes = sorted(G.nodes())
    if root is None:
        root = nodes[-1]
    nodes = [v for v in nodes if v != root] + [root]
    lab = {v: i for i, v in enumerate(nodes)}
    edges = sorted((min(lab[a], lab[b]), max(lab[a], lab[b])) for a, b in G.edges())
    return Graph(len(nodes), edges, name), lab


# ---------------------------------------------------------------------- named graphs

def polygon(m):
    """m-gon 0 -> 1 -> ... -> m-1 -> 0, edges e_k = (k, k+1 mod m); root = m-1.
    Nonzero for D = m - 1 (S4 check_single_polygon, S3 heptagon_check: m = 7)."""
    return Graph(m, [(k, (k + 1) % m) for k in range(m)], "polygon%d" % m)


def theta(p, q, r):
    """Theta graph of S3 (forms_check.theta): N = 0, N' = root, paths A, B, C from N to N' with
    p, q, r edges, interior vertices numbered along A, then B, then C; edge order A, B, C.
    Returns (Graph, interior) with interior = [list of interior vertices of A, B, C]."""
    V = p + q + r - 1
    edges, interior, nxt = [], [], 1
    for n_ in (p, q, r):
        prev, ints = 0, []
        for j in range(n_):
            if j == n_ - 1:
                cur = V - 1
            else:
                cur = nxt
                nxt += 1
                ints.append(cur)
            edges.append((prev, cur))
            prev = cur
        interior.append(ints)
    return Graph(V, edges, "theta(%d,%d,%d)" % (p, q, r)), interior


def banana(k):
    """Two vertices joined by k edges 0 -> 1."""
    return Graph(2, [(0, 1)] * k, "banana%d" % k)


def fan(n):
    """n-loop fan of S1 (Theorem thm-fan-kernel): path 0 -> 1 -> ... -> n, spokes (i, n+1) towards
    the insertion n+1 = root.  S1 vertex i is vertex i-1 here.  Edge order of S1 (Theorem b-01142
    applied to the eliminations of vertices 1, 2, ...): (i,i+1), (i,root) for i = 0..n-1, then
    (n, root)."""
    root = n + 1
    edges = []
    for i in range(n):
        edges += [(i, i + 1), (i, root)]
    edges.append((n, root))
    return Graph(n + 2, edges, "fan%d" % n)


def box():
    """Box of S2 Section 6: O = 3 (root), v_1, v_2, v_3 = 0, 1, 2;
    edges (O->1), (1->2), (2->3), (3->O)."""
    return Graph(4, [(3, 0), (0, 1), (1, 2), (2, 3)], "box")


def graph_A():
    """Double box (A) of S2 Section 7: vertices 1..5 -> 0..4, O = 5;
    edges (O->1), (1->2), (2->3), (O->3), (3->4), (4->5), (5->O)."""
    return Graph(6, [(5, 0), (0, 1), (1, 2), (5, 2), (2, 3), (3, 4), (4, 5)], "A")


def graph_B():
    """Spider (B) of S2 Section 7: u, w, q, a, b -> 0, 1, 2, 3, 4, O = 5;
    edges (O->u), (u->q), (O->w), (w->q), (q->a), (a->b), (b->O)."""
    return Graph(6, [(5, 0), (0, 2), (5, 1), (1, 2), (2, 3), (3, 4), (4, 5)], "B")


def ladder(n):
    """n-box ladder on C^3 of S4 (check_ladder.ladder): a_i = 2i, b_i = 2i+1, i = 0..n.
    Returns (Graph, cotree of S4, a, b)."""
    a = lambda i: 2 * i
    b = lambda i: 2 * i + 1
    edges, cotree = [], []
    for i in range(n - 1):
        cotree.append(len(edges))
        edges += [(a(i + 1), a(i)), (a(i), b(i)), (b(i), b(i + 1))]
    edges += [(a(n - 1), b(n - 1)), (b(n - 1), b(n)), (b(n), a(n)), (a(n), a(n - 1))]
    cotree.append(len(edges) - 3)
    return Graph(2 * n + 2, edges, "ladder%d" % n), cotree, a, b


def glued_polygons(d):
    """(d+1)-gon chain glued on an edge of a base (d+1)-gon (S4 check_cd.glued_polygons).
    Returns (Graph, cotree of S4, chain, base)."""
    u, w = d - 1, d
    chain = [u] + list(range(d - 1)) + [w]
    base = [u, w] + list(range(d + 1, 2 * d))
    edges = [(chain[k], chain[k + 1]) for k in range(d)]
    edges.append((u, w))
    edges += [(base[k], base[(k + 1) % (d + 1)]) for k in range(1, d + 1)]
    return Graph(2 * d, edges, "glued%d" % d), [0, d + 1], chain, base
