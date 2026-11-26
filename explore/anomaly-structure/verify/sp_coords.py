"""Own series-parallel coordinates for k parallel paths between N = 0 and N' = root (engine labelling
of theta(p,q,r) and banana(k)), and exact Grundmann-Moeller integration of the density
g(xi) = D(t(xi)) |det[t; d_xi t]|  (T = 0: exact if g is a polynomial of degree <= 2s+1 per block).
Coordinates: w in Delta^{k-1}, x^P in Delta^{n_P-1}, t^P_j = x^P_j / w_P.
Free coordinates: w_2..w_k (w_1 = 1 - sum), x^P_2..x^P_{n_P} (x^P_1 = 1 - sum).
"""
from fractions import Fraction as Fr
import itertools
import vk


def paths_graph(ns):
    """vertices: N = 0, interior along path 1, path 2, ..., N' = last. edges along paths N -> N'."""
    V = sum(ns) - len(ns) + 2
    root = V - 1
    edges, interior, nxt = [], [], 1
    for n_ in ns:
        prev, ints = 0, []
        for j in range(n_):
            if j == n_ - 1:
                cur = root
            else:
                cur = nxt
                nxt += 1
                ints.append(cur)
            edges.append((prev, cur))
            prev = cur
        interior.append(ints)
    return V, edges, interior


def tmap(ns, w, xs):
    """w: list of k barycentric coords, xs[P]: list of n_P barycentric coords."""
    return [xs[P][j] / w[P] for P in range(len(ns)) for j in range(ns[P])]


def jac(ns, w, xs):
    """N x N matrix [t; d t/d xi_1; ...], xi = (w_2..w_k, x^1_2.., x^2_2.., ...)."""
    k = len(ns)
    t = tmap(ns, w, xs)
    rows = [t]
    # d/d w_Q, Q = 1..k-1 (0-based), w_0 = 1 - sum
    for Q in range(1, k):
        row = []
        for P in range(k):
            dP = (1 if P == Q else 0) - (1 if P == 0 else 0)
            for j in range(ns[P]):
                row.append(-xs[P][j] / w[P] ** 2 * dP)
        rows.append(row)
    for P0 in range(k):
        for kk in range(1, ns[P0]):
            row = []
            for P in range(k):
                for j in range(ns[P]):
                    if P != P0:
                        row.append(Fr(0))
                    else:
                        row.append(Fr((1 if j == kk else 0) - (1 if j == 0 else 0)) / w[P])
            rows.append(row)
    return t, rows


def g_value(ns, w, xs, lam, H, T, C):
    V, edges, interior = paths_graph(ns)
    t, J = jac(ns, w, xs)
    dJ = vk.det(J)
    q, U = vk.density(V, edges, t, lam, H, T, C)
    return q, U, dJ


def integrate_paths(nv, edges, paths, lam, H, s, C=None):
    """T = 0 kernel of an arbitrarily labelled/oriented graph made of k parallel paths;
    paths[P] = list of edge ids in order from one trivalent vertex to the other.
    t_e = x^P_j / w_P for e = paths[P][j];  g = D(t) |det[t; d_xi t]| (columns permuted: |det| unchanged)."""
    if C is None:
        C = vk.positive_cycles(nv, edges, vk.dfs_tree(nv, edges))
    ns = [len(P) for P in paths]
    k = len(ns)
    order = [e for P in paths for e in P]
    blocks = [vk.gm_rule(k - 1, s)]
    for n_ in ns:
        blocks.append(vk.gm_rule(n_ - 1, s) if n_ >= 2 else [(Fr(1), (Fr(1),))])
    total = Fr(0)
    for combo in itertools.product(*blocks):
        wt = Fr(1)
        for (ww, _) in combo:
            wt *= ww
        w = list(combo[0][1])
        xs = [list(c[1]) for c in combo[1:]]
        tl, J = jac(ns, w, xs)
        t = [None] * len(edges)
        for pos, e in enumerate(order):
            t[e] = tl[pos]
        q, U = vk.density(nv, edges, t, lam, H, 0, C)
        total += wt * q * abs(vk.det(J))
    return total


def integrate_T0(ns, lam, H, s, C=None):
    V, edges, interior = paths_graph(ns)
    if C is None:
        C = vk.positive_cycles(V, edges, vk.dfs_tree(V, edges))
    k = len(ns)
    blocks = [vk.gm_rule(k - 1, s)] if k >= 2 else [[(Fr(1), (Fr(1),))]]
    for n_ in ns:
        blocks.append(vk.gm_rule(n_ - 1, s) if n_ >= 2 else [(Fr(1), (Fr(1),))])
    total = Fr(0)
    for combo in itertools.product(*blocks):
        wt = Fr(1)
        for (ww, _) in combo:
            wt *= ww
        w = list(combo[0][1])
        xs = [list(c[1]) for c in combo[1:]]
        q, U, dJ = g_value(ns, w, xs, lam, H, 0, C)
        total += wt * q * abs(dJ)
    return total
