"""Theta colour networks.

Standard labelling of theta(p,q,r) (htk.graph.theta): N = 0, N' = root = nv-1, paths A, B, C with interior
vertices numbered along A, B, C; all edges oriented N -> N'; local order L' (externals first, then incident
half-edges by (edge, end)):
    N : f(A_1, B_1, C_1),  N' : f(A_p, B_q, C_r),  interior v on path P : f(X_v, toward N, toward N').
Contraction (kappa^{-1} on every edge):
    C_std(X) = sum_{a,b,c} f(R_A T^a, R_B T^b, R_C T^c) f(T_a, T_b, T_c),   R_P = ad_{X_{P,1}} ... ad_{X_{P,m}}
(X_{P,1} adjacent to N); the sign (-1)^{#interior} = +1 of the path transport is included (8 interior vertices).

gl_N, kappa = tr: ribbon expansion of C_std in the trace basis T_pi (pi = successor map on the external
letters), coefficients polynomial in N; on traceless inputs this is the sl_N contraction (kappa = tr).
"""
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict
from common import parity, perm_from_lists


def theta_paths(shape):
    """interior vertex lists of A, B, C in the standard labelling."""
    p, q, r = shape
    nxt, out = 1, []
    for n_ in (p, q, r):
        out.append(list(range(nxt, nxt + n_ - 1)))
        nxt += n_ - 1
    return out


def theta_edges(shape):
    p, q, r = shape
    V = p + q + r - 1
    edges, nxt = [], 1
    for n_ in (p, q, r):
        prev = 0
        for j in range(n_):
            cur = V - 1 if j == n_ - 1 else nxt
            if j < n_ - 1:
                nxt += 1
            edges.append((prev, cur))
            prev = cur
    return V, edges


# ------------------------------------------------------------------ Wick sign (flagship Prop. relative Wick sign)
def wick_sign(nv, edges, ext_vertices, local_orders):
    """sigma = sgn(source -> target); source = local orders, vertices in vertex order;
    target = (e,0),(e,1) in edge order, then ('x',v) for external v in vertex order."""
    source = []
    for v in range(nv):
        source.extend(local_orders[v])
    target = []
    for e in range(len(edges)):
        target += [(e, 0), (e, 1)]
    target += [('x', v) for v in sorted(ext_vertices)]
    assert sorted(map(str, source)) == sorted(map(str, target))
    return perm_from_lists(source, target)


def local_orders_Lprime(nv, edges, ext_vertices):
    lo = {}
    for v in range(nv):
        hs = [(e, end) for e, (a, b) in enumerate(edges) for end, w in ((0, a), (1, b)) if w == v]
        hs.sort()
        lo[v] = ([('x', v)] if v in ext_vertices else []) + hs
    return lo


def standard_sigma(shape):
    nv, edges = theta_edges(shape)
    ext = [v for v in range(1, nv - 1)]
    return wick_sign(nv, edges, ext, local_orders_Lprime(nv, edges, ext))


# ------------------------------------------------------------------ ribbon expansion, gl_N
def _canon_word(w):
    if not w:
        return w
    k = min(range(len(w)), key=lambda i: w[i:] + w[:i])
    return w[k:] + w[:k]


def _fierz_gl(traces, g):
    """contract the pair of generator letters ('g',g,0), ('g',g,1) (Casimir of gl_N, kappa = tr).
    traces: list of tuples (cyclic words). Returns list of (traces, Npower)."""
    a, b = ('g', g, 0), ('g', g, 1)
    ia = next(i for i, w in enumerate(traces) if a in w)
    ib = next(i for i, w in enumerate(traces) if b in w)
    if ia != ib:
        wa, wb = traces[ia], traces[ib]
        ka, kb = wa.index(a), wb.index(b)
        # tr(U T^e) tr(T_e V) = tr(U V): U = word after a (cyclically), V = word after b
        U = wa[ka + 1:] + wa[:ka]
        V = wb[kb + 1:] + wb[:kb]
        rest = [w for i, w in enumerate(traces) if i not in (ia, ib)]
        return [(rest + [U + V], 0)]
    w = traces[ia]
    ka = w.index(a)
    w = w[ka:] + w[:ka]
    kb = w.index(b)
    U = w[1:kb]
    V = w[kb + 1:]
    rest = [x for i, x in enumerate(traces) if i != ia]
    # tr(T^e U T_e V) = tr(U) tr(V)
    return [(rest + [U, V], 0)]


def _ad_expand(Xs, gen):
    """R(Y) = ad_{X_1}...ad_{X_m}(Y) = sum_S (-1)^{m-|S|} X_{S up} Y X_{S^c down}; returns [(word, sign)]."""
    m = len(Xs)
    out = []
    for k in range(m + 1):
        for S in combinations(range(m), k):
            left = tuple(Xs[i] for i in S)
            right = tuple(Xs[i] for i in reversed(range(m)) if i not in S)
            out.append((left + (gen,) + right, (-1) ** (m - k)))
    return out


def ribbon_theta_gl(shape):
    """C_std for gl_N in the trace basis: dict {canonical trace structure: {Npower: coefficient}}.
    External letters are the vertex labels 1..8."""
    A, B, C = theta_paths(shape)
    gens = [('g', i, 0) for i in range(3)]            # T^a, T^b, T^c in the first f
    gens2 = [('g', i, 1) for i in range(3)]           # T_a, T_b, T_c in the second f
    exA, exB, exC = _ad_expand(A, gens[0]), _ad_expand(B, gens[1]), _ad_expand(C, gens[2])
    out = defaultdict(lambda: defaultdict(int))
    # f(U,V,W) = tr(UVW) - tr(VUW)
    for wa, sa in exA:
        for wb, sb in exB:
            for wc, sc in exC:
                for t1, s1 in (((wa, wb, wc), 1), ((wb, wa, wc), -1)):
                    for t2, s2 in (((gens2[0], gens2[1], gens2[2]), 1), ((gens2[1], gens2[0], gens2[2]), -1)):
                        traces = [t1[0] + t1[1] + t1[2], tuple(t2)]
                        stack = [(traces, 0)]
                        for g in range(3):
                            new = []
                            for tr, npw in stack:
                                for tr2, d in _fierz_gl(tr, g):
                                    new.append((tr2, npw + d))
                            stack = new
                        for tr, npw in stack:
                            empties = sum(1 for w in tr if len(w) == 0)
                            key = tuple(sorted((_canon_word(w) for w in tr if w), key=lambda w: (len(w), w)))
                            out[key][npw + empties] += sa * sb * sc * s1 * s2
    return {k: {n: c for n, c in v.items() if c} for k, v in out.items() if any(v.values())}


def traces_to_perm(key, letters):
    """trace structure (tuple of cyclic words on `letters`) -> successor tuple on positions 0..k-1."""
    pos = {x: i for i, x in enumerate(letters)}
    succ = [None] * len(letters)
    for w in key:
        for i, x in enumerate(w):
            succ[pos[x]] = pos[w[(i + 1) % len(w)]]
    return tuple(succ)


def colour_perm_dict(shape, N, traceless=True):
    """C_std as {perm (successor tuple on 0..k-1, letter i = vertex i+1): coefficient at numeric N}."""
    rib = ribbon_theta_gl(shape)
    k = sum(shape) - 3
    out = {}
    for key, poly in rib.items():
        if traceless and any(len(w) == 1 for w in key):
            continue
        val = sum(Fr(c) * Fr(N) ** n for n, c in poly.items())
        if val:
            perm = traces_to_perm(key, list(range(1, k + 1)))
            out[perm] = out.get(perm, 0) + val
    return {k: v for k, v in out.items() if v}


# ------------------------------------------------------------------ explicit evaluation (any Lie algebra, dense numpy)
def eval_theta_dense(shape, f_low, kinv, ad_of, Xs, mod=None):
    """C_std(X) for a Lie algebra given by dense f_low[a,b,c] = kappa([e_a,e_b],e_c), kinv = kappa^{-1},
    ad_of(X) = matrix of ad_X (column b -> [X, e_b]).  Xs: dict vertex -> coordinate vector.
    With mod: integer arithmetic mod p (kinv given as integer matrix s*kappa^{-1}, result times s^3)."""
    import numpy as np
    A, B, C = theta_paths(shape)

    def red(M):
        return M % mod if mod else M

    def R(path):
        d = f_low.shape[0]
        M = np.eye(d, dtype=f_low.dtype)
        for v in path:
            M = red(M @ ad_of(Xs[v]))
        return M
    MA, MB, MC = (red(R(P) @ kinv) for P in (A, B, C))   # M^{x a} = (R kappa^{-1})^{x a}
    T = red(np.einsum('xyz,xa->ayz', f_low, MA))
    T = red(np.einsum('ayz,yb->abz', T, MB))
    T = red(np.einsum('abz,zc->abc', T, MC))
    val = (T * f_low).sum()
    return red(val) if mod else val
