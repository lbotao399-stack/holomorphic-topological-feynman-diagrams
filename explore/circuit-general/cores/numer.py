"""Exact Schwinger numerators of holomorphic kernels (T = 0).

K_Gamma = int_simplex D(t) prod_ins y_e^i dt = int P(t) / U^q,  q = n + H + m,
P = U^(n+H) D * prod_ins (U y_e^i).  U^(n+H) D is interpolated exactly (Newton
differences on a downward-closed grid; degree <= n-1 in each time, total n^2-1 after
fixing t_0 = 1) from htk.brute.density; U y_e^i = sum over spanning trees of the cotree
monomial times the tree flow.  Conventions of code/htk/CONVENTIONS.md (root = last vertex).
"""
import os
import sys
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../../../code'))
from htk import brute                     # noqa: E402
from htk.graph import Graph               # noqa: E402


def lower_grid(size, remaining, degree):
    if size == 0:
        yield ()
        return
    for i in range(min(degree, remaining) + 1):
        for rest in lower_grid(size - 1, remaining - i, degree):
            yield (i,) + rest


def newton_interpolate(samples, nvar, degree_limit, total_limit):
    """samples: dict exps -> value f(1 + exps) on the downward-closed grid
    {e in [0,degree_limit]^nvar, sum e <= total_limit}; returns monomial coefficients
    of the interpolating polynomial in (t_1..t_nvar)."""
    order = sorted(samples, key=sum)
    vals = dict(samples)
    for axis in range(nvar):
        for ex in reversed(order):
            d = ex[axis]
            if not d:
                continue
            acc = 0
            for j in range(d + 1):
                acc += (-1) ** (d - j) * comb(d, j) * vals[ex[:axis] + (j,) + ex[axis + 1:]]
            vals[ex] = F(acc, factorial(d))
    # now vals[ex] = coefficient of prod_i binom(t_i - 1, ex_i) * ex_i!  i.e. falling factorials
    # (t-1)(t-2)...(t-k) basis with Newton coefficients; convert per axis to monomials in t
    # falling[k] = coefficients of prod_{j=1}^{k} (t - j) in powers of t
    falling = [[F(1)]]
    for k in range(1, degree_limit + 1):
        prev = falling[-1]
        cur = [F(0)] * (k + 1)
        for i, c in enumerate(prev):
            cur[i + 1] += c
            cur[i] -= k * c
        falling.append(cur)
    poly = {ex: v for ex, v in vals.items() if v}
    for axis in range(nvar):
        new = {}
        for ex, c in poly.items():
            k = ex[axis]
            for p, fc in enumerate(falling[k]):
                if fc:
                    key = ex[:axis] + (p,) + ex[axis + 1:]
                    new[key] = new.get(key, 0) + c * fc
        poly = {ex: v for ex, v in new.items() if v}
    return poly


def tree_current_polys(G, lam_full):
    """U y_e^i = sum_T prod_{f not in T} t_f * flow_T(e)^i ;  returns dict e -> list over i of
    dict(exponent tuple -> int)."""
    N, nv = G.ne, G.nv
    H = len(lam_full[0])
    out = {e: [dict() for _ in range(H)] for e in range(N)}
    for tree in combinations(range(N), nv - 1):
        # check spanning tree
        parent = list(range(nv))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        ok = True
        for e in tree:
            a, b = find(G.edges[e][0]), find(G.edges[e][1])
            if a == b:
                ok = False
                break
            parent[a] = b
        if not ok:
            continue
        # flows on the tree: solve B^T y = lam on the tree (unique)
        adj = {v: [] for v in range(nv)}
        for e in tree:
            u, w = G.edges[e]
            adj[u].append((w, e))
            adj[w].append((u, e))
        root = nv - 1
        order, par = [root], {root: (None, None)}
        for v in order:
            for w, e in adj[v]:
                if w not in par:
                    par[w] = (v, e)
                    order.append(w)
        flow = {e: [0] * H for e in tree}
        sub = {v: list(lam_full[v]) for v in range(nv)}
        for v in reversed(order[1:]):
            p, e = par[v]
            # net source in subtree of v flows out of v along e toward p
            sgn = 1 if G.edges[e][0] == v else -1
            for i in range(H):
                flow[e][i] += sgn * sub[v][i]
                sub[p][i] += sub[v][i]
        cot = tuple(1 if f not in tree else 0 for f in range(N))
        for e in tree:
            for i in range(H):
                if flow[e][i]:
                    out[e][i][cot] = out[e][i].get(cot, 0) + flow[e][i]
    return out


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def bare_numerator(G, lam, verbose=False):
    """U^(n+H) D as dict over full exponent tuples (homogeneous of degree n^2 - 1)."""
    n, N = G.loops, G.ne
    lam_full = brute.sources(G, lam)
    H = len(lam_full[0])
    if n * H != N - 1:
        raise ValueError('not a T=0 holomorphic graph: n H != N - 1')
    total = n * n - 1
    deg = n - 1
    samples = {}
    C = G.cycle_matrix()
    cnt = 0
    for ex in lower_grid(N - 1, total, deg):
        t = [1] + [1 + x for x in ex]
        D = brute.density(G, t, lam_full, C=C).rational()
        U = brute.symanzik_U(G, t, C)
        samples[ex] = D * U ** (n + H)
        cnt += 1
    if verbose:
        print('  %d density samples' % cnt, flush=True)
    poly = newton_interpolate(samples, N - 1, deg, total)
    full = {}
    for ex, c in poly.items():
        s = sum(ex)
        if s > total:
            raise ArithmeticError('degree overflow')
        full[(total - s,) + ex] = c
    # independent check at an unequal point
    import random
    rng = random.Random(7)
    t = [rng.randint(2, 9) for _ in range(N)]
    D = brute.density(G, t, lam_full, C=C).rational()
    U = brute.symanzik_U(G, t, C)
    val = sum(c * _mono(t, ex) for ex, c in full.items())
    assert val == D * U ** (n + H), 'interpolation check failed'
    return full


def _mono(t, ex):
    r = 1
    for x, k in zip(t, ex):
        if k:
            r *= x ** k
    return r


def numerator(G, lam, insertions=(), verbose=False):
    """returns (poly, q) with K = int poly(t)/U^q (projective)."""
    poly = bare_numerator(G, lam, verbose)
    if insertions:
        lam_nr = brute.sources(G, lam)
        H = len(lam_nr[0])
        lam_full = [list(v) for v in lam_nr] + [[-sum(v[i] for v in lam_nr) for i in range(H)]]
        Y = tree_current_polys(G, lam_full)
        # check the tree-flow currents against htk at an unequal point
        import random
        rng = random.Random(11)
        t = [rng.randint(2, 9) for _ in range(G.ne)]
        U = brute.symanzik_U(G, t, G.cycle_matrix())
        y = brute.currents(G, t, lam_nr)
        for e, i in set(insertions):
            val = sum(c * _mono(t, ex) for ex, c in Y[e][i].items())
            assert F(val) == y[e][i] * U, 'tree-flow current mismatch'
        for e, i in insertions:
            poly = pmul(poly, Y[e][i])
        # check one insertion product numerically
    n = G.loops
    H = len(brute.sources(G, lam)[0])
    q = n + H + len(insertions)
    return poly, q


def numerator_T1(G, lam, insertions=(), verbose=False):
    """T = 1, zero topological source: D = P / U^q with half-integer q.  The polynomial
    P = coeff * U^(q - 1/2) (htk Density: D = coeff * U^(-1/2)) is interpolated on the
    simplex grid of its total degree; q is the least half-integer making it polynomial
    (checked against an independent unequal point)."""
    n, N = G.loops, G.ne
    lam_nr = brute.sources(G, lam)
    H = len(lam_nr[0])
    if n * (H + 1) != N - 1:
        raise ValueError('not a T=1 graph: n (H+1) != N - 1')
    C = G.cycle_matrix()
    import random
    rng = random.Random(13)
    t_chk = [rng.randint(2, 9) for _ in range(N)]
    for k in range(0, 3 * n + 2 * H + 4):
        q = F(2 * k + 1, 2)
        total = n * q - N
        if total < 0 or total.denominator != 1:
            continue
        total = int(total)
        samples = {}
        for ex in lower_grid(N - 1, total, total):
            t = [1] + [1 + x for x in ex]
            d = brute.density(G, t, lam_nr, T=1, C=C)
            U = brute.symanzik_U(G, t, C)
            samples[ex] = d.coeff * U ** k
        poly = newton_interpolate(samples, N - 1, total, total)
        full = {(total - sum(ex),) + ex: c for ex, c in poly.items()}
        d = brute.density(G, t_chk, lam_nr, T=1, C=C)
        U = brute.symanzik_U(G, t_chk, C)
        if sum(c * _mono(t_chk, ex) for ex, c in full.items()) == d.coeff * U ** k:
            if verbose:
                print('  T=1 numerator: q = %s, %d samples' % (q, len(samples)), flush=True)
            break
    else:
        raise ArithmeticError('no half-integer q found')
    if insertions:
        H = len(lam_nr[0])
        lam_full = [list(v) for v in lam_nr] + [[-sum(v[i] for v in lam_nr) for i in range(H)]]
        Y = tree_current_polys(G, lam_full)
        for e, i in insertions:
            full = pmul(full, Y[e][i])
        q = q + len(insertions)
    return full, q
