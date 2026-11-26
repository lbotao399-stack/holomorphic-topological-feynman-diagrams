"""Evaluation-based Chevalley-Eilenberg layer (own code).

A k-cochain modulo conservation is represented by its values on k-tuples x_i = (c_i, p_i),
c_i an integer coordinate vector in the basis of the algebra, p_i in Z^H, with sum_i p_i = 0
(polynomials in p modulo the ideal (sum_i p_i) = polynomial functions on that hyperplane).
d: (d w)(x_0..x_k) = sum_{i<j} (-1)^{i+j} w([x_i,x_j], x_0, .., ^i, .., ^j, .., x_k),
[(c,p),(c',q)] = ([c,c'], p+q).
Cochain types:
  TraceCochain(mu, tau): alt over S_k of mu(p_1..p_{k-1}) tr_tau(X_1..X_k), tau a permutation of
    range(k) given by its cycles (tr of products along cycles), mu an exponent tuple over the
    (k-1) x H momentum components.
  TensorCochain(P, W): alt over S_k of P(p_1..p_k) W(c_1..c_k), W an integer numpy tensor.
Ranks: exact over Q (Fractions).
"""
from fractions import Fraction as Fr
import itertools
import numpy as np
import random


def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


PERMS = {k: [(s, perm_sign(s)) for s in itertools.permutations(range(k))] for k in range(1, 7)}


def cycles_of(s):
    k = len(s)
    seen, cyc = set(), []
    for i in range(k):
        if i in seen:
            continue
        c, j = [], i
        while j not in seen:
            seen.add(j)
            c.append(j)
            j = s[j]
        cyc.append(tuple(c))
    return tuple(cyc)


def trace_types(k):
    """all products of traces over the permutations of S_k (as cycle tuples)."""
    return [cycles_of(s) for s in itertools.permutations(range(k))]


def tr_eval(cyc, mats):
    val = 1
    for c in cyc:
        m = mats[c[0]]
        for j in c[1:]:
            m = m @ mats[j]
        val *= int(np.trace(m))
        if val == 0:
            return 0
    return val


def monomials(nmom, weights):
    """exponent tuples e[(i, comp)] over i < nmom, comp < H, degree weights[comp] in component comp."""
    H = len(weights)
    per = []
    for c in range(H):
        per.append([e for e in itertools.product(range(weights[c] + 1), repeat=nmom) if sum(e) == weights[c]])
    out = []
    for combo in itertools.product(*per):
        out.append(tuple(tuple(combo[c][i] for c in range(H)) for i in range(nmom)))
    return out


def mono_eval(mu, ps):
    v = 1
    for i, ex in enumerate(mu):
        for c, e in enumerate(ex):
            if e:
                v *= ps[i][c] ** e
    return v


class TraceCochain:
    def __init__(self, alg, mu, cyc, k):
        self.alg, self.mu, self.cyc, self.k = alg, mu, cyc, k

    def __call__(self, args):
        # args: list of (c, p, M) with M the matrix
        tot = 0
        for s, sg in PERMS[self.k]:
            ps = [args[s[i]][1] for i in range(self.k)]
            m = mono_eval(self.mu, ps)
            if m == 0:
                continue
            tot += sg * m * tr_eval(self.cyc, [args[s[i]][2] for i in range(self.k)])
        return tot


class TensorCochain:
    def __init__(self, P, W, k=4):
        self.P, self.W, self.k = P, W, k

    def __call__(self, args):
        k = self.k
        cs = [np.array(a[0], dtype=object) for a in args]
        # suffix contraction: suff[(s_0..s_{k-1})] = W(c_{s_0}, ..., c_{s_{k-1}})
        suff = {(): self.W}
        for L in range(1, k + 1):
            new = {}
            for key, T in suff.items():
                if len(key) != L - 1:
                    continue
                for i in range(k):
                    if i in key:
                        continue
                    new[(i,) + key] = np.tensordot(T, cs[i], axes=([T.ndim - 1], [0]))
            suff.update(new)
        tot = 0
        for s, sg in PERMS[k]:
            ps = [args[s[i]][1] for i in range(k)]
            pv = self.P(ps)
            if pv == 0:
                continue
            tot += sg * pv * int(suff[tuple(s)])
        return tot


def make_arg(alg, c, p):
    return (list(c), tuple(p), alg.mat(c))


def bracket_arg(alg, x, y):
    c = alg.bracket(x[0], y[0])
    p = tuple(a + b for a, b in zip(x[1], y[1]))
    return (c, p, alg.mat(c))


def d_eval(alg, cochain, args):
    """(d w)(x_0..x_k), args of length k+1."""
    n = len(args)
    tot = 0
    for i, j in itertools.combinations(range(n), 2):
        rest = [args[r] for r in range(n) if r not in (i, j)]
        tot += (-1) ** (i + j) * cochain([bracket_arg(alg, args[i], args[j])] + rest)
    return tot


def random_args(alg, rng, npts, H, crange=2, prange=3):
    cs = [[rng.randint(-crange, crange) for _ in range(alg.dim)] for _ in range(npts)]
    ps = [[rng.randint(-prange, prange) for _ in range(H)] for _ in range(npts - 1)]
    ps.append([-sum(p[c] for p in ps) for c in range(H)])
    return [make_arg(alg, c, p) for c, p in zip(cs, ps)]


def rank_Q(rows):
    """exact rank of a list of integer rows (Fraction elimination)."""
    A = [[Fr(x) for x in r] for r in rows]
    if not A:
        return 0
    n, m = len(A), len(A[0])
    r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        piv = A[r][c]
        for i in range(r + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / piv
                Ai, Ar = A[i], A[r]
                A[i] = [x - f * y for x, y in zip(Ai, Ar)]
        r += 1
        if r == n:
            break
    return r


def solve_Q(rows, rhs):
    """least-squares-free exact solve: find x with rows @ x = rhs if consistent (else None)."""
    n, m = len(rows), len(rows[0])
    A = [[Fr(x) for x in r] + [Fr(b)] for r, b in zip(rows, rhs)]
    piv_cols = []
    r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
        if r == n:
            break
    for i in range(r, n):
        if A[i][m] != 0:
            return None
    x = [Fr(0)] * m
    for i, c in enumerate(piv_cols):
        x[c] = A[i][m]
    return x


def pivots_Q(rows):
    """pivot column indices of the row-echelon form (exact)."""
    A = [[Fr(x) for x in r] for r in rows]
    n, m = len(A), len(A[0])
    r = 0
    piv = []
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        for i in range(r + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / pv
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
        if r == n:
            break
    return piv
