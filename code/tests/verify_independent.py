#!/usr/bin/env python3
"""Independent verification of code/htk.

Part I   (no htk) exact: Kirchhoff (nodal and mesh), random positive cycle bases, period matrix,
         T = 0 density (exterior algebra and determinant), T = 1 Pfaffian form from
         rho_n = det(M)^(-1/2) Pf(-1/2 dM M^-1 dM) (exterior algebra; n = 2 also eq-rho2 of S3),
         exact Dirichlet integration of the one-loop densities (polynomial on the simplex).
Part II  (no htk) double precision: densities S = D prod t, trapezoidal rule in log coordinates
         on the Schwinger simplex (flood fill), and for theta graphs a second rule
         (log coordinates of the three path resistances x Gauss-Jacobi on the path simplices).
Part III comparison with htk: brute.density, reduced densities g(xi), kernel_exact, transformation rules.

Definitions: root = nv-1, B_ev = delta_{v,e(0)} - delta_{v,e(1)} (v < root), t y = B phi, B^T y = lam,
C integer cycle basis with det[B^T; C] > 0, M = C t C^T, U = det M, psi^i_r = sum_e C_re y^i_e dt_e,
Omega = U^-H rho_n(M)^T ^ prod_{i=1..H} prod_{r=1..n} psi^i_r = D(t) iota_E dt_0 ^ ... ^ dt_{N-1},
K = int_{t > 0, sum t = 1} D(t) dt_0 ... dt_{N-2}.

usage: python3 verify_independent.py [section ...]
       sections: self pointwise reduced theta166 exact1 integrals integrals6 semianalytic0 semianalytic1
                 lattice6 transform n4;  'all' runs every section in this order (about 75 min in one process;
                 the log verify_independent.log was produced by one process per section).
"""
import sys
import os
import time
import math
import random
import itertools
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

P_ROOT = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + ""

# =====================================================================================
# Part I.  exact
# =====================================================================================


def fdet(A):
    A = [[Fr(x) for x in row] for row in A]
    n = len(A)
    d = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return d


def fsolve(A, Bm):
    n = len(A)
    m = len(Bm[0]) if n else 0
    W = [[Fr(x) for x in A[r]] + [Fr(x) for x in Bm[r]] for r in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if W[r][c] != 0), None)
        if p is None:
            raise ZeroDivisionError("singular")
        W[c], W[p] = W[p], W[c]
        inv = 1 / W[c][c]
        W[c] = [x * inv for x in W[c]]
        for r in range(n):
            if r != c and W[r][c] != 0:
                f = W[r][c]
                W[r] = [x - f * y for x, y in zip(W[r], W[c])]
    return [W[r][n:n + m] for r in range(n)]


def perm_sign(p):
    p = list(p)
    s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]
            p[i], p[j] = p[j], p[i]
            s = -s
    return s


class Net:
    def __init__(self, nv, edges):
        self.nv = nv
        self.edges = [tuple(e) for e in edges]
        self.N = len(self.edges)
        self.n = self.N - nv + 1
        self.g = nv - 1

    def B(self):
        Bm = [[0] * self.g for _ in range(self.N)]
        for e, (a, b) in enumerate(self.edges):
            if a < self.g:
                Bm[e][a] += 1
            if b < self.g:
                Bm[e][b] -= 1
        return Bm

    def tree(self, rng=None):
        order = list(range(self.N))
        if rng is not None:
            rng.shuffle(order)
        par = list(range(self.nv))

        def f(x):
            while par[x] != x:
                x = par[x]
            return x
        T = []
        for e in order:
            a, b = f(self.edges[e][0]), f(self.edges[e][1])
            if a != b:
                par[a] = b
                T.append(e)
        assert len(T) == self.nv - 1
        return sorted(T)

    def fund_cycles(self, T):
        """rows: cotree edge c with coefficient +1, closed through the tree path head(c) -> tail(c)."""
        adj = {v: [] for v in range(self.nv)}
        for e in T:
            a, b = self.edges[e]
            adj[a].append((b, e, 1))
            adj[b].append((a, e, -1))
        rows = []
        for c in range(self.N):
            if c in T:
                continue
            a, b = self.edges[c]
            prev = {b: None}
            st = [b]
            while st:
                x = st.pop()
                for (u, e, s) in adj[x]:
                    if u not in prev:
                        prev[u] = (x, e, s)
                        st.append(u)
            row = [0] * self.N
            row[c] = 1
            x = a
            while prev[x] is not None:
                px, e, s = prev[x]
                row[e] += s
                x = px
            rows.append(row)
        return rows

    def ntrees(self):
        Bm = self.B()
        return fdet([[sum(Bm[e][v] * Bm[e][w] for e in range(self.N)) for w in range(self.g)]
                     for v in range(self.g)])

    def detBC(self, C):
        Bm = self.B()
        return fdet([[Bm[e][v] for e in range(self.N)] for v in range(self.g)] + [list(r) for r in C])

    def fund_basis(self):
        """Fundamental cycles (entries 0, +-1) of the tree taken in edge order, first row sign fixed by
        det[B^T; C] > 0 (used for the double-precision integrations: well conditioned)."""
        C = self.fund_cycles(self.tree())
        if self.detBC(C) < 0:
            C[0] = [-x for x in C[0]]
        assert self.detBC(C) == self.ntrees()
        return C

    def cycle_basis(self, rng):
        """Random tree, random order and signs of the rows, random unimodular mixing,
        then the sign of the first row fixed by det[B^T; C] > 0."""
        T = self.tree(rng)
        C = self.fund_cycles(T)
        rng.shuffle(C)
        sg = [rng.choice((1, -1)) for _ in C]
        C = [[x * s for x in r] for r, s in zip(C, sg)]
        n = len(C)
        for _ in range(3 * n if n > 1 else 0):
            i, j = rng.sample(range(n), 2)
            k = rng.choice((-2, -1, 1, 2))
            C[i] = [a + k * b for a, b in zip(C[i], C[j])]
        d = self.detBC(C)
        if d < 0:
            C[0] = [-x for x in C[0]]
        assert self.detBC(C) == self.ntrees(), "not a unimodular basis"
        Bm = self.B()
        for r in C:
            for v in range(self.g):
                assert sum(r[e] * Bm[e][v] for e in range(self.N)) == 0
        return C


def kirchhoff_nodal(net, t, lam):
    Bm = net.B()
    g, N = net.g, net.N
    H = len(lam[0]) if lam else 0
    L = [[sum(Fr(Bm[e][v] * Bm[e][w]) / t[e] for e in range(N)) for w in range(g)] for v in range(g)]
    phi = fsolve(L, lam)
    return [[sum(Bm[e][v] * phi[v][i] for v in range(g)) / t[e] for i in range(H)] for e in range(N)]


def kirchhoff_mesh(net, t, lam, C, rng):
    Bm = net.B()
    T = net.tree(rng)
    H = len(lam[0])
    BT = [[Bm[e][v] for e in T] for v in range(net.g)]       # B_T^T
    y0T = fsolve(BT, lam)
    y0 = [[Fr(0)] * H for _ in range(net.N)]
    for k, e in enumerate(T):
        y0[e] = y0T[k]
    M = [[sum(Fr(a * b) * t[e] for e, (a, b) in enumerate(zip(ra, rb))) for rb in C] for ra in C]
    rhs = [[-sum(C[r][e] * t[e] * y0[e][i] for e in range(net.N)) for i in range(H)] for r in range(len(C))]
    ell = fsolve(M, rhs)
    return [[y0[e][i] + sum(C[r][e] * ell[r][i] for r in range(len(C))) for i in range(H)] for e in range(net.N)]


def check_laws(net, t, lam, y, C):
    Bm = net.B()
    H = len(lam[0])
    ok = all(sum(Bm[e][v] * y[e][i] for e in range(net.N)) == lam[v][i] for v in range(net.g) for i in range(H))
    ok &= all(sum(r[e] * t[e] * y[e][i] for e in range(net.N)) == 0 for r in C for i in range(H))
    return ok


def period(net, t, C):
    return [[sum(Fr(ra[e] * rb[e]) * t[e] for e in range(net.N)) for rb in C] for ra in C]


# ---- exterior algebra: {bitmask: Fraction}
def _wsign(ma, mb):
    s = 0
    m = mb
    while m:
        low = m & -m
        j = low.bit_length() - 1
        s += bin(ma >> (j + 1)).count("1")
        m ^= low
    return -1 if s & 1 else 1


def wedge(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            if ma & mb:
                continue
            m = ma | mb
            out[m] = out.get(m, 0) + _wsign(ma, mb) * ca * cb
    return {m: c for m, c in out.items() if c != 0}


def fadd(a, b, cb=1):
    out = dict(a)
    for m, c in b.items():
        out[m] = out.get(m, 0) + cb * c
    return {m: c for m, c in out.items() if c != 0}


def oneform(coeffs):
    return {1 << e: Fr(c) for e, c in enumerate(coeffs) if c != 0}


def pfaff(A, idx):
    """Pf(A) = sum_{j=2..n} (-1)^j A_{1j} Pf(A without 1, j)  (1-based j), entries 2-forms."""
    if not idx:
        return {0: Fr(1)}
    i0, rest = idx[0], idx[1:]
    out = {}
    for k, j in enumerate(rest):
        sub = rest[:k] + rest[k + 1:]
        out = fadd(out, wedge(A[i0][j], pfaff(A, sub)), 1 if k % 2 == 0 else -1)
    return out


def rho_form(net, t, C):
    """Pf(-1/2 dM M^-1 dM) as a form (rho_n = U^-1/2 times this); also the matrix A."""
    n = len(C)
    M = period(net, t, C)
    Mi = fsolve(M, [[Fr(int(i == j)) for j in range(n)] for i in range(n)])
    dM = [[oneform([C[r][e] * C[a][e] for e in range(net.N)]) for a in range(n)] for r in range(n)]
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Mi[a][b] != 0:
                        acc = fadd(acc, wedge(dM[r][a], dM[b][s]), -Mi[a][b] / 2)
            A[r][s] = acc
    for r in range(n):                      # antisymmetry of the matrix of 2-forms
        assert A[r][r] == {}
        for s in range(n):
            assert fadd(A[r][s], A[s][r]) == {}
    if n % 2:
        return {}, A, M
    return pfaff(A, list(range(n))), A, M


def rho2_S3(net, t, C):
    """S3 eq-rho2 numerator: -(M11 dM12 dM22 - M12 dM11 dM22 + M22 dM11 dM12)/(2 det M)."""
    M = period(net, t, C)
    d = {(r, s): oneform([C[r][e] * C[s][e] for e in range(net.N)]) for r in range(2) for s in range(2)}
    num = fadd(fadd(fadd({}, wedge(d[0, 1], d[1, 1]), M[0][0]), wedge(d[0, 0], d[1, 1]), -M[0][1]),
               wedge(d[0, 0], d[0, 1]), M[1][1])
    U = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return {m: -c / (2 * U) for m, c in num.items()}


def density_exact(net, t, lam, T, C, rng=None, y=None, full_check=True):
    """D(t) = q * U^(-1/2 * T); returns (q, U)."""
    t = [Fr(x) for x in t]
    N, n = net.N, len(C)
    H = len(lam[0]) if lam else 0
    assert n * (H + T) == N - 1
    if y is None:
        y = kirchhoff_nodal(net, t, lam) if H else []
    M = period(net, t, C)
    U = fdet(M)
    psis = [oneform([C[r][e] * y[e][i] for e in range(N)]) for i in range(H) for r in range(n)]
    form = {0: Fr(1)}
    if T == 1:
        form, _A, _M = rho_form(net, t, C)
    for p in psis:
        form = wedge(form, p)
    full = (1 << N) - 1
    vals = []
    e0s = range(N) if full_check else [N - 1]
    for e0 in e0s:
        c = form.get(full ^ (1 << e0), Fr(0))
        vals.append((-1) ** e0 * c / (t[e0] * U ** H))
    assert all(v == vals[0] for v in vals), "D depends on e0"
    q = vals[0]
    if T == 0 and full_check:
        # determinant formula: det[a; Psi] = (a.t) D U^H
        Psi = [[C[r][e] * y[e][i] for e in range(N)] for i in range(H) for r in range(n)]
        r_ = rng or random.Random(1)
        a = [Fr(r_.randint(-5, 9), r_.randint(1, 4)) for _ in range(N)]
        at = sum(x * te for x, te in zip(a, t))
        if at != 0:
            assert fdet([a] + Psi) == at * q * U ** H, "determinant identity"
    return q, U


def rand_fr(rng, lo=-9, hi=9, dmax=5, nonzero=False):
    while True:
        x = Fr(rng.randint(lo, hi), rng.randint(1, dmax))
        if x != 0 or not nonzero:
            return x


def rand_lam(rng, g, H):
    if H == 0:
        return [[] for _ in range(g)]
    while True:
        lam = [[rand_fr(rng, nonzero=True) for _ in range(H)] for _ in range(g)]
        if any(x != 0 for r in lam for x in r):
            return lam


def rand_t(rng, N):
    return [Fr(rng.randint(1, 40), rng.randint(1, 12)) for _ in range(N)]


# ---- exact integration of one-loop densities (polynomial on sum t = 1)
def symbolic_density_T0(net, lam, C, powers=None):
    """D(t) * prod (y^i_e)^p as a sympy rational function of t (T = 0), own Kirchhoff (nodal)."""
    N, g = net.N, net.g
    H = len(lam[0])
    ts = sp.symbols("t0:%d" % N, positive=True)
    Bm = sp.Matrix(net.B())
    Tinv = sp.diag(*[1 / x for x in ts])
    L = Bm.T * Tinv * Bm
    lamM = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in lam])
    phi = L.LUsolve(lamM)
    Y = (Tinv * Bm * phi).applyfunc(sp.cancel)
    Cm = sp.Matrix(C)
    M = Cm * sp.diag(*ts) * Cm.T
    U = sp.expand(M.det())
    n = len(C)
    Psi = sp.Matrix([[C[r][e] * Y[e, i] for e in range(N)] for i in range(H) for r in range(n)])
    W = sp.Matrix.vstack(sp.ones(1, N), Psi)
    D = sp.cancel(W.det(method="berkowitz") / (sum(ts) * U ** H))
    if powers:
        for (e, i), p in powers.items():
            D = D * Y[e, i] ** p
        D = sp.cancel(D)
    return D, ts, U


def simplex_integral_exact(D, ts, U):
    """int_{sum t = 1} D dt_0..dt_{N-2} for D = P / (c * U^k) with U = sum t (one loop)."""
    N = len(ts)
    s = sum(ts)
    assert sp.expand(U - s) == 0
    num, den = sp.fraction(sp.cancel(D))
    k = sp.degree(sp.Poly(den, *ts), ts[0])
    c = sp.cancel(den / s ** k)
    assert c.is_number, "denominator is not a power of sum t"
    poly = sp.Poly(sp.expand(num), *ts)
    tot = sp.Integer(0)
    for mon, coef in poly.terms():
        val = sp.Integer(1)
        for a in mon:
            val *= sp.factorial(a)
        tot += coef * val / sp.factorial(N - 1 + sum(mon))
    r = sp.Rational(tot / c)
    return Fr(int(r.p), int(r.q))


# =====================================================================================
# Part II.  double precision
# =====================================================================================


class FDens:
    """Vectorized S(t) = D(t) prod_e t_e (scale invariant) for a fixed positive basis C."""

    def __init__(self, net, lam, T, C, powers=None):
        self.net, self.T = net, T
        self.C = np.array(C, dtype=float)
        self.n, self.N = len(C), net.N
        self.H = len(lam[0])
        assert self.n * (self.H + T) == self.N - 1
        assert T == 0 or self.n == 2
        Bm = net.B()
        Tr = net.tree(random.Random(7))
        BT = [[Bm[e][v] for e in Tr] for v in range(net.g)]
        y0T = fsolve(BT, lam)
        y0 = np.zeros((self.N, self.H))
        for k, e in enumerate(Tr):
            y0[e] = [float(x) for x in y0T[k]]
        self.y0 = y0
        self.powers = powers or {}
        self.nbad = 0

    def currents(self, t):
        C = self.C
        M = np.einsum("re,pe,se->prs", C, t, C)
        rhs = -np.einsum("re,pe,ei->pri", C, t, self.y0)
        try:
            ell = np.linalg.solve(M, rhs)
        except np.linalg.LinAlgError:
            # numerically singular M (extreme Schwinger ratios): those points are set to zero
            ell = np.zeros_like(rhs)
            for p in range(M.shape[0]):
                try:
                    ell[p] = np.linalg.solve(M[p], rhs[p])
                except np.linalg.LinAlgError:
                    ell[p] = np.nan
        y = self.y0[None, :, :] + np.einsum("re,pri->pei", C, ell)
        return y, M

    def currents_unused(self, t):
        C = self.C
        M = np.einsum("re,pe,se->prs", C, t, C)
        rhs = -np.einsum("re,pe,ei->pri", C, t, self.y0)
        ell = np.linalg.solve(M, rhs)
        y = self.y0[None, :, :] + np.einsum("re,pri->pei", C, ell)
        return y, M

    def __call__(self, t):
        t = t / t.max(axis=1, keepdims=True)
        P, N, n, H = t.shape[0], self.N, self.n, self.H
        np.seterr(all="ignore")
        y, M = self.currents(t)
        U = np.linalg.det(M)
        Psi = np.einsum("re,pei->pire", self.C, y).reshape(P, H * n, N)
        st = t.sum(axis=1)
        if self.T == 0:
            W = np.concatenate([np.ones((P, 1, N)), Psi], axis=1)
            D = np.linalg.det(W) / (st * U ** H)
        else:
            try:
                Mi = np.linalg.inv(M)
            except np.linalg.LinAlgError:
                Mi = np.full_like(M, np.nan)
                for p in range(M.shape[0]):
                    try:
                        Mi[p] = np.linalg.inv(M[p])
                    except np.linalg.LinAlgError:
                        pass
            Pm = np.einsum("ae,pab,bf->pef", self.C, Mi, self.C)
            beta = self.C[0][None, :, None] * Pm * self.C[1][None, None, :]
            alpha = -0.5 * (beta - np.transpose(beta, (0, 2, 1)))
            acc = np.zeros(P)
            for e in range(N):
                for f in range(e + 1, N):
                    a = alpha[:, e, f]
                    if not np.any(a):
                        continue
                    W = np.zeros((P, N, N))
                    W[:, 0, :] = 1.0
                    W[:, 1, e] = 1.0
                    W[:, 2, f] = 1.0
                    W[:, 3:, :] = Psi
                    acc += a * np.linalg.det(W)
            D = acc / (st * U ** H * np.sqrt(U))
        for (e, i), p in self.powers.items():
            D = D * y[:, e, i] ** p
        with np.errstate(all="ignore"):
            out = D * np.prod(t, axis=1)
        bad = ~np.isfinite(out)
        if bad.any():
            self.nbad += int(bad.sum())
            out[bad] = 0.0
        return out


def flood_lattice(func, dim, h, tol, shift=0.0, batch=400000, maxpts=2e8):
    """h^dim * sum over h(Z^dim + shift) of func(s), restricted by flood fill from 0 to
    {|f| > tol * running max} and its nearest neighbours."""
    off = np.vstack([np.eye(dim, dtype=np.int64), -np.eye(dim, dtype=np.int64)])
    base = np.int64(512)
    half = np.int64(256)

    def key(k):
        kk = k + half
        out = np.zeros(k.shape[0], dtype=np.int64)
        for j in range(dim):
            out = out * base + kk[:, j]
        return out
    frontier = np.zeros((1, dim), dtype=np.int64)
    seen = set(key(frontier).tolist())
    total, comp, vmax, npts = 0.0, 0.0, 0.0, 0
    while frontier.shape[0]:
        keep_all = []
        for b0 in range(0, frontier.shape[0], batch):
            fr = frontier[b0:b0 + batch]
            v = func(h * (fr + shift))
            # Kahan summation
            ysum = float(np.sum(v)) - comp
            tt = total + ysum
            comp = (tt - total) - ysum
            total = tt
            vmax = max(vmax, float(np.max(np.abs(v))))
            keep_all.append(fr[np.abs(v) > tol * vmax])
        npts += frontier.shape[0]
        if npts > maxpts:
            raise RuntimeError("too many points")
        keep = np.concatenate(keep_all)
        if keep.shape[0] == 0:
            break
        cand = (keep[:, None, :] + off[None, :, :]).reshape(-1, dim)
        assert np.abs(cand).max() < half - 1
        ks = key(cand)
        ks, idx = np.unique(ks, return_index=True)
        new = np.fromiter((k not in seen for k in ks.tolist()), dtype=bool, count=ks.shape[0])
        seen.update(ks[new].tolist())
        frontier = cand[idx[new]]
    return total * h ** dim, npts


def flood_lattice_np(func, dim, h, tol, shift=0.0, batch=400000, maxpts=2e8):
    """Same rule as flood_lattice, visited set kept as a sorted int64 array (memory 8 bytes/point)."""
    off = np.vstack([np.eye(dim, dtype=np.int64), -np.eye(dim, dtype=np.int64)])
    base, half = np.int64(512), np.int64(256)

    def key(k):
        kk = k + half
        out = np.zeros(k.shape[0], dtype=np.int64)
        for j in range(dim):
            out = out * base + kk[:, j]
        return out
    frontier = np.zeros((1, dim), dtype=np.int64)
    seen = key(frontier)
    total, comp, vmax, npts = 0.0, 0.0, 0.0, 0
    while frontier.shape[0]:
        keep_all = []
        for b0 in range(0, frontier.shape[0], batch):
            fr = frontier[b0:b0 + batch]
            v = func(h * (fr + shift))
            ysum = float(np.sum(v)) - comp
            tt = total + ysum
            comp = (tt - total) - ysum
            total = tt
            vmax = max(vmax, float(np.max(np.abs(v))))
            keep_all.append(fr[np.abs(v) > tol * vmax])
        npts += frontier.shape[0]
        if npts > maxpts:
            raise RuntimeError("too many points")
        keep = np.concatenate(keep_all)
        if keep.shape[0] == 0:
            break
        cand = (keep[:, None, :] + off[None, :, :]).reshape(-1, dim)
        assert np.abs(cand).max() < half - 1
        ks = key(cand)
        ks, idx = np.unique(ks, return_index=True)
        pos = np.searchsorted(seen, ks)
        pos[pos == seen.shape[0]] = 0
        new = seen[pos] != ks
        seen = np.union1d(seen, ks[new])
        frontier = cand[idx[new]]
    return total * h ** dim, npts


def log_integrand(fd):
    N = fd.N

    def f(s):
        t = np.exp(np.concatenate([s, np.zeros((s.shape[0], 1))], axis=1))
        return fd(t)
    return f


def simplex_rule(m, k):
    """Nodes (Q, m+1) and weights for int over {x >= 0, sum x = 1} in Lebesgue measure of
    x_1..x_m (volume 1/m!), collapsed Gauss-Jacobi with k points per direction."""
    from scipy.special import roots_jacobi
    if m == 0:
        return np.ones((1, 1)), np.ones(1)
    xs, ws = roots_jacobi(k, m - 1, 0)
    u = (1 + xs) / 2
    w = ws / 2 ** m
    Xs, Ws = simplex_rule(m - 1, k)
    nodes, weights = [], []
    for ui, wi in zip(u, w):
        nodes.append(np.concatenate([(1 - ui) * Xs, np.full((Xs.shape[0], 1), ui)], axis=1))
        weights.append(wi * Ws)
    return np.concatenate(nodes), np.concatenate(weights)


def theta_chain_integral(fd, paths, h, tol, k):
    """Theta graph with paths (lists of edges): t_{P,j} = R_P x_{P,j}, R_last = 1 (gauge),
    K = int ds_1 ds_2 R_1^{a_1} R_2^{a_2} int_{simplices} D(t) dx, R_P = e^{s_P}."""
    rules = [simplex_rule(len(pth) - 1, k) for pth in paths]
    grids = list(itertools.product(*[range(len(r[1])) for r in rules]))
    N = fd.N
    Q = len(grids)
    X = np.zeros((Q, N))
    Wq = np.ones(Q)
    for qi, combo in enumerate(grids):
        for (nodes, wts), pth, j in zip(rules, paths, combo):
            for col, e in enumerate(pth):
                X[qi, e] = nodes[j, col]
            Wq[qi] *= wts[j]
    a = [len(p) for p in paths]
    mask = [np.array([e in pth for e in range(N)]) for pth in paths]

    def f(s):
        P = s.shape[0]
        R = np.concatenate([np.exp(s), np.ones((P, 1))], axis=1)       # (P, 3)
        scale = np.zeros((P, N))
        for pi in range(3):
            scale[:, mask[pi]] = R[:, pi:pi + 1]
        t = (scale[:, None, :] * X[None, :, :]).reshape(P * Q, N)
        S = fd(t).reshape(P, Q)
        D = S / np.prod(t, axis=1).reshape(P, Q)
        # undo the normalization of fd: S is scale invariant, D(c t) = c^-N D(t); fd returns
        # D(t) prod t evaluated consistently, so D here is D at the unnormalized t
        inner = D @ Wq
        return inner * R[:, 0] ** a[0] * R[:, 1] ** a[1]
    return flood_lattice(f, 2, h, tol)


def digits(a, b):
    if a == b:
        return 16.0
    return -math.log10(abs(a - b) / max(abs(a), abs(b)))


# =====================================================================================
# Part III.  comparison with htk
# =====================================================================================
LOG = []
DISC = []


def say(s=""):
    print(s, flush=True)
    LOG.append(s)


def record(ok, label, detail=""):
    say(("PASS " if ok else "FAIL ") + label + ("  " + detail if detail else ""))
    if not ok:
        DISC.append(label + ("  " + detail if detail else ""))


def import_htk():
    sys.path.insert(0, os.path.join(P_ROOT, "code"))
    import htk
    return htk


def cases(htk):
    """(label, Graph, H, T, interior or None)."""
    out = []
    out.append(("triangle", htk.polygon(3), 2, 0))
    out.append(("fan2", htk.fan(2), 2, 0))
    out.append(("box", htk.box(), 3, 0))
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        out.append(("theta%s H3T0" % (pqr,), htk.theta(*pqr)[0], 3, 0))
    out.append(("theta(1, 2, 2) H1T1", htk.theta(1, 2, 2)[0], 1, 1))
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        out.append(("theta%s H2T1" % (pqr,), htk.theta(*pqr)[0], 2, 1))
    out.append(("theta(1, 6, 6) H5T1", htk.theta(1, 6, 6)[0], 5, 1))
    return out


def extra_graphs(htk):
    """Non-series-parallel and n = 4 graphs for pointwise checks."""
    prism = htk.Graph(6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)], "prism")
    out = [("prism H2T0", prism, 2, 0), ("prism H1T1", prism, 1, 1),
           ("graph_A H3T0", htk.graph_A(), 3, 0), ("graph_B H3T0", htk.graph_B(), 3, 0),
           ("ladder(2) H3T0", htk.ladder(2)[0], 3, 0),
           ("banana(3) H0T1", htk.banana(3), 0, 1), ("banana(5) H0T1", htk.banana(5), 0, 1)]
    return out


def net_of(G):
    return Net(G.nv, G.edges)


def sec_self(rng):
    say("== self-tests of Part I/II (no htk)")
    # Pfaffian of scalars: Pf^2 = det for random antisymmetric 4x4, 6x6
    for n in (2, 4, 6):
        A = [[Fr(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                A[i][j] = rand_fr(rng)
                A[j][i] = -A[i][j]
        Af = [[{0: A[i][j]} if A[i][j] != 0 else {} for j in range(n)] for i in range(n)]
        pf = pfaff(Af, list(range(n))).get(0, Fr(0))
        record(pf ** 2 == fdet(A), "Pf^2 = det, n=%d" % n, "Pf = %s" % pf)
    A = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]
    Af = [[{0: Fr(x)} if x else {} for x in r] for r in A]
    record(pfaff(Af, [0, 1, 2, 3]).get(0) == 1, "Pf(J) = +1 (a12 a34 term +)")
    # wedge sign
    a = oneform([1, 0, 0]); b = oneform([0, 1, 0])
    record(wedge(a, b) == {3: 1} and wedge(b, a) == {3: -1}, "dt0^dt1 = -dt1^dt0")
    # simplex rule
    for m in (1, 2, 3, 5):
        X, W = simplex_rule(m, 6)
        ok = abs(W.sum() - 1 / math.factorial(m)) < 1e-14
        ex = [2, 1] + [0] * (m - 1)
        val = (W * np.prod(X ** np.array(ex), axis=1)).sum()
        exact = math.factorial(2) * 1 / math.factorial(m + 3)
        record(ok and abs(val - exact) < 1e-14, "simplex rule m=%d" % m, "%.15g vs %.15g" % (val, exact))
    # lattice rule on a known integral: int_{sum t=1, N=3} t0^1 t1^0 t2^2 = 1! 2!/5! via log coords
    def f(s):
        t = np.exp(np.concatenate([s, np.zeros((s.shape[0], 1))], axis=1))
        tt = t / t.sum(axis=1, keepdims=True)
        D = tt[:, 0] * tt[:, 2] ** 2 / t.sum(axis=1) ** 3
        return D * np.prod(t, axis=1)
    v, npts = flood_lattice(f, 2, 0.4, 1e-13)
    record(abs(v - 2 / 120) < 1e-12, "log-lattice Dirichlet test", "%.15g vs %.15g (%d pts)" % (v, 2 / 120, npts))


def sec_pointwise(htk, rng, which=None, npts=4):
    say("== pointwise densities: own exact D(t) vs htk.brute.density (random rational t, random positive bases)")
    allc = cases(htk) + extra_graphs(htk)
    for (label, G, H, T) in allc:
        if which and label not in which:
            continue
        net = net_of(G)
        t0 = time.time()
        bad = []
        vals = []
        for k in range(npts):
            lam = rand_lam(rng, net.g, H)
            t = rand_t(rng, net.N)
            C1 = net.cycle_basis(rng)
            C2 = net.cycle_basis(rng)
            y = kirchhoff_nodal(net, t, lam)
            y2 = kirchhoff_mesh(net, t, lam, C1, rng)
            if y != y2 or not check_laws(net, t, lam, y, C1):
                bad.append("Kirchhoff")
            q1, U1 = density_exact(net, t, lam, T, C1, rng, y, full_check=(net.N <= 9))
            q2, U2 = density_exact(net, t, lam, T, C2, rng, y, full_check=False)
            if (q1, U1) != (q2, U2):
                bad.append("basis dependence %s %s" % (q1, q2))
            if T == 1 and len(C1) == 2:
                rf, _A, _M = rho_form(net, t, C1)
                if rf != rho2_S3(net, t, C1):
                    bad.append("eq-rho2")
            De = htk.brute.density(G, t, lam if H else None, T)
            okv = (De.coeff == q1) and (De.radicand == (U1 if T else 1))
            if not okv:
                bad.append("t=%s: htk %r, own %s * U^-%d/2 (U=%s)" % ([str(x) for x in t], De, q1, T, U1))
            vals.append(q1)
            # float density vs exact at the same point
            if T == 0 or len(C1) == 2:
                fd = FDens(net, lam, T, C1)
                tf = np.array([[float(x) for x in t]])
                Sf = fd(tf)[0]
                Sx = float(q1) / (math.sqrt(float(U1)) if T else 1.0) * float(np.prod([float(x) for x in t]))
                if abs(Sf - Sx) > 1e-9 * max(1.0, abs(Sx)):
                    bad.append("float density %g vs %g" % (Sf, Sx))
        record(not bad, "pointwise %s (%d points)" % (label, npts),
               ("; ".join(bad[:3]) if bad else "e.g. D = %s%s" % (vals[0], " * U^-1/2" if T else "")) +
               " [%.1fs]" % (time.time() - t0))


def ratio_jacobian(red):
    """|det[t; dt/dxi]| as a function of the independent ratios, own sympy differentiation of t_expr."""
    texpr = red.t_expr()
    subs = {}
    for blk in red.blocks:
        subs[sp.Symbol(blk.bary[0])] = 1 - sum(sp.Symbol(b) for b in blk.indep)
    coords = [sp.Symbol(c) for c in red.coords]
    tx = [sp.sympify(e).subs(subs) for e in texpr]
    rows = [tx] + [[sp.diff(e, c) for e in tx] for c in coords]
    J = sp.Matrix(rows)
    return tx, J, coords


def sec_reduced(htk, rng, which=None, npts=3, prefer=None):
    say("== reduced densities: htk g(xi) vs own D(t(xi)) |det[t; dt/dxi]| (own Jacobian), both plans;"
        " K for both plans and three source formats")
    for (label, G, H, T) in cases(htk):
        if which and label not in which:
            continue
        if label.startswith("theta(1, 6, 6)") and not (which and label in which):
            continue
        net = net_of(G)
        if prefer is None:
            lam = rand_lam(rng, net.g, H)
            Ks = [htk.sp_reduce(G, prefer=p).integrate(lam, T) for p in ("small", "large")]
            la = src_all(lam, net.g)
            K_all = htk.kernel_exact(G, la, T)
            order = list(range(net.g))
            rng.shuffle(order)
            K_dict = htk.kernel_exact(G, {v: lam[v] for v in order}, T)
            record(Ks[0] == Ks[1] and K_all == Ks[0] and K_dict == Ks[0] and Ks[0] != 0,
                   "plans/source formats %s" % label, "K(small) = %s, K(large) = %s, K(all-vertex list) = %s,"
                   " K(dict) = %s" % (Ks[0], Ks[1], K_all, K_dict))
        for pf in ([prefer] if prefer else ["small", "large"]):
            sec_reduced_one(htk, rng, label, G, H, T, pf, npts)


def sec_reduced_one(htk, rng, label, G, H, T, pf, npts):
        net = net_of(G)
        t0 = time.time()
        red = htk.sp_reduce(G, prefer=pf)
        lam = rand_lam(rng, net.g, H)
        dens = red.density(lam, T, check=0)
        tx, J, coords = ratio_jacobian(red)
        bad = []
        C = net.cycle_basis(rng)
        for k in range(npts):
            pt = red.random_point(rng)
            vals = tuple(pt[c] for c in red.coords)
            sub = {c: sp.Rational(v.numerator, v.denominator) for c, v in zip(coords, vals)}
            t = [Fr(int(sp.Rational(e.subs(sub)).p), int(sp.Rational(e.subs(sub)).q)) for e in tx]
            t_htk = red.tmap(pt)
            if t != t_htk:
                bad.append("t(xi) mismatch")
            Jr = sp.Rational(J.subs(sub).det())
            Jv = Fr(int(Jr.p), int(Jr.q))
            q, U = density_exact(net, t, lam, T, C, rng, full_check=False)
            rat, rad = dens.evaluate(vals)
            if T == 0:
                ok = (rad == 1 and rat == q * abs(Jv))
            else:
                gm = q * abs(Jv)
                ok = ((rat > 0) - (rat < 0) == (gm > 0) - (gm < 0)) and rat * rat * rad == gm * gm / U
            if not ok:
                bad.append("xi=%s: htk %s*sqrt(%s), own %s*|J|=%s (U=%s)" % ([str(v) for v in vals], rat, rad, q, Jv, U))
        record(not bad, "reduced density %s (%s plan, %d points)" % (label, pf, npts),
               ("; ".join(bad[:2]) if bad else "exps %s" % {k: str(v) for k, v in dens.exponents().items()}) +
               " [%.1fs]" % (time.time() - t0))


def sec_exact1(htk, rng):
    say("== one-loop cases: exact integral of own polynomial density on the simplex vs htk.kernel_exact")
    for (label, G, H) in [("triangle", htk.polygon(3), 2), ("box", htk.box(), 3), ("square H3", htk.polygon(4), 3),
                          ("pentagon H4", htk.polygon(5), 4)]:
        net = net_of(G)
        for trial in range(2):
            lam = rand_lam(rng, net.g, H)
            C = net.cycle_basis(rng)
            D, ts, U = symbolic_density_T0(net, lam, C)
            Kown = simplex_integral_exact(D, ts, U)
            Kh = htk.kernel_exact(G, lam, 0)
            record(Kown == Kh, "exact K %s trial %d" % (label, trial), "own %s, htk %s" % (Kown, Kh))
        # insertions
        lam = rand_lam(rng, net.g, H)
        C = net.cycle_basis(rng)
        for powers in [{(0, 0): 1}, {(1, H - 1): 2}, {(0, 0): 1, (net.N - 1, 1): 1}]:
            D, ts, U = symbolic_density_T0(net, lam, C, powers)
            Kown = simplex_integral_exact(D, ts, U)
            Kh = htk.kernel_exact(G, lam, 0, powers=powers)
            record(Kown == Kh, "exact moment %s powers %s" % (label, powers), "own %s, htk %s" % (Kown, Kh))
        # polygon formula K = (-1)^d det(lam_0..lam_{d-1})/d!  (d = H, sources at vertices 0..d-1)
        if net.N == H + 1:
            d = H
            Kf = Fr((-1) ** d) * fdet([[lam[v][i] for v in range(d)] for i in range(d)]) / math.factorial(d)
            Kh = htk.kernel_exact(G, lam, 0)
            record(Kf == Kh, "polygon formula %s" % label, "%s vs %s" % (Kf, Kh))


def integrate_case(net, lam, T, C, hs, tol):
    fd = FDens(net, lam, T, C)
    f = log_integrand(fd)
    out = []
    for k, h in enumerate(hs):
        t0 = time.time()
        v, npts = flood_lattice(f, net.N - 1, h, tol)
        vs = flood_lattice(f, net.N - 1, h, tol, shift=0.5)[0] if k == 0 else float("nan")
        out.append((h, v, vs, npts, time.time() - t0))
    out.append(fd.nbad)
    return out


def kval(K):
    return float(K)


def sec_integrals(htk, rng):
    say("== numerical K(lam; 0): own log-lattice on the simplex (dims 2..4) vs htk.kernel_exact")
    spec = [("triangle", htk.polygon(3), 2, 0, [0.5, 0.35], 1e-14),
            ("fan2", htk.fan(2), 2, 0, [0.5, 0.4], 1e-12),
            ("box", htk.box(), 3, 0, [0.5, 0.4], 1e-13),
            ("theta(1, 2, 2) H1T1", htk.theta(1, 2, 2)[0], 1, 1, [0.5, 0.4], 1e-12),
            ("banana(3) H0T1", htk.banana(3), 0, 1, [0.5, 0.4], 1e-13)]
    for (label, G, H, T, hs, tol) in spec:
        net = net_of(G)
        for trial in range(2):
            lam = rand_lam(rng, net.g, H)
            C = net.fund_basis()
            Kh = htk.kernel_exact(G, lam if H else None, T)
            res = integrate_case(net, lam, T, C, hs, tol)
            nbad = res.pop()
            h, v, vs, npts, dt = res[-1]
            h0, v0, vs0 = res[0][0], res[0][1], res[0][2]
            dg = digits(v, kval(Kh))
            say("   %s: htk %s = %.12g; lattice h=%.2f: %.12g (shifted by h/2: %.12g), h=%.2f: %.12g (%d pts, %.0fs); %d non-finite set to 0"
                % (label, Kh, kval(Kh), h0, v0, vs0, h, v, npts, dt, nbad))
            record(dg >= 6, "integral %s trial %d" % (label, trial), "%.1f digits agree" % dg)


def theta_paths(G, pqr):
    p, q, r = pqr
    return [list(range(0, p)), list(range(p, p + q)), list(range(p + q, p + q + r))]


def sec_integrals6(htk, rng, ntrial=2):
    say("== numerical K(lam; 0) of theta graphs: own path-resistance rule (log lattice in R_A, R_B, R_C = 1;"
        " Gauss-Jacobi on path simplices) vs htk.kernel_exact")
    spec = [((1, 2, 2), 1, 1), ((1, 3, 3), 3, 0), ((2, 2, 3), 3, 0), ((1, 3, 3), 2, 1), ((2, 2, 3), 2, 1)]
    for (pqr, H, T) in spec:
        G = htk.theta(*pqr)[0]
        net = net_of(G)
        paths = theta_paths(G, pqr)
        for trial in range(ntrial):
            lam = rand_lam(rng, net.g, H)
            C = net.fund_basis()
            t0 = time.time()
            Kh = htk.kernel_exact(G, lam, T)
            tk = time.time() - t0
            fd = FDens(net, lam, T, C)
            res = []
            for (h, k) in [(0.3, 4), (0.25, 5)]:
                t1 = time.time()
                v, npts = theta_chain_integral(fd, paths, h, 1e-13, k)
                res.append((h, k, v, npts, time.time() - t1))
            nbad = fd.nbad
            dg = digits(res[-1][2], kval(Kh))
            say("   theta%s H=%d T=%d: htk %s = %.12g [%.1fs]; rule (h,k)=(%.2f,%d): %.12g, (%.2f,%d): %.12g (%d outer pts, %.0fs, %d non-finite set to 0)"
                % (pqr, H, T, Kh, kval(Kh), tk, res[0][0], res[0][1], res[0][2], res[1][0], res[1][1], res[1][2],
                   res[1][3], res[1][4], nbad))
            record(dg >= 6, "integral theta%s H%dT%d trial %d" % (pqr, H, T, trial), "%.1f digits agree" % dg)


def theta_semianalytic(net, lam, T, C, paths, dps=30, check_rng=None, maxdegree=8):
    """Theta graph, n = 2: t_{P,j} = R_P x_{P,j} (R_last = 1).  Own exact density in the sparse
    polynomial ring QQ[R_0, R_1, x]: with a = sum_e dt_e,
        a ^ Omega' = G dt_0..dt_{N-1},  Omega' = (U alpha)^T ^ prod (U psi),  alpha = -1/2 dM M^-1 dM [1,2],
        D = G / (sum t * U^(1 + nH + H) * U^(T/2)),
    exact Dirichlet integration of G over the path simplices, then the 2-dim integral over
    R_0, R_1 in (0, inf) by mpmath (tanh-sinh, log variables):
        K = int dR_0 dR_1 R_0^{a_0-1} R_1^{a_1-1} int D dx."""
    import mpmath as mp
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    N, n = net.N, len(C)
    H = len(lam[0])
    assert n == 2
    names = ["R0", "R1"]
    xidx = []
    for pi, pth in enumerate(paths):
        xidx.append([])
        for j in range(1, len(pth)):
            xidx[-1].append(len(names))
            names.append("x%d_%d" % (pi, j))
    Rg, *gens = ring(",".join(names), QQ)
    Rv = [gens[0], gens[1], Rg.one]
    t = [None] * N
    for pi, pth in enumerate(paths):
        xv = [gens[k] for k in xidx[pi]]
        x1 = Rg.one - sum(xv, Rg.zero)
        for j, e in enumerate(pth):
            t[e] = Rv[pi] * (x1 if j == 0 else xv[j - 1])
    Bm = net.B()
    Tr = net.tree(random.Random(7))
    BT = [[Bm[e][v] for e in Tr] for v in range(net.g)]
    y0T = fsolve(BT, lam)
    y0 = [[QQ(0)] * H for _ in range(N)]
    for k, e in enumerate(Tr):
        y0[e] = [QQ(x.numerator, x.denominator) for x in y0T[k]]
    M = [[sum((t[e] * (C[r][e] * C[s][e]) for e in range(N)), Rg.zero) for s in range(n)] for r in range(n)]
    U = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    adj = [[M[1][1], -M[0][1]], [-M[1][0], M[0][0]]]
    rhs = [[-sum((t[e] * (C[r][e] * y0[e][i]) for e in range(N)), Rg.zero) for i in range(H)] for r in range(n)]
    ell = [[sum((adj[r][s] * rhs[s][i] for s in range(n)), Rg.zero) for i in range(H)] for r in range(n)]   # U * l
    Uy = [[U * y0[e][i] + sum((ell[r][i] * C[r][e] for r in range(n)), Rg.zero) for i in range(H)] for e in range(N)]
    psis = [{1 << e: Uy[e][i] * C[r][e] for e in range(N) if C[r][e] != 0 and Uy[e][i] != 0}
            for i in range(H) for r in range(n)]
    form = {1 << e: Rg.one for e in range(N)}                  # a = sum_e dt_e
    if T == 1:
        Pp = [[sum((adj[a_][b_] * (C[a_][e] * C[b_][f]) for a_ in range(n) for b_ in range(n)), Rg.zero)
               for f in range(N)] for e in range(N)]
        alpha = {}
        for e in range(N):
            for f in range(e + 1, N):
                v = (Pp[e][f] * (C[0][e] * C[1][f]) - Pp[f][e] * (C[0][f] * C[1][e])) * QQ(-1, 2)
                if v != 0:
                    alpha[(1 << e) | (1 << f)] = v
        form = wedge(form, alpha)
    for p in psis:
        form = wedge(form, p)
    G = form.get((1 << N) - 1, Rg.zero)
    k_U = (1 if T else 0) + n * H + H
    stR = Rv[0] + Rv[1] + Rg.one
    if check_rng is not None:
        for _ in range(3):
            vals = [QQ(check_rng.randint(1, 20), check_rng.randint(1, 7)) for _ in range(2)]
            vx = [None] * (len(names) - 2)
            for pi in range(len(paths)):
                ws = [check_rng.randint(1, 9) for _ in range(len(xidx[pi]) + 1)]
                for j, k in enumerate(xidx[pi]):
                    vx[k - 2] = QQ(ws[j + 1], sum(ws))
            pt = vals + vx
            tq = [te(*pt) if te.degree() >= 0 and te != 0 else QQ(0) for te in t]
            tf = [Fr(int(v.numerator), int(v.denominator)) for v in tq]
            q, Uq = density_exact(net, tf, lam, T, C, full_check=False)
            Gv, Uv, sv = G(*pt), U(*pt), stR(*pt)
            lhs = Gv / (sv * Uv ** k_U)
            assert Fr(int(lhs.numerator), int(lhs.denominator)) == q, "semianalytic density"
            assert Fr(int(Uv.numerator), int(Uv.denominator)) == Uq
    # exact integration over the path simplices: int prod_{j>=2} x_j^b_j = prod b_j! / (m + sum b)!
    I = {}
    for mon, coef in G.terms():
        val = coef
        for pi in range(len(paths)):
            bs = [mon[k] for k in xidx[pi]]
            m = len(bs)
            num = 1
            for b in bs:
                num *= math.factorial(b)
            val = val * QQ(num, math.factorial(m + sum(bs)))
        key = (mon[0], mon[1])
        I[key] = I.get(key, QQ(0)) + val
    Ipoly = [(a, b, Fr(int(c.numerator), int(c.denominator))) for (a, b), c in I.items() if c != 0]
    Upoly = [(mon[0], mon[1], Fr(int(c.numerator), int(c.denominator))) for mon, c in U.terms()]
    a0, a1 = len(paths[0]) - 1, len(paths[1]) - 1
    mp.mp.dps = dps
    Ipm = [(a, b, mp.mpf(c.numerator) / c.denominator) for a, b, c in Ipoly]
    Upm = [(a, b, mp.mpf(c.numerator) / c.denominator) for a, b, c in Upoly]

    def F(u, v):
        r0, r1 = mp.exp(u), mp.exp(v)
        num = mp.fsum(c * r0 ** a * r1 ** b for a, b, c in Ipm)
        Uv = mp.fsum(c * r0 ** a * r1 ** b for a, b, c in Upm)
        val = num / ((r0 + r1 + 1) * Uv ** k_U)
        if T:
            val /= mp.sqrt(Uv)
        return val * r0 ** (a0 + 1) * r1 ** (a1 + 1)
    val = mp.quad(F, [-mp.inf, 0, mp.inf], [-mp.inf, 0, mp.inf], maxdegree=maxdegree)
    return val, len(Ipoly), len(G.terms())


def sec_semianalytic(htk, rng, ntrial=2, Ts=(0, 1)):
    say("== theta graphs, T = %s: own density with exact Dirichlet integration over the path simplices"
        " and mpmath (tanh-sinh, 30 digits, maxdegree 8, log variables) over the two free path resistances,"
        " vs htk.kernel_exact" % (Ts,))
    spec = [((1, 2, 2), 1, 1), ((1, 3, 3), 2, 1), ((2, 2, 3), 2, 1), ((1, 3, 3), 3, 0), ((2, 2, 3), 3, 0)]
    for (pqr, H, T) in spec:
        if T not in Ts:
            continue
        G = htk.theta(*pqr)[0]
        net = net_of(G)
        paths = theta_paths(G, pqr)
        for trial in range(ntrial):
            lam = rand_lam(rng, net.g, H)
            C = net.cycle_basis(rng)
            t0 = time.time()
            Kh = htk.kernel_exact(G, lam, T)
            val, nI, nG = theta_semianalytic(net, lam, T, C, paths, check_rng=rng)
            import mpmath as mp
            if T:
                Khm = sum(mp.mpf(q.numerator) / q.denominator * mp.pi ** k for k, q in Kh.c.items())
            else:
                Khm = mp.mpf(Kh.numerator) / Kh.denominator
            rel = abs(val - Khm) / abs(Khm)
            dg = float(-mp.log10(rel)) if rel != 0 else 25.0
            say("   theta%s H=%d T=%d: htk %s = %s; own %s (G: %d terms, int G dx: %d terms) [%.0fs]"
                % (pqr, H, T, Kh, mp.nstr(Khm, 22), mp.nstr(val, 22), nG, nI, time.time() - t0))
            record(dg >= 12, "semianalytic theta%s H%dT%d trial %d" % (pqr, H, T, trial), "%.1f digits agree" % dg)


def sec_lattice6(htk, rng, hs=(1.2, 1.0, 0.9), tol=1e-8):
    say("== structure-free 6-dim log-lattice on the simplex (all N - 1 = 6 log Schwinger times, no path"
        " coordinates) for |I| = 7 thetas (low accuracy check)")
    spec = [((2, 2, 3), 3, 0), ((2, 2, 3), 2, 1)]
    for (pqr, H, T) in spec:
        G = htk.theta(*pqr)[0]
        net = net_of(G)
        lam = rand_lam(rng, net.g, H)
        C = net.fund_basis()
        Kh = htk.kernel_exact(G, lam, T)
        fd = FDens(net, lam, T, C)
        f = log_integrand(fd)
        out = []
        for h in hs:
            t0 = time.time()
            v, npts = flood_lattice_np(f, 6, h, tol, batch=200000, maxpts=1.2e8)
            out.append("h=%.2f: %.10g (%d pts, %.0fs)" % (h, v, npts, time.time() - t0))
            dg = digits(v, kval(Kh))
        say("   theta%s H=%d T=%d: htk %s = %.10g; %s" % (pqr, H, T, Kh, kval(Kh), "; ".join(out)))
        record(dg >= 5, "6-dim lattice theta%s H%dT%d" % (pqr, H, T), "%.1f digits agree at h=%.2f" % (dg, hs[-1]))


def src_all(lam, g):
    """sources at all vertices (root = minus the sum)."""
    H = len(lam[0])
    return [list(r) for r in lam] + [[-sum(lam[v][i] for v in range(g)) for i in range(H)]]


def sec_transform(htk, rng, ntrial=4):
    say("== transformation rules (D = H + T): reversal (-1)^{Dk} (z_e -> -z_e), edge order sgn^{D+1},"
        " vertex relabelling sgn^D")
    pool = [("triangle", htk.polygon(3), 2, 0), ("fan2", htk.fan(2), 2, 0), ("box", htk.box(), 3, 0),
            ("theta(1,3,3) H3T0", htk.theta(1, 3, 3)[0], 3, 0), ("theta(2,2,3) H3T0", htk.theta(2, 2, 3)[0], 3, 0),
            ("graph_A", htk.graph_A(), 3, 0),
            ("theta(1,2,2) H1T1", htk.theta(1, 2, 2)[0], 1, 1), ("theta(2,2,3) H2T1", htk.theta(2, 2, 3)[0], 2, 1),
            ("theta(1,3,3) H2T1", htk.theta(1, 3, 3)[0], 2, 1)]
    pool_pw = pool + [("prism H2T0", extra_graphs(htk)[0][1], 2, 0), ("prism H1T1", extra_graphs(htk)[0][1], 1, 1)]
    # pointwise (own densities; also htk.brute)
    for (label, G, H, T) in pool_pw:
        net = net_of(G)
        D_ = H + T
        bad = []
        for trial in range(ntrial):
            lam = rand_lam(rng, net.g, H)
            t = rand_t(rng, net.N)
            q, U = density_exact(net, t, lam, T, net.cycle_basis(rng), rng, full_check=False)
            # reversal of a random subset, with exp(-sum y.z)
            which = [e for e in range(net.N) if rng.random() < 0.5]
            Gr = G.reverse_edges(which)
            netr = net_of(Gr)
            qr, Ur = density_exact(netr, t, lam, T, netr.cycle_basis(rng), rng, full_check=False)
            s = (-1) ** (D_ * len(which))
            if not (qr == s * q and Ur == U):
                bad.append("reversal %s: %s vs %s" % (which, qr, s * q))
            z = [[rand_fr(rng) for _ in range(H)] for _ in range(net.N)]
            zr = [[-x for x in z[e]] if e in which else z[e] for e in range(net.N)]
            I1 = htk.brute.integrand(G, t, lam, T, z=z)
            I2 = htk.brute.integrand(Gr, t, lam, T, z=zr)
            if not (I2.coeff == s * I1.coeff and I2.expo == I1.expo and I2.radicand == I1.radicand):
                bad.append("htk integrand reversal")
            # edge permutation: new edge j = old edge order[j], t'_j = t_{order[j]}
            order = list(range(net.N))
            rng.shuffle(order)
            Gp = G.permute_edges(order)
            netp = net_of(Gp)
            tp = [t[e] for e in order]
            qp, Up = density_exact(netp, tp, lam, T, netp.cycle_basis(rng), rng, full_check=False)
            s = perm_sign(order) ** (D_ + 1)
            if not (qp == s * q and Up == U):
                bad.append("permutation %s: %s vs %s" % (order, qp, s * q))
            Dp = htk.brute.density(Gp, tp, lam, T)
            if Dp.coeff != qp:
                bad.append("htk permuted density")
            # vertex relabelling (root may move)
            nl = list(range(net.nv))
            rng.shuffle(nl)
            Gv = G.relabel(nl)
            netv = net_of(Gv)
            la = src_all(lam, net.g)
            lav = [None] * net.nv
            for v in range(net.nv):
                lav[nl[v]] = la[v]
            lamv = lav[:net.g]
            qv, Uv = density_exact(netv, t, lamv, T, netv.cycle_basis(rng), rng, full_check=False)
            s = perm_sign(nl) ** D_
            if not (qv == s * q and Uv == U):
                bad.append("relabel %s: %s vs %s" % (nl, qv, s * q))
            Dv = htk.brute.density(Gv, t, lav, T)
            if Dv.coeff != qv:
                bad.append("htk relabelled density")
        record(not bad, "pointwise transformation rules %s (%d trials)" % (label, ntrial), "; ".join(bad[:3]))
    # kernel level (htk.kernel_exact on transformed graphs)
    for (label, G, H, T) in pool:
        net = net_of(G)
        D_ = H + T
        bad = []
        Ks = []
        t0 = time.time()
        for trial in range(2):
            lam = rand_lam(rng, net.g, H)
            K = htk.kernel_exact(G, lam, T)
            Ks.append(K)
            if K == 0:
                bad.append("K = 0 (test vacuous)")
            which = [e for e in range(net.N) if rng.random() < 0.5]
            Kr = htk.kernel_exact(G.reverse_edges(which), lam, T)
            if Kr != K * (-1) ** (D_ * len(which)):
                bad.append("reversal %s: %s vs %s" % (which, Kr, K))
            order = list(range(net.N))
            rng.shuffle(order)
            Kp = htk.kernel_exact(G.permute_edges(order), lam, T)
            if Kp != K * perm_sign(order) ** (D_ + 1):
                bad.append("permutation %s: %s vs %s" % (order, Kp, K))
            nl = list(range(net.nv))
            rng.shuffle(nl)
            la = src_all(lam, net.g)
            lav = [None] * net.nv
            for v in range(net.nv):
                lav[nl[v]] = la[v]
            Kv = htk.kernel_exact(G.relabel(nl), lav, T)
            if Kv != K * perm_sign(nl) ** D_:
                bad.append("relabel %s: %s vs %s" % (nl, Kv, K))
            if T == 0 and trial == 0:
                # moments: reversing edge e with insertion (y^i_e)^p gives (-1)^p in addition
                e1 = rng.randrange(net.N)
                powers = {(e1, 0): 1, ((e1 + 1) % net.N, H - 1): 2}
                Km = htk.kernel_exact(G, lam, T, powers=powers)
                which = [e1]
                Kmr = htk.kernel_exact(G.reverse_edges(which), lam, T, powers=powers)
                sgn_extra = (-1) ** sum(p for (e, i), p in powers.items() if e in which)
                if Kmr != Km * (-1) ** (D_ * len(which)) * sgn_extra:
                    bad.append("reversal with insertion: %s vs %s" % (Kmr, Km))
        record(not bad, "kernel_exact transformation rules %s" % label,
               ("; ".join(bad[:3]) if bad else "K = %s, %s" % (Ks[0], Ks[1])) + " [%.1fs]" % (time.time() - t0))


def sec_n4(htk, rng, npts=4):
    say("== four-loop series-parallel graphs at (D, n) = (2, 4) (H, T) = (1, 1), (2, 0) and ladder(4) (H, T) = (2, 1):"
        " htk reduced density vs own D(t(xi)) |J| (pointwise; own Pfaffian with n = 4)")
    gl = []
    for G in htk.catalog.tight_graphs(2, 4):
        if htk.catalog.is_series_parallel(G):
            gl.append((G.name, G, 1, 1))
            gl.append((G.name, G, 2, 0))
    gl.append(("ladder(4)", htk.ladder(4)[0], 2, 1))
    for (label, G, H, T) in gl:
        net = net_of(G)
        t0 = time.time()
        red = htk.sp_reduce(G, prefer="large" if T else "small")
        lam = rand_lam(rng, net.g, H)
        dens = red.density(lam, T, check=0)
        tx, J, coords = ratio_jacobian(red)
        C = net.cycle_basis(rng)
        bad, own = [], []
        for k in range(npts):
            pt = red.random_point(rng)
            vals = tuple(pt[c] for c in red.coords)
            sub = {c: sp.Rational(v.numerator, v.denominator) for c, v in zip(coords, vals)}
            t = [Fr(int(sp.Rational(e.subs(sub)).p), int(sp.Rational(e.subs(sub)).q)) for e in tx]
            Jr = sp.Rational(J.subs(sub).det())
            Jv = Fr(int(Jr.p), int(Jr.q))
            q, U = density_exact(net, t, lam, T, C, rng, full_check=False)
            rat, rad = dens.evaluate(vals)
            gm = q * abs(Jv)
            own.append(gm)
            if T == 0:
                ok = rad == 1 and rat == gm
            else:
                ok = ((rat > 0) - (rat < 0) == (gm > 0) - (gm < 0)) and rat * rat * rad == gm * gm / U
            if not ok:
                bad.append("htk %s*sqrt(%s), own %s (U=%s)" % (rat, rad, gm, U))
        K = red.integrate(lam, T, dens=dens)
        record(not bad, "reduced density %s %s H%dT%d" % (label, sorted(G.degree(v) for v in range(G.nv)), H, T),
               ("; ".join(bad[:2]) if bad else "own g(xi) nonzero at %d/%d points; htk K = %s" %
                (sum(1 for x in own if x != 0), len(own), K)) + " [%.1fs]" % (time.time() - t0))


def sec_theta166(htk, rng):
    say("== theta(1,6,6), H = 5, T = 1: htk reduced density (prefer='large') vs own D(t(xi)) |J|")
    sec_reduced(htk, rng, which=["theta(1, 6, 6) H5T1"], npts=3, prefer="large")


def main(argv):
    ALL = ["self", "pointwise", "reduced", "theta166", "exact1", "integrals", "integrals6", "semianalytic0",
           "semianalytic1", "lattice6", "transform", "n4"]
    secs = argv[1:] or ["all"]
    if "all" in secs:
        secs = ALL
    rng = random.Random(20261120)
    htk = None if secs == ["self"] else import_htk()
    T0 = time.time()
    for s in secs:
        t0 = time.time()
        if s == "self":
            sec_self(rng)
        elif s == "pointwise":
            sec_pointwise(htk, rng)
        elif s == "reduced":
            sec_reduced(htk, rng)
        elif s == "exact1":
            sec_exact1(htk, rng)
        elif s == "integrals":
            sec_integrals(htk, rng)
        elif s == "integrals6":
            sec_integrals6(htk, rng)
        elif s == "lattice6":
            sec_lattice6(htk, rng)
        elif s == "transform":
            sec_transform(htk, rng)
        elif s == "semianalytic":
            sec_semianalytic(htk, rng)
        elif s == "semianalytic0":
            sec_semianalytic(htk, rng, Ts=(0,))
        elif s == "semianalytic1":
            sec_semianalytic(htk, rng, Ts=(1,))
        elif s == "n4":
            sec_n4(htk, rng)
        elif s == "theta166":
            sec_theta166(htk, rng)
        else:
            raise SystemExit("unknown section %s" % s)
        say("   [section %s: %.1fs]" % (s, time.time() - t0))
    say("discrepancies: %d" % len(DISC))
    for d in DISC:
        say("  " + d)
    say("total %.1fs" % (time.time() - T0))
    return 0 if not DISC else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
