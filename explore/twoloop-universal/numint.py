"""Direct numerical integration (double precision) of the brute-force density over the Schwinger simplex.

Integrand: D(t) from the Kirchhoff solution (Laplacian solve) and the forms (cycle matrix C3), at the
points t(xi) of the chart xi = (w_A, w_B, x^A_2.., x^B_2.., x^C_2..),
    t_{P_j} = x^P_j prod_{P' != P} w_{P'}     (a section of RP^{N-1}_{>0} -> chart; D |J| is scale free),
times |J| = |det[t; dt/dxi]| computed from the analytic derivatives of this section.  No closed form
and no reduction formula enters the integrand.  Cube coordinates: stick-breaking on each simplex,
u = sin^2(pi s / 2) on each cube coordinate, Gauss-Legendre in s (m nodes per direction).
"""
import itertools
import math
import sys

import numpy as np

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twoloop-universal")
from twoloop import theta_graph, Theta  # noqa: E402  (graph data only)


class Setup:
    def __init__(self, pqr):
        nv, edges = theta_graph(*pqr)
        self.nv, self.edges, self.pqr = nv, edges, pqr
        self.th = Theta(nv, edges)
        self.N = len(edges)
        g = nv - 1
        B = np.zeros((self.N, g))
        for e, (a, b) in enumerate(edges):
            if a < g:
                B[e, a] += 1
            if b < g:
                B[e, b] -= 1
        self.B = B
        self.C = np.array(self.th.C3(), dtype=float)
        self.paths = {P: [e for e, _ in self.th.paths[P]] for P in "ABC"}


def density_batch(st, t, lam, H, T):
    """D(t) (w.r.t. iota_E vol, edge order), t: (n, N), lam: (nv-1, H) sources at the non-root vertices."""
    n, N = t.shape
    B, C = st.B, st.C
    if H > 0:
        Lap = np.einsum("ea,ne,eb->nab", B, 1.0 / t, B)
        phi = np.linalg.solve(Lap, np.broadcast_to(lam, (n,) + lam.shape))
        y = np.einsum("ea,nai->nei", B, phi) / t[:, :, None]          # (n, N, H)
        Psi = (C[None, None, :, :] * y.transpose(0, 2, 1)[:, :, None, :]).reshape(n, 2 * H, N)
    M = np.einsum("re,ne,se->nrs", C, t, C)
    U = M[:, 0, 0] * M[:, 1, 1] - M[:, 0, 1] * M[:, 1, 0]
    e0 = N - 1
    cols = [e for e in range(N) if e != e0]
    sg = -1.0 if e0 % 2 else 1.0
    if T == 0:
        d = np.linalg.det(Psi[:, :, cols])
        return sg * d / (t[:, e0] * U ** H)
    Mi = np.empty_like(M)
    Mi[:, 0, 0] = M[:, 1, 1] / U
    Mi[:, 1, 1] = M[:, 0, 0] / U
    Mi[:, 0, 1] = -M[:, 0, 1] / U
    Mi[:, 1, 0] = -M[:, 1, 0] / U
    top = np.zeros(n)
    for i, e in enumerate(cols):
        for j in range(i + 1, len(cols)):
            f = cols[j]
            a12 = np.zeros(n)
            for a in range(2):
                for b in range(2):
                    k = C[0, e] * C[a, e] * C[b, f] * C[1, f] - C[0, f] * C[a, f] * C[b, e] * C[1, e]
                    if k:
                        a12 += Mi[:, a, b] * k
            a12 *= -0.5
            if not np.any(a12):
                continue
            rest = [x for x in cols if x not in (e, f)]
            dd = np.linalg.det(Psi[:, :, rest]) if H > 0 else np.ones(n)
            top += (-1.0 if (i + j - 1) % 2 else 1.0) * a12 * dd
    return sg * top / (t[:, e0] * U ** H * np.sqrt(U))


def stick(u, v):
    """stick-breaking: cube coords u (n, m) with v = 1 - u -> barycentric (n, m+1), Jacobian (n,)"""
    n, m = u.shape
    out = np.empty((n, m + 1))
    rem = np.ones(n)
    jac = np.ones(n)
    for k in range(m):
        out[:, k] = rem * u[:, k]
        jac *= rem
        rem = rem * v[:, k]
    out[:, m] = rem
    # jacobian of (b_0..b_{m-1}) wrt u: prod_{k} rem_k, rem_k = prod_{i<k} v_i
    return out, jac


def chart_eval(st, s_all, lam, H, T):
    """integrand times cube Jacobian at Gauss nodes s_all (n, dim) in (0,1)^dim."""
    th = st.th
    nP = {P: th.n[P] for P in "ABC"}
    u = np.sin(np.pi * s_all / 2) ** 2
    v = np.cos(np.pi * s_all / 2) ** 2
    du = (np.pi / 2) * np.sin(np.pi * s_all)
    jac_cube = np.prod(du, axis=1)
    off = 0
    w, jw = stick(u[:, off:off + 2], v[:, off:off + 2])
    off += 2
    xs, jx = {}, {}
    for P in "ABC":
        m = nP[P] - 1
        if m == 0:
            xs[P] = np.ones((u.shape[0], 1))
            jx[P] = np.ones(u.shape[0])
        else:
            xs[P], jx[P] = stick(u[:, off:off + m], v[:, off:off + m])
            off += m
    n = u.shape[0]
    N = st.N
    wd = {"A": w[:, 0], "B": w[:, 1], "C": w[:, 2]}
    prod_other = {"A": wd["B"] * wd["C"], "B": wd["A"] * wd["C"], "C": wd["A"] * wd["B"]}
    t = np.empty((n, N))
    for P in "ABC":
        for j, e in enumerate(st.paths[P]):
            t[:, e] = xs[P][:, j] * prod_other[P]
    # derivatives of t wrt xi = (w_A, w_B, x^A_2.., x^B_2.., x^C_2..)
    rows = [t]
    dA = {"A": -wd["B"], "B": wd["C"] - wd["A"], "C": wd["B"]}
    dB = {"A": wd["C"] - wd["B"], "B": -wd["A"], "C": wd["A"]}
    for dd in (dA, dB):
        r = np.empty((n, N))
        for P in "ABC":
            for j, e in enumerate(st.paths[P]):
                r[:, e] = xs[P][:, j] * dd[P]
        rows.append(r)
    for P in "ABC":
        es = st.paths[P]
        for j in range(1, len(es)):
            r = np.zeros((n, N))
            r[:, es[j]] = prod_other[P]
            r[:, es[0]] = -prod_other[P]
            rows.append(r)
    Jm = np.stack(rows, axis=1)
    J = np.linalg.det(Jm)
    D = density_batch(st, t, lam, H, T)
    jac = jac_cube * jw
    for P in "ABC":
        jac = jac * jx[P]
    return D * np.abs(J) * jac


def gl_rule(m):
    x, wt = np.polynomial.legendre.leggauss(m)
    return (x + 1) / 2, wt / 2


def integrate(pqr, lam, H, T, m, chunk=400000):
    st = Setup(pqr)
    dim = st.N - 1
    s1, w1 = gl_rule(m)
    lam = np.asarray(lam, dtype=float).reshape(st.nv - 1, H) if H else np.zeros((st.nv - 1, 0))
    total = 0.0
    # iterate over the first (dim - 4) coordinates, vectorize the last 4
    head = max(0, dim - 4)
    grids = np.meshgrid(*([s1] * (dim - head)), indexing="ij")
    wgrid = np.meshgrid(*([w1] * (dim - head)), indexing="ij")
    tail_s = np.stack([g.ravel() for g in grids], axis=1)
    tail_w = np.prod(np.stack([g.ravel() for g in wgrid], axis=1), axis=1)
    for idx in itertools.product(range(m), repeat=head):
        hs = np.array([s1[i] for i in idx])
        hw = float(np.prod([w1[i] for i in idx])) if head else 1.0
        S = np.concatenate([np.broadcast_to(hs, (tail_s.shape[0], head)), tail_s], axis=1)
        vals = chart_eval(st, S, lam, H, T)
        total += hw * float(np.sum(tail_w * vals))
    return total, m ** dim
