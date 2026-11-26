"""(c) Numerical integration of the brute-force (1,1) density (no reduction formula used in the integrand).

 N1  n = 2: htk oracle, lattice rule in log Schwinger coordinates (no chart, no reduction).
 N2  n = 2, 4: Gauss-Legendre product rule in the S1 chart r(t) with x = s^2 (3 - 2 s) in every
     variable (removes x^{k/2}, (1-x)^{k/2} endpoint behaviour) and u_n = (1 - v_n) w on the simplex;
     integrand = brute-force density D(t(r)) |det[t; dt/dr]| = S(t) |det[1; d log t/dr]|, S = D prod t_e
     from htk.oracle.scaled_density_batch (Kirchhoff + Pfaffian of -1/2 dM M^-1 dM, numpy).
     Brute-force values at a few nodes compared with the exact rational brute force of fanlib.
 N3  n = 4: htk oracle, scrambled Sobol QMC on the Schwinger simplex (no chart).
Usage: python3 check_numeric.py [m_list for n=4, default 5,6,7]
"""
import sys
import time
import math
import itertools
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import fanlib as F
from htk import fan, oracle


def exact_value(n, lam):
    """exact K from the closed forms: endpoint (4 pi l1 lN)^{n/2}/((n+1) n!^2); n = 2 general:
    pi l1 l3 / 3; n = 4 general: pi^2 l1 l5 Q / 180."""
    l = [float(x) for x in lam]
    if n == 2:
        return math.pi * l[0] * l[2] / 3
    if n == 4:
        l1, l2, l3, l4, l5 = l
        Q = (2*l1*l3 + 2*l1*l4 + l1*l5 + 4*l2*l3 + 4*l2*l4 + 2*l2*l5 + 12*l3**2 + 4*l3*l4 + 2*l3*l5)
        return math.pi ** 2 * l1 * l5 * Q / 180
    raise ValueError


def chart_functions(n):
    us = [sp.Symbol("u%d" % k) for k in range(1, n + 1)]
    vs = [sp.Symbol("v%d" % k) for k in range(1, n + 1)]
    u = {k: us[k - 1] for k in range(1, n + 1)}
    v = {k: vs[k - 1] for k in range(1, n + 1)}
    t = F.t_of_r(n, u, v)
    order = []
    for k in range(1, n + 1):
        order += [u[k], v[k]]
    logt = [sp.log(te) for te in t]
    rows = [[sp.simplify(sp.diff(lt, x)) for lt in logt] for x in order]
    ft = sp.lambdify(order, t, "numpy")
    fd = sp.lambdify(order, rows, "numpy")
    return ft, fd


def smooth(s):
    return s * s * (3 - 2 * s), 6 * s * (1 - s)


def gl_chart(n, lam, m, chunk=200000):
    G = fan(n)
    st = oracle._Setup(G)
    ft, fd = chart_functions(n)
    x, w = np.polynomial.legendre.leggauss(m)
    x = (x + 1) / 2
    w = w / 2
    X, Wx = smooth(x)
    W = w * Wx
    dim = 2 * n
    lamf = np.array([[float(l)] for l in lam])
    total = 0.0
    npts = m ** dim
    for a in range(0, npts, chunk):
        idx = np.stack(np.unravel_index(np.arange(a, min(a + chunk, npts)), (m,) * dim), axis=1)
        vals = X[idx]                          # columns: u1, v1, ..., u_{n-1}, v_{n-1}, w, v_n (chart order)
        wts = np.prod(W[idx], axis=1)
        args = [vals[:, j] for j in range(dim)]
        vn = args[dim - 1]
        wv = args[dim - 2]
        un = (1 - vn) * wv
        args[dim - 2] = un
        wts = wts * (1 - vn)
        T = np.stack([np.broadcast_to(np.asarray(c, dtype=float), (len(wts),)) for c in ft(*args)], axis=1)
        Drows = fd(*args)
        Mat = np.empty((len(wts), dim + 1, dim + 1))
        Mat[:, 0, :] = 1.0
        for j in range(dim):
            for e in range(dim + 1):
                Mat[:, j + 1, e] = np.broadcast_to(np.asarray(Drows[j][e], dtype=float), (len(wts),))
        J = np.abs(np.linalg.det(Mat))
        S = oracle.scaled_density_batch(G, T, lamf, 1, None, None, None, st)
        total += float(np.sum(wts * S * J))
    return total


def spot_check(n, lam, npts=3):
    """oracle brute-force scaled density vs exact rational brute force (fanlib) at rational chart points."""
    import random
    rng = random.Random(5)
    G = fan(n)
    for _ in range(npts):
        u, v = F.random_r(n, rng)
        t = F.t_of_r(n, u, v)
        Dc, U = F.density_brute(n, t, [[Fr(x) for x in lam]], T=1)
        ex = float(Dc) / math.sqrt(float(U)) * float(np.prod([float(x) for x in t]))
        nu = oracle.scaled_density_batch(G, np.array([[float(x) for x in t]]), np.array([[float(l)] for l in lam]), 1)[0]
        assert abs(ex - nu) <= 1e-9 * abs(ex), (ex, nu)


if __name__ == "__main__":
    mlist = [int(x) for x in sys.argv[1].split(",") if x] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else [5, 6, 7]
    cases = {2: [[1, 0, 1], [3, -2, 5]], 4: [[1, 0, 0, 0, 1], [2, -1, 3, 1, -2]]}
    # N1
    for lam in cases[2]:
        t0 = time.time()
        val, err = oracle.kernel(fan(2), [[x] for x in lam], T=1)
        ex = exact_value(2, lam)
        print("N1 n=2 lam=%s: lattice %.10f (est. err %.1e), exact %.10f, rel.diff %.1e, %.0f s" % (
            lam, val, err, ex, abs(val - ex) / abs(ex), time.time() - t0))
    # N2
    for n, ms in ((2, [8, 10, 12]), (4, mlist)):
        if not ms:
            continue
        for lam in cases[n]:
            spot_check(n, lam)
            ex = exact_value(n, lam)
            for m in ms:
                t0 = time.time()
                val = gl_chart(n, lam, m)
                print("N2 n=%d lam=%s m=%d (%d nodes): %.12f exact %.12f rel.diff %.2e (%.0f s)" % (
                    n, lam, m, m ** (2 * n), val, ex, abs(val - ex) / abs(ex), time.time() - t0))
    # N3
    if "--qmc" in sys.argv:
        from scipy.stats import qmc
        for lam in cases[4]:
            t0 = time.time()
            G = fan(4)
            st = oracle._Setup(G)
            dim = G.ne - 1
            lamf = np.array([[float(l)] for l in lam])
            vals = []
            dropped = 0
            for rep in range(8):
                uq = qmc.Sobol(dim, scramble=True, seed=100 + rep).random(2 ** 17)
                vq = uq * uq * (3.0 - 2.0 * uq)
                jac = np.prod(6.0 * uq * (1.0 - uq), axis=1)
                for k in range(dim):
                    jac = jac * (1.0 - vq[:, k]) ** (dim - 1 - k)
                tt = oracle._collapse(vq)
                Mf = np.einsum("re,be,se->brs", st.Cf, tt / tt.max(axis=1)[:, None], st.Cf)
                ok = np.linalg.cond(Mf) < 1e12
                dropped += int((~ok).sum())
                d = oracle.scaled_density_batch(G, tt[ok], lamf, 1, None, None, None, st) / np.prod(tt[ok], axis=1)
                vals.append(float(np.sum(jac[ok] * d)) / 2 ** 17)
            vals = np.array(vals)
            ex = exact_value(4, lam)
            print("N3 n=4 lam=%s: QMC (8 x 2^17 Sobol points on the Schwinger simplex, %d points with cond(M) > 1e12 dropped) %.6f +- %.1e, exact %.6f, rel.diff %.1e (%.0f s)" % (
                lam, dropped, vals.mean(), vals.std(ddof=1) / math.sqrt(8), ex, abs(vals.mean() - ex) / abs(ex), time.time() - t0))
