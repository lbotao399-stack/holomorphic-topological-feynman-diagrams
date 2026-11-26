"""V1 (claim a0) and V2 (claim a1): chart, U, Jacobian sign, grounded Green's function.
Own code (vlib); brute-force period matrix in the path-tree cycle basis (b_r - b_{n+1} - sum_{k>=r} a_k)."""
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

from vlib import (fan_cycles, t_of_r, r_of_t, rand_r, jac_det, jac_symbolic, period, fdet, laplacian,
                  finv, claim_gh, claim_R, ia, ib)

rng = random.Random(20261120)
NS = [int(x) for x in sys.argv[1:]] or [2, 4, 6]

for n in NS:
    t0 = time.time()
    # --- symbolic bijectivity of the chart (n <= 4)
    us = {k: sp.Symbol("u%d" % k, positive=True) for k in range(1, n + 1)}
    vs = {k: sp.Symbol("v%d" % k, positive=True) for k in range(1, n + 1)}
    te = t_of_r(n, us, vs)
    fac = [sp.factor(x) for x in te]
    print("n=%d: t_e(r) =" % n, fac)
    if n <= 4:
        uu, vv = r_of_t(n, te)
        ok1 = all(sp.cancel(sp.together(uu[k] - us[k])) == 0 and sp.cancel(sp.together(vv[k] - vs[k])) == 0 for k in range(1, n + 1))
        ts = [sp.Symbol("t%d" % e, positive=True) for e in range(2 * n + 1)]
        uu2, vv2 = r_of_t(n, ts)
        tt2 = t_of_r(n, uu2, vv2)
        ratios = [sp.cancel(sp.together(tt2[e] / ts[e])) for e in range(2 * n + 1)]
        ok2 = all(sp.cancel(x - ratios[0]) == 0 for x in ratios)
        print("   symbolic: r_of_t(t_of_r(r)) = r: %s;  t_of_r(r_of_t(t)) = c(t) t with one common factor c(t): %s" % (ok1, ok2))
    # --- random rational points
    G, C = fan_cycles(n, "pathtree")
    for p in range(3):
        u, v = rand_r(n, rng)
        t = t_of_r(n, u, v)
        assert all(x > 0 for x in t)
        uu, vv = r_of_t(n, t)
        assert uu == u and vv == v
        sJ = jac_det(n, u, v)
        Jf = Fr((-1) ** (n - 1))
        for k in range(1, n):
            Jf *= u[k + 1] ** (2 * k) / (v[k] ** 2 * (1 - v[k]) ** (2 * k + 1))
        assert sJ == Jf, (sJ, Jf)
        # U: brute det M (path-tree basis) vs claim prod u_{k+1}^k (1-v_k)^{-k} v_k^{-1}; also prod R_i / v_i
        U = fdet(period(G, C, t))
        Uc = Fr(1)
        for k in range(1, n):
            Uc *= u[k + 1] ** k * (1 - v[k]) ** (-k) / v[k]
        R = claim_R(n, u, v)
        Uc2 = Fr(1)
        for i in range(1, n):
            Uc2 *= R[i] / v[i]
        assert U == Uc == Uc2, (U, Uc, Uc2)
        # Green's function
        Lap = laplacian(G, t)
        Gm = finv(Lap)               # rows/cols = path vertices 1..n+1
        g, h, P, R, c = claim_gh(n, u, v)
        for a_ in range(1, n + 2):
            for b_ in range(a_, n + 2):
                assert Gm[a_ - 1][b_ - 1] == g[a_] * h[b_], (a_, b_)
        for i in range(1, n + 1):
            assert g[i] * h[i + 1] - g[i + 1] * h[i] == -t[ia(i)]
        # triangle entries
        assert Gm[n - 1][n - 1] == u[n] * (1 - u[n]) and Gm[n - 1][n] == u[n] * v[n] and Gm[n][n] == v[n] * (1 - v[n])
        print("   point %d: t>0, r_of_t(t_of_r(r)) = r, det[t;dt/dr] = (-1)^{n-1} prod_k u_{k+1}^{2k} v_k^-2 (1-v_k)^-(2k+1) (sign %+d), U = det M (path-tree basis) = claim,"
              " G = Lambda^-1 semiseparable g_min h_max, Wronskian = -a_i" % (p, 1 if sJ > 0 else -1))
    print("   n=%d done (%.1f s)" % (n, time.time() - t0))
