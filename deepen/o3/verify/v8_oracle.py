"""Numerical cross-check (htk.oracle, lattice rule, no reduction code) of two J_1-vertex kernels not checked by
the explorer: theta(1,3,3) v=1 (external + internal half-edge) and theta(2,2,3) v=5 (three internal half-edges),
local order A, against the exact spred values with the own insertion (vgraph.insertion)."""
import time, sys
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from vgraph import *
import htk

t0 = time.time()
lamnum = {1: (Fr(3, 2), Fr(-1)), 2: (Fr(1, 2), Fr(2)), 3: (Fr(-2), Fr(1, 3)), 4: (Fr(1), Fr(-3, 2)), 0: (Fr(0), Fr(0))}
cases = [((1, 3, 3), 1), ((2, 2, 3), 5)]
for pqr, v in cases:
    G, ext, lam = theta_data(pqr)
    loc = local_orders(G, ext, 'A')
    h1, h2 = loc[v][0], loc[v][1]
    sub = {PS[ext.index(u)][c]: sp.Rational(lamnum[u][c].numerator, lamnum[u][c].denominator) for u in ext for c in range(2)}
    exact = sp.nsimplify(kernel_over_pi(G, lam, insertion(G, loc, v)).subs(sub)) * sp.pi
    # oracle: express [q_h1 q_h2] through monomials in the currents; external q via Kirchhoff (numeric lambda)
    def qlin(h):
        """q = sum_e c_e y_e + const vector: returns ({e: coeff}, const)."""
        if h[0] == 'x':
            return ({}, (float(lamnum[h[1]][0]), float(lamnum[h[1]][1])))
        e, s = h[1], h[2]
        return ({e: (-1 if s == 0 else 1)}, (0.0, 0.0))
    (ca, ka), (cb, kb) = qlin(h1), qlin(h2)
    lamlist = [[float(lamnum[u][0]), float(lamnum[u][1])] for u in range(G.nv - 1)]
    total, err = 0.0, 0.0
    # [qa qb] = qa^1 qb^2 - qa^2 qb^1
    terms = []
    # [qa qb] = qa^1 qb^2 - qa^2 qb^1,  q = sum_e c_e y_e + k
    for (i, j, sgn) in [(0, 1, 1), (1, 0, -1)]:
        A = [(c, (e, i)) for e, c in ca.items()] + ([(ka[i], None)] if ka[i] else [])
        Bq = [(c, (e, j)) for e, c in cb.items()] + ([(kb[j], None)] if kb[j] else [])
        for xa, pa in A:
            for xb, pb in Bq:
                powers = {}
                for pp in (pa, pb):
                    if pp is not None:
                        powers[pp] = powers.get(pp, 0) + 1
                terms.append((sgn * xa * xb, powers))
    for coef, powers in terms:
        val, er = htk.oracle.kernel(G, lamlist, T=1, powers=powers or None)
        total += coef * val
        err += abs(coef) * er
    ex = float(exact)
    print("theta%s v=%d: exact %s = %.10f;  oracle %.10f +- %.1e;  agreeing digits %.1f  (%.0fs)"
          % (pqr, v, exact, ex, total, err, -np.log10(abs(total - ex) / abs(ex)), time.time() - t0), flush=True)
