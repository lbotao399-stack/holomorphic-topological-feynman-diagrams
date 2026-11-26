"""Pointwise exact checks of the three-loop kernels against the brute-force layer of q3lib.

For each graph:  currents (polynomial in r) = Kirchhoff currents at t(r);
density  s * prod(local factors)  =  eps_T det(B_T) det[d y^a_c / d r_k] sgn(omega_m(d_r))  at random
rational r, with random integer sources at all seven non-root vertices, for the cotree of the moves
and for two further spanning trees;  -sum_e y_e.z_e = sum_v lam_v.X_v at random rational z.
"""
import random
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations

import sympy as sp

from q3lib import (GRAPHS, GRAPHS_S2, ALL, LAM, Z, R, Reduction, Brute, fdet, random_point, lam_sub, outputs)

rng = random.Random(20261120)
NPTS = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def spanning_trees(edges, k, rng, ROOT, nloop):
    """k random spanning trees (as cotrees)."""
    N = len(edges)
    out = []
    tries = 0
    while len(out) < k and tries < 1000:
        tries += 1
        cot = sorted(rng.sample(range(N), nloop))
        tree = [e for e in range(N) if e not in cot]
        BT = [[0] * len(tree) for _ in range(ROOT)]
        for j, e in enumerate(tree):
            a, b = edges[e]
            if a < ROOT:
                BT[a][j] += 1
            if b < ROOT:
                BT[b][j] -= 1
        if fdet(BT) != 0 and cot not in out:
            out.append(cot)
    return out


results = {}
names = sys.argv[2:] if len(sys.argv) > 2 else list(GRAPHS_S2) + list(GRAPHS)
for name in names:
    t0 = time.time()
    red = Reduction(name)
    ROOT = red.root
    nloop = len(red.steps)
    nr = 3 * nloop
    RR = R[:nr]
    E = red.edges
    N = len(E)
    facs = red.factors()
    dom = red.domain()
    dt = [[sp.diff(te, r) for r in RR] for te in red.t]
    yfun = red.y
    X, expo = red.displacement()
    cot0 = red.cotree()
    trees = [cot0] + spanning_trees(E, 2, rng, ROOT, nloop)
    brutes = [Brute(E, c, ROOT) for c in trees]
    sgn = None
    nchk = 0
    for k in range(NPTS):
        lam = [[Fr(rng.randint(-6, 6)) for _ in range(3)] for _ in range(ROOT)]
        pt = random_point(dom, rng)
        subs = {r: sp.Rational(pt[r].numerator, pt[r].denominator) for r in RR}
        t = [Fr(sp.Rational(te.subs(subs)).p, sp.Rational(te.subs(subs)).q) for te in red.t]
        assert all(x > 0 for x in t)
        # inverse map: ratios from t (definitions of S2)
        J_t = [[sp.Rational(d.subs(subs)) for d in row] for row in dt]
        Jt = [[Fr(int(x.p), int(x.q)) for x in row] for row in J_t]      # N x 9: dt_e / dr_k
        omega = fdet([t] + [[Jt[e][k] for e in range(N)] for k in range(nr)])
        ls = lam_sub(lam)
        ypoly = [[sp.Rational(yfun[e][a].subs(ls).subs(subs)) for a in range(3)] for e in range(N)]
        ypoly = [[Fr(int(x.p), int(x.q)) for x in row] for row in ypoly]
        dens_poly = sp.Rational(sp.prod([f.subs(ls).subs(subs) for f in facs]))
        dens_poly = Fr(int(dens_poly.p), int(dens_poly.q))
        vals = []
        for B in brutes:
            y, D = B.dy_dt(t, lam)
            B.check_laws(t, lam, y)
            assert y == ypoly, (name, 'currents differ')
            rows = []
            for c in B.cot:
                for a in range(3):
                    rows.append([sum((D[c][a][j] * Jt[j][kk] for j in range(N)), Fr(0)) for kk in range(nr)])
            Jy = fdet(rows)
            vals.append(B.epsT * B.detBT * Jy * (1 if omega > 0 else -1))
        assert all(v == vals[0] for v in vals), (name, 'tree dependence', vals)
        if dens_poly != 0:
            s = vals[0] / dens_poly
            assert s in (1, -1), (name, s)
            if sgn is None:
                sgn = s
            assert s == sgn, (name, 'sign not constant')
        else:
            assert vals[0] == 0
        # displacement map at random z
        nz = len(red.spokes)
        zs = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(3)] for _ in range(nz)]
        zsub = {Z[j][a]: sp.Rational(zs[j][a].numerator, zs[j][a].denominator) for j in range(nz) for a in range(3)}
        lhs = Fr(0)
        for j, (e, s) in enumerate(red.spokes):
            lhs -= sum(ypoly[e][a] * s * zs[j][a] for a in range(3))
        rhs = 0
        for v in range(ROOT):
            for a in range(3):
                rhs += lam[v][a] * X[v][a].subs(subs).subs(zsub)
        rhs = sp.Rational(rhs)
        assert lhs == Fr(int(rhs.p), int(rhs.q)), (name, 'displacement')
        nchk += 1
    results[name] = sgn
    print('%s: %d points, currents = Kirchhoff, density = %+d * prod(factors) for 3 spanning trees '
          '(cotrees %s), displacement map OK; eps_T det B_T (move cotree %s) = %+d  [%.1fs]'
          % (name, nchk, sgn, trees, cot0, brutes[0].epsT * brutes[0].detBT, time.time() - t0), flush=True)

print('signs s_Gamma:', results)
