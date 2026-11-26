"""Kernels of the Q_2 graphs at T = 1:
  5d (H = 2): graph (A) = theta(1,3,3), graph (B) = theta(2,2,3) with the data of S2 Section 7
              (htk.graph_A, htk.graph_B: vertex order, edge order, orientations, root O);
  7d (H = 3): theta(1,4,4), theta(2,3,4), theta(3,3,3) with O = N (root);
  3d (H = 1): theta(1,2,2) with O = N.
Checks: (i) pointwise density and currents against htk.brute (Kirchhoff at rational t, no reduction);
(ii) K(lam; 0) against the engine reduction htk.spred (symbolic lam); (iii) displacement maps
against brute currents; (iv) Taylor coefficients in z against spred.taylor; (v) numerical oracle.
"""
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import brute, oracle
from theta_kernels import ThetaKernel, to_fr

T00 = time.time()
DO_ORACLE = '--oracle' in sys.argv


def theta_rootN(p, q, r):
    """theta graph with N = root: interior vertices along A, B, C first, N' = nv - 2, N = nv - 1;
    edges oriented along the paths from N, edge order A, B, C."""
    nv = p + q + r - 1
    N, Np = nv - 1, nv - 2
    edges, paths, nxt = [], [], 0
    for n_ in (p, q, r):
        prev, pe = N, []
        for j in range(n_):
            cur = Np if j == n_ - 1 else nxt
            if j < n_ - 1:
                nxt += 1
            pe.append(len(edges))
            edges.append((prev, cur))
            prev = cur
        paths.append(pe)
    return htk.Graph(nv, edges, "thetaN(%d,%d,%d)" % (p, q, r)), paths


CASES = []
CASES.append(("5d (A) theta(1,3,3)", htk.graph_A(), [[3], [0, 1, 2], [6, 5, 4]], 2,
              {0: 'z1', 3: 'z2', 6: '-z3'}))
CASES.append(("5d (B) theta(2,2,3)", htk.graph_B(), [[0, 1], [2, 3], [6, 5, 4]], 2,
              {0: 'z1', 2: 'z2', 6: '-z3'}))
for pqr in [(1, 4, 4), (2, 3, 4), (3, 3, 3)]:
    G, paths = theta_rootN(*pqr)
    CASES.append(("7d theta%s" % (pqr,), G, paths, 3, {paths[0][0]: 'z1', paths[1][0]: 'z2', paths[2][0]: 'z3'}))
G, paths = theta_rootN(1, 2, 2)
CASES.append(("3d theta(1,2,2)", G, paths, 1, {paths[0][0]: 'z1', paths[1][0]: 'z2', paths[2][0]: 'z3'}))

RESULTS = {}
for name, G, paths, H, zspec in (CASES if __name__ == "__main__" else []):
    t0 = time.time()
    TK = ThetaKernel(G, paths, H, name)
    s = TK.check_pointwise(npts=6, seed=hash(name) % 1000)
    K0 = TK.K0()
    D = sp.factor(TK.Dtheta())
    assert set(K0) == {1}
    c = sp.simplify(K0[1] / TK.Dtheta())
    # (ii) engine, symbolic sources
    red = htk.sp_reduce(G, prefer="large")
    Ke = red.integrate(TK.lam, T=1, check=2)
    ok_engine = sp.expand(Ke - sp.pi * K0[1]) == 0
    # (iii) displacement against brute currents at random rational t
    zs = {k: [sp.Symbol('%s_%d' % (v.strip('-'), i)) for i in range(H)] for k, v in zspec.items()}
    zedge = {e: [(-1 if zspec[e].startswith('-') else 1) * c_ for c_ in zs[e]] for e in zspec}
    X = TK.displacement(zedge)
    rng = random.Random(5)
    ok_disp = True
    tex, _ = TK.t_of()
    for _ in range(3):
        wv, xv = TK.random_point(rng)
        lamn = [[Fr(rng.randint(-6, 6)) for _ in range(H)] for _ in range(G.nv - 1)]
        zn = {e: [Fr(rng.randint(-5, 5)) for _ in range(H)] for e in zedge}
        sub = {TK.w[0]: wv[0], TK.w[1]: wv[1], TK.w[2]: wv[2]}
        for P in range(3):
            for j in range(TK.nP[P]):
                sub[TK.x[P][j]] = xv[P][j]
        subs = {k: sp.Rational(v.numerator, v.denominator) for k, v in sub.items()}
        t = [Fr(1) * to_fr(sp.Rational(xv[P][j].numerator, xv[P][j].denominator) / sp.Rational(wv[P].numerator, wv[P].denominator))
             for P in range(3) for j in range(TK.nP[P])]
        tt = [None] * G.ne
        k = 0
        for P in range(3):
            for j, e in enumerate(TK.paths[P]):
                tt[e] = t[k]; k += 1
        yb = brute.currents(G, tt, lamn)
        lhs = -sum(yb[e][i] * (zn[e][i] * (-1 if zspec[e].startswith('-') else 1)) for e in zedge for i in range(H))
        zsub = {}
        for e in zs:
            for i in range(H):
                zsub[zs[e][i]] = sp.Rational(zn[e][i].numerator, zn[e][i].denominator)
        rhs = sum(sp.Rational(lamn[v][i].numerator, lamn[v][i].denominator) * X[v][i].subs(subs).subs(zsub)
                  for v in range(G.nv - 1) for i in range(H))
        ok_disp &= (sp.Rational(lhs.numerator, lhs.denominator) == sp.nsimplify(rhs))
    # (iv) Taylor coefficients in z (orders 1 and 2) against spred.taylor, numeric lam
    rng = random.Random(9)
    lamn = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(G.nv - 1)]
    lams = {TK.lam[v][i]: lamn[v][i] for v in range(G.nv - 1) for i in range(H)}
    E = sp.expand(sum(TK.lam[v][i] * X[v][i] for v in range(G.nv - 1) for i in range(H)).subs(lams))
    Dn = TK.Dtheta().subs(lams)
    mine = []
    for m in (1, 2):
        r = TK.integrate_poly(E ** m)
        mine.append(sp.expand(TK.s * Dn / 2 * sp.pi * r.get(1, 0) / sp.factorial(m)))
        assert set(r) <= {1}
    zlist = {e: zedge[e] for e in zedge}
    tay = red.taylor(lamn, zlist, 2, T=1, check=1)
    ok_tay = all(sp.expand(sp.sympify(tay[m].sympy() if hasattr(tay[m], 'sympy') else tay[m]) - mine[m - 1]) == 0
                 for m in (1, 2))
    RESULTS[name] = dict(s=s, c=c, D=D, TK=TK, X=X, zs=zs)
    print("%s: |I|=%d |V|=%d, paths n_P=%s, exponents alpha_P=%s" % (name, G.ne, G.nv, TK.nP, [str(a) for a in TK.alpha]))
    print("   pointwise density (6 points, sign and square) and currents vs brute: OK, s = %+d" % s)
    print("   D_Gamma = %s" % D)
    print("   K(lam;0) = %s * pi * D_Gamma;  engine (spred, symbolic lam) equal: %s" % (c, ok_engine))
    print("   displacement -sum y.z = sum lam.X at 3 rational points: %s;  Taylor orders 1,2 = spred.taylor: %s  (%.1fs)"
          % (ok_disp, ok_tay, time.time() - t0))
    if DO_ORACLE:
        t1 = time.time()
        lamf = [[float(v) for v in row] for row in lamn]
        val, err = htk.kernel_numeric(G, lamn, T=1)
        exact = float(c * Dn) * float(sp.pi)
        import math
        dig = -math.log10(abs(val - exact) / abs(exact)) if val != exact else 16
        print("   oracle: K(lam;0) = %.10g +- %.2g, exact %.10g, agreeing digits %.1f  (%.0fs)"
              % (val, err, exact, dig, time.time() - t1))
if __name__ == "__main__":
    print("total %.1fs" % (time.time() - T00))
