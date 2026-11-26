"""Cross-check of the index forms (claim c3) used here against the run-1 obstruction cochains
omega_G = alt((60/pi) K_G W_G^{orientation rule}) (run-1 code, read-only import: lie, cochain, c_obstruction, htk).
Claim to test: omega_133 = -3 F_133, omega_223 = 2 F_223 as functions (exact rationals, random points)."""
import sys
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import lie, cochain as cc
from c_obstruction import kernel_poly, mom
import valg, vcoch as vc
import sympy as sp

rng = np.random.default_rng(11)
P = mom.p
for N in (2, 3):
    g1 = lie.sl(N)
    om = {}
    for pqr in [(1, 3, 3), (2, 2, 3)]:
        G, ext, K = kernel_poly(pqr)
        W, sc = cc.colour_weight(g1, G.nv, G.edges, ext)
        om[pqr] = (cc.alternate(mom, K, W, 4), sc)
        print("sl%d theta%s: (60/pi) K = %s ; ext = %s" % (N, pqr, sp.factor(K), ext))
    A = valg.sl(N)
    S = vc.Setup(A)
    coordsA = valg.coords_solver([m // 1 for m in A.mats])
    res = {'133': [], '223': []}
    for trial in range(3):
        mats = []
        for _ in range(4):
            m = rng.integers(-3, 4, (N, N))
            m[N - 1, N - 1] = -sum(m[i, i] for i in range(N - 1))
            mats.append(m.astype(np.int64))
        Pv = vc.config4(rng)
        X1 = [g1.coords(m) for m in mats]                  # run-1 coordinates
        XA = [[int(t) for t in coordsA(m)] for m in mats]  # my coordinates (integral for this basis)
        for pqr, G in [((1, 3, 3), '133'), ((2, 2, 3), '223')]:
            cochain, sc = om[pqr]
            val = Fr(0)
            for mono, T in cochain.items():
                mv = 1
                idx = 0
                for i in range(3):
                    for c_ in range(2):
                        mv *= Pv[i][c_] ** mono[idx]
                        idx += 1
                tv = int(np.einsum('abcd,a,b,c,d->', T.astype(object), *[np.array(x, dtype=object) for x in X1]))
                val += mv * tv
            val = val / Fr(sc)
            F = Fr(vc.alt_point(lambda Y: S.W_point(G, Y), XA, Pv), S.s ** 7)
            res[G].append((val, F, val / F if F else None))
    print("sl%d: kappa(run-1) = tr, kappa(here) = K/%d (= tr for sl_N iff g = 2N: %s)" % (N, S.g, S.g == 2 * N))
    for G in ('133', '223'):
        print("   omega_%s / F_%s at 3 random points: %s" % (G, G, [str(r[2]) for r in res[G]]))
