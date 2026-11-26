"""My d (t_closed.d_of, full mode) against run-1 cochain.d_CE (read-only import) on a non-invariant cochain
X = alt([p0p1][p2p3] R), R a random integer 4-tensor, structure constants of run-1 lie.sl(2), lie.sl(3).
Both outputs evaluated at 6 random 5-momentum configurations, all entries compared exactly."""
import sys, itertools
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import lie, cochain as cc
import vcoch as vc
from t_closed import d_of, alt_gen, poly_F

mom = cc.Momenta(6, 2)
P = mom.p
rng = np.random.default_rng(3)
for N in (2, 3):
    g1 = lie.sl(N)
    n = g1.dim
    R = rng.integers(-2, 3, (n,) * 4).astype(np.int64)
    expr = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) * (P[2][0] * P[3][1] - P[2][1] * P[3][0])
    X1 = cc.alternate(mom, expr, R, 4)
    D1 = cc.d_CE(mom, g1, X1, 4)
    ok = True
    nz = False
    for _ in range(6):
        p = [rng.integers(-9, 10, 2) for _ in range(4)]
        p.append(-(p[0] + p[1] + p[2] + p[3]))
        Q = [tuple(int(t) for t in x) for x in p]
        V1 = np.zeros((n,) * 5, dtype=object)
        for mono, T in D1.items():
            mv = 1
            idx = 0
            for i in range(4):
                for c_ in range(2):
                    mv *= Q[i][c_] ** mono[idx]
                    idx += 1
            V1 = V1 + mv * T.astype(object)
        V2 = d_of(g1.fu, lambda q: alt_gen(R, q, poly_F), Q, 'full')
        ok &= bool(np.array_equal(V1.astype(np.int64), V2))
        nz |= bool(V2.any())
    print("sl%d (run-1 basis): my d = run-1 d_CE on all 5-tensor entries at 6 configs: %s ; output nonzero: %s" % (N, ok, nz))
