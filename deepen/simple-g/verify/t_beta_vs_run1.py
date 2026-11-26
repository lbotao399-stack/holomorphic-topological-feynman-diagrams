"""beta (vcoch.beta_tensor, derived formula) against run-1 cochain.d_CE(alt([q0 q1]^2 f)) (read-only import),
run-1 lie.sl(2), lie.sl(3) structure constants and trace form; expected ratio 6 (alt over S3 of an S3-invariant
eta_0 gives 6 eta_0).  Values at 4 random configurations, all tensor entries."""
import sys
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import lie, cochain as cc
import vcoch as vc

mom = cc.Momenta(6, 2)
P = mom.p
rng = np.random.default_rng(8)
for N in (2, 3):
    g1 = lie.sl(N)
    n = g1.dim
    expr = (P[0][0] * P[1][1] - P[0][1] * P[1][0]) ** 2
    eta = cc.alternate(mom, expr, g1.fl, 3)
    D = cc.d_CE(mom, g1, eta, 3)
    TT = np.einsum('abe,cde->abcd', g1.fu, g1.fl)
    ratios = set()
    for _ in range(4):
        Q = vc.config4(rng)
        V = np.zeros((n,) * 4, dtype=object)
        for mono, T in D.items():
            mv = 1
            idx = 0
            for i in range(3):
                for c_ in range(2):
                    mv *= Q[i][c_] ** mono[idx]
                    idx += 1
            V = V + mv * T.astype(object)
        B = vc.beta_tensor(TT, Q).astype(object)
        nz = np.flatnonzero(B.reshape(-1) != 0)
        r = Fr(int(V.reshape(-1)[nz[0]]), int(B.reshape(-1)[nz[0]]))
        ok = bool(np.array_equal(V * r.denominator, B * r.numerator))
        ratios.add((str(r), ok))
    print("sl%d (run-1 basis, tr form): d_CE(alt([q0q1]^2 f)) / beta at 4 configs: %s" % (N, sorted(ratios)))
