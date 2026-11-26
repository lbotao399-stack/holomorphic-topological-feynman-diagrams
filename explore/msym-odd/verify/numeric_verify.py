"""Direct numerical integration of the T = 1 two-loop Schwinger integral from the definition
(no reduction coordinates):  K = int_{R^{N-1}} S(e^s) d^{N-1}s,  S(t) = D(t) prod_e t_e,  t_{N-1} = 1,
D(t) from a vectorized Kirchhoff solve, M = C t C^T, rho_2 = U^{-1/2}(-1/2 dM M^-1 dM)_{12}, psi^i_r,
and the coefficient of dt_0..dt_{N-2}.  Scrambled Sobol points mapped by s = b logit(u); R independent
scrambles give mean and standard error.  Compared with the exact values (Dirichlet integrals of the
pointwise-verified density): K(lam;0) and the insertion int D(t) y^0_e (e = spoke of path B).
Usage: python3 numeric_verify.py case m R b
"""
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations

import numpy as np
from scipy.stats import qmc
import sympy as sp

from vkern import Theta, theta_rootN, S2_A, S2_B, dirichlet

case = sys.argv[1]
m = int(sys.argv[2]) if len(sys.argv) > 2 else 20
R = int(sys.argv[3]) if len(sys.argv) > 3 else 8
bsc = float(sys.argv[4]) if len(sys.argv) > 4 else 2.5
Lcut = float(sys.argv[5]) if len(sys.argv) > 5 else 18.0     # integration box |s_e| <= Lcut (tail ~ exp(-Lcut/2))
DATA = {'A': (S2_A, 2), 'B': (S2_B, 2), '122': (theta_rootN(1, 2, 2), 1),
        '144': (theta_rootN(1, 4, 4), 3), '234': (theta_rootN(2, 3, 4), 3), '333': (theta_rootN(3, 3, 3), 3)}
(nv, edges, paths), H = DATA[case]
TH = Theta(nv, edges, paths, H)
net = TH.net
N, g = net.N, nv - 1
rng = np.random.default_rng(7)
lam_int = [[int(v) for v in rng.integers(-3, 4, size=H)] for _ in range(g)]
if TH.Np != net.root:
    lam_int[TH.Np] = [0] * H
lam = np.array(lam_int, dtype=float)
B = np.array(net.B, dtype=float)          # N x g
C = np.array(net.C, dtype=float)          # 2 x N


def density(t):
    """t: P x N  ->  D(t), P-vector."""
    P = t.shape[0]
    L = np.einsum('eu,ev,pe->puv', B, B, 1.0 / t)
    phi = np.linalg.solve(L, np.broadcast_to(lam, (P, g, H)))
    y = np.einsum('ev,pvi->pei', B, phi) / t[:, :, None]
    M = np.einsum('re,se,pe->prs', C, C, t)
    U = M[:, 0, 0] * M[:, 1, 1] - M[:, 0, 1] * M[:, 1, 0]
    Mi = np.empty_like(M)
    Mi[:, 0, 0] = M[:, 1, 1] / U; Mi[:, 1, 1] = M[:, 0, 0] / U
    Mi[:, 0, 1] = -M[:, 0, 1] / U; Mi[:, 1, 0] = -M[:, 1, 0] / U
    u = np.einsum('re,se->rse', C, C)        # dM_rs = sum_e u_rse dt_e
    a = np.einsum('re,pei->pire', C, y).reshape(P, 2 * H, N)   # rows k = (i, r), i outer
    top = np.zeros(P)
    idx = list(range(N - 1))
    for e, f in combinations(idx, 2):
        r2 = np.zeros(P)
        for aa in range(2):
            for bb in range(2):
                r2 += Mi[:, aa, bb] * (u[0, aa, e] * u[bb, 1, f] - u[0, aa, f] * u[bb, 1, e])
        r2 *= -0.5
        rest = [x for x in idx if x != e and x != f]
        sgn = (-1) ** (e + f - 1)
        top += sgn * r2 * np.linalg.det(a[:, :, rest])
    return top / ((-1) ** (N - 1) * t[:, N - 1]) / U ** H / np.sqrt(U), y


# consistency of the vectorized density with the exact one at a rational point
rr = np.random.default_rng(1)
tq = [Fr(int(rr.integers(1, 20)), int(rr.integers(1, 9))) for _ in range(N)]
qD, Uq, _ = net.density(tq, [[Fr(c) for c in v] for v in lam_int])
Dex = float(qD) / float(Uq) ** 0.5
Dnum, _ = density(np.array([[float(v) for v in tq]]))
print("case %s: H=%d N=%d lam=%s; vectorized vs exact density at a rational point: %.15g vs %.15g"
      % (case, H, N, lam_int, Dnum[0], Dex))

# exact values: K0 = s (D/2) pi c0;  insertion int D y^0_{e*}, e* = first edge of path B
s = TH.check(npts=3, seed=3)
assert len(s) == 1
s = s.pop()
Dl = TH.Dtheta([[Fr(c) for c in v] for v in lam_int])
q, k = dirichlet(TH.alpha)
c0 = q
for n in TH.nP:
    c0 *= dirichlet([0] * n)[0]
K0 = float(Fr(s, 2) * Dl * c0) * np.pi
# insertion via the verified X rule: -sg * y^0_{e*} = sum_v lam^0_v X^0_v with z_B = unit vector e_0
w = sp.symbols('wA wB wC', positive=True)
xs = [sp.symbols('x%s1:%d' % ("ABC"[P], TH.nP[P] + 1), positive=True) for P in range(3)]
zP = [[sp.Integer(0)] * H for _ in range(3)]
zP[1][0] = sp.Integer(1)
Xs = TH.X_rule(list(w), [list(v) for v in xs], zP)
poly = sp.expand(sum(lam_int[v][0] * Xs[v][0] for v in Xs))
gens = list(w) + [v for P in range(3) for v in xs[P]]
Pl = sp.Poly(poly, *gens)
tot = Fr(0)
for mon, cf in Pl.terms():
    qq, kk = dirichlet([TH.alpha[P] + mon[P] for P in range(3)])
    assert kk == 2
    off = 3
    for P in range(3):
        qq *= dirichlet(list(mon[off:off + TH.nP[P]]))[0]
        off += TH.nP[P]
    tot += qq * Fr(int(sp.Rational(cf).p), int(sp.Rational(cf).q))
sgB = TH.sg[1][0]
Ins = -sgB * float(Fr(s, 2) * Dl * tot) * np.pi
eB = paths[1][0]

if __name__ == '__main__':
    # QMC
    d = N - 1
    t0 = time.time()
    vals0, vals1 = [], []
    nbad = 0
    for rep in range(R):
        sob = qmc.Sobol(d, scramble=True, seed=1000 + rep)
        acc0 = acc1 = 0.0
        nb = 0
        for chunk in range(1 << max(0, m - 18)):
            u = sob.random(1 << min(m, 18))
            u = np.clip(u, 1e-15, 1 - 1e-15)
            sv = bsc * np.log(u / (1 - u))
            logp = np.sum(-np.log(bsc) + np.log(u) + np.log(1 - u), axis=1)
            inside = np.all(np.abs(sv) <= Lcut, axis=1)
            nb += len(sv)
            sv, logp = sv[inside], logp[inside]
            t = np.concatenate([np.exp(sv), np.ones((sv.shape[0], 1))], axis=1)
            with np.errstate(all='ignore'):
                Dv, y = density(t)
                Sv = Dv * np.prod(t, axis=1)
                wgt = Sv / np.exp(logp)
            bad = ~np.isfinite(wgt) | ~np.isfinite(y[:, eB, 0])
            nbad += int(bad.sum())
            wgt[bad] = 0.0
            yy = np.where(bad, 0.0, y[:, eB, 0])
            acc0 += wgt.sum(); acc1 += (wgt * yy).sum()
        vals0.append(acc0 / nb); vals1.append(acc1 / nb)
    vals0, vals1 = np.array(vals0), np.array(vals1)
    for nm, v, ex in [("K(lam;0)", vals0, K0), ("int D y^0_{e=%d}" % eB, vals1, Ins)]:
        mean, err = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
        rel = abs(mean - ex) / abs(ex)
        print("   %-18s QMC %.8g +- %.2g (2^%d points x %d scrambles, b=%.1f)   exact %.10g   rel. diff %.1e (%.1f digits), |diff|/err = %.1f"
              % (nm, mean, err, m, R, bsc, ex, rel, -np.log10(rel), abs(mean - ex) / err))
    print("   box |s| <= %.0f, non-finite evaluations set to 0: %d;  time %.0fs" % (Lcut, nbad, time.time() - t0))
