"""b1: prism, direct integration of the brute density (vk) in a chart built here.
Chart (composition of bijections of positive orthants):
  K4 {0,1,2,3}: R01,R02,R03,R12,R13,R23 (R03 = 1);  fractions f13, f23, f15 in (0,1).
  t01=R01, t02=R02, t03=R03;  t23 = R23/f23, r23'' = R23/(1-f23);  r13 = R13/f13, r13'' = R13/(1-f13);  r12'' = R12;
  Delta->Y on {1,2,3} (r12'', r23'', r13''), centre 5:  S2 = r12''+r23''+r13'',
     rho1 = r12'' r13''/S2 (5-1), rho2 = r12'' r23''/S2 (5-2), rho3 = r13'' r23''/S2 (5-3);
  t25 = rho2;  t15 = rho1/f15, r15' = rho1/(1-f15);  r35 = rho3;
  Delta->Y on {1,3,5} (r13, r15', r35), centre 4:  S1 = r13+r15'+r35,
     t14 = r13 r15'/S1, t34 = r13 r35/S1, t45 = r15' r35/S1.
K = int D(t(xi)) |det[t; d_xi t]| dxi, xi = (log R01, log R02, log R12, log R13, log R23, f13, f23, f15):
log-lattice (trapezoid) in the 5 log R, Gauss-Legendre (m points) in each fraction.
Sources lam0 = e1, lam1 = e2, lam2 = (SA, SB), lam3 = -(sum), lam4 = lam5 = 0;
19 configurations (SA, SB): the principal lattice -1 + (i,j)/2, i+j <= 4, plus (1,1), (1/2,1), (1,-1/2), (0,0); K(SA,SB) fitted with a bivariate quartic;
c_{4-j-k,j,k} = (-1)^k [SA^k SB^j] K  (<l0,l1> = 1, <l0,l2> = SB, <l1,l2> = -SA)."""
import sys, time, itertools
import numpy as np
import sympy as sp
from vk import G, NG

PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]


def chart():
    R01, R02, R12, R13, R23, f13, f23, f15 = sp.symbols('R01 R02 R12 R13 R23 f13 f23 f15', positive=True)
    R03 = sp.Integer(1)
    t01, t02, t03 = R01, R02, R03
    t23 = R23 / f23; r23b = R23 / (1 - f23)
    r13 = R13 / f13; r13b = R13 / (1 - f13)
    r12b = R12
    S2 = r12b + r23b + r13b
    rho1 = r12b * r13b / S2; rho2 = r12b * r23b / S2; rho3 = r13b * r23b / S2
    t25 = rho2
    t15 = rho1 / f15; r15p = rho1 / (1 - f15)
    r35 = rho3
    S1 = r13 + r15p + r35
    t14 = r13 * r15p / S1; t34 = r13 * r35 / S1; t45 = r15p * r35 / S1
    t = [t01, t02, t03, t14, t15, t23, t25, t34, t45]
    xi = [R01, R02, R12, R13, R23, f13, f23, f15]
    J = [[sp.diff(ti, x) for ti in t] for x in xi]
    return xi, t, J


class PrismInt:
    def __init__(self, lams, m):
        xi, t, J = chart()
        self.tf = sp.lambdify(xi, t, 'numpy')
        self.Jf = sp.lambdify(xi, J, 'numpy')
        self.ng = NG(G(6, PRISM), lams)
        x, w = np.polynomial.legendre.leggauss(m)
        x = (x + 1) / 2; w = w / 2
        F = np.array(list(itertools.product(x, repeat=3))); W = np.prod(np.array(list(itertools.product(w, repeat=3))), 1)
        self.F, self.W = F, W

    def __call__(self, Ex):
        """Ex: (P,5) = (R01,R02,R12,R13,R23) -> (P,S): sum_f W D |det[t; dt]| * prod R (log measure)."""
        P = len(Ex); Q = len(self.W)
        a = [np.repeat(Ex[:, k], Q) for k in range(5)] + [np.tile(self.F[:, k], P) for k in range(3)]
        t = np.stack([np.broadcast_to(np.asarray(v, dtype=float), (P * Q,)) for v in self.tf(*a)], 1)
        Jl = self.Jf(*a)
        J = np.stack([np.stack([np.broadcast_to(np.asarray(v, dtype=float), (P * Q,)) for v in row], 1) for row in Jl], 1)  # (PQ, 8, 9)
        M = np.concatenate([t[:, None, :], J], 1)
        jac = np.abs(np.linalg.det(M))
        D = self.ng.D(t)
        val = (D * jac[:, None]).reshape(P, Q, -1)
        return np.einsum('pqs,q->ps', val, self.W) * np.prod(Ex, 1)[:, None]


def configs():
    grid = [(-1 + i / 2, -1 + j / 2) for i in range(5) for j in range(5) if i + j <= 4] + [(1, 1), (0.5, 1), (1, -0.5), (0, 0)]
    lams = [[[1, 0], [0, 1], [a, b], [-1 - a, -1 - b], [0, 0]] for (a, b) in grid]
    return grid, lams


def fit(grid, vals):
    mons = [(j, k) for j in range(5) for k in range(5) if j + k <= 4]          # SB^j SA^k
    A = np.array([[b ** j * a ** k for (j, k) in mons] for (a, b) in grid], dtype=float)
    co, *_ = np.linalg.lstsq(A, vals, rcond=None)
    return {(4 - j - k, j, k): (-1) ** k * c for (j, k), c in zip(mons, co)}, np.max(np.abs(A @ co - vals))


if __name__ == "__main__":
    import mpmath as mp
    from lattice2 import lattice
    h = float(sys.argv[1]); tol = float(sys.argv[2]); m = int(sys.argv[3])
    z3 = float(mp.zeta(3))
    CL = {(4, 0, 0): (-1 / 34560, 0), (3, 1, 0): (79 / 134400, -3 / 5600), (3, 0, 1): (1667 / 1209600, -3 / 2800),
          (2, 2, 0): (-79 / 134400, 3 / 5600), (2, 1, 1): (19 / 16800, -1 / 1400), (2, 0, 2): (-193 / 403200, 1 / 2800),
          (1, 3, 0): (-1 / 11520, 0), (1, 2, 1): (41 / 44800, -3 / 2800), (1, 1, 2): (-29 / 134400, -1 / 2800),
          (1, 0, 3): (299 / 1209600, -1 / 2800)}
    grid, lams = configs()
    I = PrismInt(lams, m)
    t0 = time.time()
    v0, b0, n0 = lattice(I, 5, h, tol, 0.0)
    v1, b1, n1 = lattice(I, 5, h, tol, h / 2)
    print("prism h=%.2f tol=%.0e m=%d  points %d  %.0fs  boundary-layer mass (max over configs) %.1e" % (h, tol, m, n0 + n1, time.time() - t0, np.max(b0)))
    c0, r0 = fit(grid, v0); c1, r1 = fit(grid, v1)
    print("fit residuals %.1e %.1e" % (r0, r1))
    for key in sorted(c0, key=lambda k: (-k[0], -k[1])):
        q, z = CL.get(key, (0, 0)); cl = q + z * z3
        print("  c_%d%d%d  shift0 %+.10e  shift h/2 %+.10e  claim %+.10e  diff %.1e" % (key + (c0[key], c1[key], cl, c0[key] - cl)))
    k = grid.index((0, 0))
    print("K at lam2 = 0: %.10e (claim c_400 = %.10e)" % (v0[k], -1 / 34560))
