"""b3/b4: K33 [(0,1),(0,3),(0,5),(1,2),(1,4),(2,3),(2,5),(3,4),(4,5)], brute density (vk) in a chart built here.
Chart: K4 {1,2,3,5}: R12,R23,R25,R13,R15,R35 (R35 = 1), fractions f13, f15, f35 in (0,1):
  t12=R12, t23=R23, t25=R25;  r_ab = R_ab/f_ab (triangle of the star at 4), r'_ab = R_ab/(1-f_ab) (star at 0), ab in {13,15,35};
  Delta->Y: t01 = r'13 r'15/S', t03 = r'13 r'35/S', t05 = r'15 r'35/S';  t14 = r13 r15/S, t34 = r13 r35/S, t45 = r15 r35/S.
K = int D(t(xi)) |det[t; d_xi t]| dxi; log-lattice in log R12..log R15, quadrature in the fractions with
f = phi(u), phi(u) = u^p/(u^p + (1-u)^p) (endpoint clustering), Gauss-Legendre in u."""
import sys, time, itertools
import numpy as np
import sympy as sp
from vk import G, NG

K33 = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5)]


def chart():
    R12, R23, R25, R13, R15, f13, f15, f35 = sp.symbols('R12 R23 R25 R13 R15 f13 f15 f35', positive=True)
    R35 = sp.Integer(1)
    r13, r15, r35 = R13 / f13, R15 / f15, R35 / f35
    q13, q15, q35 = R13 / (1 - f13), R15 / (1 - f15), R35 / (1 - f35)
    Sq = q13 + q15 + q35; Sr = r13 + r15 + r35
    t = [q13 * q15 / Sq, q13 * q35 / Sq, q15 * q35 / Sq, R12, r13 * r15 / Sr, R23, R25, r13 * r35 / Sr, r15 * r35 / Sr]
    xi = [R12, R23, R25, R13, R15, f13, f15, f35]
    J = [[sp.diff(ti, x) for ti in t] for x in xi]
    return xi, t, J


class K33Int:
    def __init__(self, lams, m, p=1):
        xi, t, J = chart()
        self.tf = sp.lambdify(xi, t, 'numpy'); self.Jf = sp.lambdify(xi, J, 'numpy')
        self.ng = NG(G(6, K33), lams)
        u, w = np.polynomial.legendre.leggauss(m); u = (u + 1) / 2; w = w / 2
        if p == 1:
            f, df = u, np.ones_like(u)
        else:
            A, B = u ** p, (1 - u) ** p
            f = A / (A + B)
            df = p * (u ** (p - 1) * (1 - u) ** (p - 1)) / (A + B) ** 2
        self.F = np.array(list(itertools.product(f, repeat=3)))
        self.W = np.prod(np.array(list(itertools.product(w * df, repeat=3))), 1)

    def __call__(self, Ex):
        P = len(Ex); Q = len(self.W)
        a = [np.repeat(Ex[:, k], Q) for k in range(5)] + [np.tile(self.F[:, k], P) for k in range(3)]
        t = np.stack([np.broadcast_to(np.asarray(v, float), (P * Q,)) for v in self.tf(*a)], 1)
        J = np.stack([np.stack([np.broadcast_to(np.asarray(v, float), (P * Q,)) for v in row], 1) for row in self.Jf(*a)], 1)
        jac = np.abs(np.linalg.det(np.concatenate([t[:, None, :], J], 1)))
        D = self.ng.D(t)
        return np.einsum('pqs,q->ps', (D * jac[:, None]).reshape(P, Q, -1), self.W) * np.prod(Ex, 1)[:, None]
