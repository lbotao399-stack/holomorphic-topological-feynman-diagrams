"""V13/V14 (claims c1, c2, b1 order 1 at n = 2): numerical integration of the brute-force density
(own vectorized float implementation of the definition: period matrix in the path-tree basis,
A_rs = -1/2 sum_ab dM_ra M^-1_ab dM_bs written as A_rs[e,f] = -1/2 Q_ef (C_re C_sf - C_rf C_se),
Q = C^T M^-1 C, nodal Kirchhoff solve, Pfaffian over matchings, top coefficient by Laplace pairing).

  lattice:  K = int_{R^{N-1}} D(e^s) prod_e e^{s_e} ds  (s_{N-1} = 0), s_e = alpha sinh(x_e),
            trapezoid in x  (chart-free; n = 2)
  chart:    K = int g(r) dr on (0,1)^{2n-2} x Delta^2, r_j = sin^2 theta_j, v_n = (1-u_n) sin^2 phi,
            Gauss-Legendre in (theta, phi); g = D(t(r)) |det[t; dt/dr]| with
            |det| = prod_k u_{k+1}^{2k} v_k^-2 (1-v_k)^-(2k+1) (v1_chart_green / jac_sym)  (n = 4)
"""
import itertools
import math
import sys
import time
from fractions import Fraction as Fr

import numpy as np

from vlib import fan_cycles, brute_density, t_of_r, rand_r, perm_sign, ia, ib
import random


class FloatFan:
    def __init__(self, n):
        self.n = n
        G, C = fan_cycles(n, "pathtree")
        self.G = G
        self.C = np.array(C, dtype=float)                   # n x N
        self.N = G.N
        self.B = np.array(G.incidence(), dtype=float)       # N x V
        m = 2 * n                                           # variables after dropping dt_{N-1}
        self.subsets = list(itertools.combinations(range(m), n))
        self.comp = [tuple(sorted(set(range(m)) - set(S))) for S in self.subsets]
        self.ssign = np.array([perm_sign(list(S) + list(Sc)) for S, Sc in zip(self.subsets, self.comp)], dtype=float)
        self.match = []
        from vlib import signed_matchings
        self.match = signed_matchings(n)

    def density(self, t, lam):
        """t: (P, N) positive; lam: (V,) sources at the path vertices. Returns D(t) (P,)."""
        n, N, C = self.n, self.N, self.C
        P = t.shape[0]
        M = np.einsum("re,pe,se->prs", C, t, C)
        Minv = np.linalg.inv(M)
        U = np.linalg.det(M)
        Q = np.einsum("ae,pab,bf->pef", C, Minv, C)
        Lap = np.einsum("ev,pe,ew->pvw", self.B, 1.0 / t, self.B)
        phi = np.linalg.solve(Lap, np.broadcast_to(lam, (P, len(lam)))[..., None])[..., 0]
        y = np.einsum("ev,pv->pe", self.B, phi) / t
        Psi = C[None, :, :] * y[:, None, :]                 # (P, n, N)
        m = 2 * n
        Qm = Q[:, :m, :m]
        Cm = C[:, :m]

        def Aform(r, s):
            return -0.5 * Qm * (np.outer(Cm[r], Cm[s]) - np.outer(Cm[s], Cm[r]))[None]
        A = {(r, s): Aform(r, s) for r in range(n) for s in range(n) if r < s}
        # Pf components on n-subsets
        pf = np.zeros((P, len(self.subsets)))
        for k, S in enumerate(self.subsets):
            acc = np.zeros(P)
            # (wedge of n/2 two-forms)_S = sum over ways to split S into ordered pairs assigned to the factors
            for mt, sg in self.match:
                pairs = [(i, j) for (i, j) in mt]
                # coefficient of dt_S in wedge_k A_{pairs[k]}: sum over set partitions of S into pairs,
                # assigned in order to the factors, with sign of the arrangement
                acc += sg * self._wedge_coeff(A, pairs, S)
            pf[:, k] = acc
        dets = np.linalg.det(np.stack([Psi[:, :, list(Sc)] for Sc in self.comp], axis=1))   # (P, nsub)
        top = (pf * dets * self.ssign[None, :]).sum(axis=1)
        e0 = N - 1
        D = top / t[:, e0] * U ** (-1.5)
        return D

    def _wedge_coeff(self, A, pairs, S):
        """coefficient of dt_S (S sorted) in A_{pairs[0]} ^ ... ^ A_{pairs[-1]} (2-forms)."""
        S = list(S)
        if not pairs:
            return 1.0
        out = 0.0
        first = S[0]
        # the factor that takes S[0]: any factor; to keep it simple, enumerate perfect matchings of S
        # assigned to the factors in all orders
        tot = 0.0
        for mt in _matchings_of(S):
            for perm in itertools.permutations(range(len(pairs))):
                seq = []
                prod = 1.0
                for fac, pr in zip(perm, mt):
                    pass
                # assign matching pair mt[q] to factor perm[q]
                order = [None] * len(pairs)
                for q, f in enumerate(perm):
                    order[f] = mt[q]
                flat = [x for pr in order for x in pr]
                sg = perm_sign([S.index(x) for x in flat])
                val = sg
                for f, (i, j) in enumerate(order):
                    val = val * A[pairs[f]][:, i, j]
                tot = tot + val
        return tot


def _matchings_of(seq):
    if not seq:
        yield []
        return
    a = seq[0]
    for k in range(1, len(seq)):
        rest = seq[1:k] + seq[k + 1:]
        for m in _matchings_of(rest):
            yield [(a, seq[k])] + m


def selftest(ff, rng):
    n = ff.n
    for p in range(2):
        t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(ff.N)]
        lam = [Fr(rng.randint(-5, 5)) for _ in range(n + 1)]
        D, U = brute_density(ff.G, [list(map(int, r)) for r in ff.C], t, [lam], 1)
        exact = float(D) / math.sqrt(float(U))
        fl = ff.density(np.array([[float(x) for x in t]]), np.array([float(x) for x in lam]))[0]
        print("   selftest n=%d: float density %.15e vs exact %.15e (rel %.1e)" % (n, fl, exact, abs(fl - exact) / abs(exact)))


RATIO_CUT = 32.0     # points with max t / min t > e^32 are not evaluated (counted in lattice.dropped)


def lattice(ff, lam, h, X, alpha, chunk=200000, ywt=None):
    """trapezoid in x, s = alpha sinh(x); returns K (and optionally int D y_e for all e)."""
    N = ff.N
    xs = np.arange(-X, X + 1e-12, h)
    s = alpha * np.sinh(xs)
    w = alpha * np.cosh(xs) * h
    grids = list(itertools.product(range(len(xs)), repeat=N - 1))
    idx = np.array(grids, dtype=np.int32)
    tot = 0.0
    toty = np.zeros(N)
    lattice.dropped = 0
    for c0 in range(0, len(idx), chunk):
        ii = idx[c0:c0 + chunk]
        sv = s[ii]
        tt = np.concatenate([np.exp(sv), np.ones((len(ii), 1))], axis=1)
        keep = np.log(tt.max(axis=1) / tt.min(axis=1)) < RATIO_CUT
        tt = tt[keep]
        ii = ii[keep]
        lattice.dropped += int((~keep).sum())
        wt = np.prod(w[ii], axis=1) * np.prod(tt, axis=1)
        D = ff.density(tt, lam)
        val = D * wt
        tot += np.nansum(val)
        if ywt is not None:
            Lap = np.einsum("ev,pe,ew->pvw", ff.B, 1.0 / tt, ff.B)
            phi = np.linalg.solve(Lap, np.broadcast_to(lam, (len(ii), len(lam)))[..., None])[..., 0]
            y = np.einsum("ev,pv->pe", ff.B, phi) / tt
            toty += np.nansum(val[:, None] * y, axis=0)
    return tot, toty


def chart_gl(ff, lam, m, chunk=40000):
    n = ff.n
    xg, wg = np.polynomial.legendre.leggauss(m)
    th = (xg + 1) * np.pi / 4           # (0, pi/2)
    wth = wg * np.pi / 4
    nd = 2 * n
    tot = 0.0
    idx = np.array(list(itertools.product(range(m), repeat=nd)), dtype=np.int32)
    for c0 in range(0, len(idx), chunk):
        ii = idx[c0:c0 + chunk]
        TH = th[ii]
        W = np.prod(wth[ii], axis=1)
        S2 = np.sin(TH) ** 2
        jac = np.prod(np.sin(2 * TH), axis=1)        # d(sin^2)/dtheta
        # variables: u_1, v_1, ..., u_{n-1}, v_{n-1}, u_n, y  with v_n = (1 - u_n) y
        u = {k: S2[:, 2 * k - 2] for k in range(1, n + 1)}
        v = {k: S2[:, 2 * k - 1] for k in range(1, n)}
        v[n] = (1 - u[n]) * S2[:, 2 * n - 1]
        jac = jac * (1 - u[n])
        t = np.stack(t_of_r(n, u, v), axis=1)
        D = ff.density(t, lam)
        absJ = np.ones(len(ii))
        for k in range(1, n):
            absJ *= u[k + 1] ** (2 * k) / (v[k] ** 2 * (1 - v[k]) ** (2 * k + 1))
        tot += np.sum(D * absJ * jac * W)
    return tot


if __name__ == "__main__":
    rng = random.Random(3)
    what = sys.argv[1]
    if what == "selftest":
        for n in (2, 4):
            selftest(FloatFan(n), rng)
    elif what == "lattice2":
        ff = FloatFan(2)
        for lam, exact in [([1.0, 0.0, 1.0], math.pi / 3), ([3.0, -2.0, 5.0], 5 * math.pi)]:
            for (h, X, alpha) in [tuple(float(a) for a in arg.split(",")) for arg in sys.argv[2:]]:
                t0 = time.time()
                K, _ = lattice(ff, np.array(lam), h, X, alpha)
                print("lattice n=2 lam=%s h=%.2f X=%.2f alpha=%.1f (%d^4 nodes, %d with max/min t > e^%g skipped): K = %.12f, exact %.12f, rel.diff %.2e (%.0f s)"
                      % (lam, h, X, alpha, len(np.arange(-X, X + 1e-12, h)), lattice.dropped, RATIO_CUT, K, exact, abs(K - exact) / exact, time.time() - t0))
    elif what == "lattice2y":
        ff = FloatFan(2)
        lam = np.array([1.0, 0.0, 1.0])
        h, X, alpha = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        t0 = time.time()
        K, Ky = lattice(ff, lam, h, X, alpha, ywt=True)
        # claim b1 n=2 at l1 = l3 = 1: order-1 coefficient of z_e is -int D y_e:
        # -(pi/30) [ (4-1) z_a1 + (1-4) z_a2 + (6+1) z_b1 + (3+3) z_b2 + (1+6) z_b3 ]
        cl = {0: 3, 2: -3, 1: 7, 3: 6, 4: 7}      # edge order a1, b1, a2, b2, b3
        print("lattice n=2 lam=(1,0,1) h=%.2f X=%.1f alpha=%.1f: K = %.12f (exact pi/3 = %.12f)  (%.0f s)" % (h, X, alpha, K, math.pi / 3, time.time() - t0))
        names = ["z_a1", "z_b1", "z_a2", "z_b2", "z_b3"]
        for e in range(5):
            ex = -math.pi / 30 * cl[e]
            print("   coefficient of %s: -int D y_e = %.10f, claim -(pi/30)*%d = %.10f, rel.diff %.2e" % (names[e], -Ky[e], cl[e], ex, abs(-Ky[e] - ex) / abs(ex)))
    elif what == "chart4":
        ff = FloatFan(4)
        for lam, exact in [([1.0, 0.0, 0.0, 0.0, 1.0], math.pi ** 2 / 180), ([2.0, -1.0, 3.0, 1.0, -2.0], -2.4 * math.pi ** 2)]:
            for m in [int(x) for x in sys.argv[2:]]:
                t0 = time.time()
                K = chart_gl(ff, np.array(lam), m)
                print("chart GL n=4 lam=%s m=%d (%d nodes): K = %.12f, exact %.12f, rel.diff %.2e (%.0f s)"
                      % (lam, m, m ** 8, K, exact, abs(K - exact) / abs(exact), time.time() - t0))
    elif what == "chart2":
        ff = FloatFan(2)
        for lam, exact in [([1.0, 0.0, 1.0], math.pi / 3), ([3.0, -2.0, 5.0], 5 * math.pi)]:
            for m in [int(x) for x in sys.argv[2:]]:
                K = chart_gl(ff, np.array(lam), m)
                print("chart GL n=2 lam=%s m=%d: K = %.13f, exact %.13f, rel.diff %.2e" % (lam, m, K, exact, abs(K - exact) / abs(exact)))
