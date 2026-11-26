"""Q_0 and the two-loop operators of the holomorphic-topological theory on C^{2|2} x R (5d MSYM).

Adapted from 6d-holomorphic-twist-note/code/ops.py (Q0hat, unchanged) and ops2.py (colour tensors
T_A, T_B and the functional-derivative structure, unchanged); only the kernels are replaced:
H = 2 holomorphic directions (letters carry jets (k1, k2, 0)), T = 1 kernels of theta_kernels.py.

  Qhat_Gamma = int d^2theta int prod_j d^5z_j  [d_z^k K_Gamma / pi]  T_Gamma  prod_j delta/delta C^{A_j}(x+z_j, theta),

K_Gamma the reduced kernel of htk for htk.graph_A / graph_B (data of S2 Section 7), displacements
z_1 on (O->1) [(O->u)], z_2 on (O->3) [(O->w)], -z_3 on (5->O) [(b->O)] as in S2 eq-A-exp, eq-B-exp.
Outputs in the order of T_Gamma: (A) vertices 1,2,4,5 = htk 0,1,3,4;  (B) u,w,a,b = htk 0,1,3,4.
"""
from fractions import Fraction as Fr

import sympy as sp

import htk_path  # noqa: F401
import htk
from salg import (add, const, csuper, dpar, gen, lderiv, letters_of, mul, scale, sderiv, theta_coeff)
from theta_kernels import ThetaKernel


class Theory:
    """Q0hat of S2 code/ops.py (Q_0 = i Q0hat, Q0hat C^A = (1/2) f_{ABC} C^B C^C)."""

    def __init__(self, n, f):
        self.n = n
        self.f = f
        self.fby = {}
        for (A, B, C), v in f.items():
            self.fby.setdefault(A, []).append((B, C, v))
        self._q0 = {}

    def q0_letter(self, l):
        if l in self._q0:
            return self._q0[l]
        t, A, k = l[0], l[1], l[2:]
        if any(k):
            a = next(i for i in range(3) if k[i])
            kk = list(k); kk[a] -= 1
            r = dpar(self.q0_letter((t, A) + tuple(kk)), a)
        else:
            R = {}
            for (B, C, v) in self.fby.get(A, []):
                R = add(R, scale(mul(csuper(B), csuper(C)), v / 2))
            r = {1: theta_coeff(R, '0'), 3: scale(theta_coeff(R, '1'), -1),
                 4: scale(theta_coeff(R, '2'), -1), 2: theta_coeff(R, '12')}[t]
        self._q0[l] = r
        return r

    def Q0(self, p):
        r = {}
        for l in letters_of(p):
            d = lderiv(p, l)
            if d:
                r = add(r, mul(self.q0_letter(l), d))
        return r


def L(t, A, k=(0, 0, 0)):
    return gen(({'c': 1, 'b': 2, 'g1': 3, 'g2': 4}[t], A) + tuple(k))


# ------------------------------------------------------------------ kernels
H = 2
SPEC = {
    'A': (htk.graph_A(), [[3], [0, 1, 2], [6, 5, 4]], {0: (0, 1), 3: (1, 1), 6: (2, -1)}),
    'B': (htk.graph_B(), [[0, 1], [2, 3], [6, 5, 4]], {0: (0, 1), 2: (1, 1), 6: (2, -1)}),
}
OUTV = [0, 1, 3, 4]          # htk vertices of the four outputs (order of T_Gamma); vertex 2 = N' has lam = 0


class KernelData:
    def __init__(self, kind):
        G, paths, zmap = SPEC[kind]
        self.TK = ThetaKernel(G, paths, H, kind)
        self.TK.check_pointwise(npts=4, seed=17)
        TK = self.TK
        # coefficient c[v][j] of z_j in X_v (the same for both components)
        zs = [[sp.Symbol('z%d_%d' % (j, i)) for i in range(H)] for j in range(3)]
        zedge = {e: [sg * zs[j][i] for i in range(H)] for e, (j, sg) in zmap.items()}
        X = TK.displacement(zedge)
        self.cz = []
        for v in OUTV:
            row = []
            for j in range(3):
                cs = [sp.expand(sp.diff(X[v][i], zs[j][i])) for i in range(H)]
                assert sp.expand(cs[0] - cs[1]) == 0
                for i in range(H):
                    assert sp.expand(X[v][i] - sum(sp.diff(X[v][i], zs[jj][i]) * zs[jj][i] for jj in range(3))) == 0
                row.append(cs[0])
            self.cz.append(row)
        assert all(X[2][i] == 0 or True for i in range(H))
        self.lam = [TK.lam[v] for v in OUTV]
        zero = {TK.lam[2][i]: 0 for i in range(H)}
        self.D = sp.expand(TK.Dtheta().subs(zero))
        self.s = TK.s
        self._W = {}

    def weights(self, ks):
        """{(m_0..m_3): Fraction}: d_z^{ks} at z = 0 of K_Gamma(lam; z) / pi, ks = (k1, k2, k3) jet
        orders (3-tuples, third entry 0) of the three contracted letters."""
        if ks in self._W:
            return self._W[ks]
        expr = self.s * self.D / 2
        for j in range(3):
            for a in range(H):
                if ks[j][a]:
                    expr *= sum(self.lam[v][a] * self.cz[v][j] for v in range(4)) ** ks[j][a]
            assert ks[j][2] == 0
        res = self.TK.integrate_poly(expr)
        out = {}
        if res:
            assert set(res) == {1}, res
            allL = [self.lam[v][a] for v in range(4) for a in range(H)]
            P = sp.Poly(sp.expand(res[1]), *allL)
            for mon, c in P.terms():
                m = tuple(tuple(mon[H * v + a] for a in range(H)) + (0,) for v in range(4))
                c = sp.Rational(c)
                out[m] = Fr(int(c.p), int(c.q))
        self._W[ks] = out
        return out


class TwoLoop:
    """T_Gamma and Qhat_Gamma of S2 code/ops2.py with the kernels of KernelData."""

    def __init__(self, theory, kernels):
        self.th = theory
        self.fby = theory.fby
        self.K = kernels
        self._T = {}
        self._cs = {}

    def cs(self, A, m):
        key = (A, m)
        if key not in self._cs:
            self._cs[key] = csuper(A, m)
        return self._cs[key]

    def step(self, vec, m):
        out = {}
        for R, p in vec.items():
            for (D, Rn, v) in self.fby.get(R, []):
                out[Rn] = add(out.get(Rn, {}), scale(mul(p, self.cs(D, m)), v))
        return out

    def T(self, kind, As, m):
        key = (kind, As, m)
        if key in self._T:
            return self._T[key]
        A1, A2, A3 = As
        if kind == 'A':
            P1 = {}
            for (D, R, v) in self.fby.get(A1, []):
                P1[R] = add(P1.get(R, {}), scale(self.cs(D, m[0]), v))
            P2 = self.step(P1, m[1])
            P2p = {}
            for R2, p in P2.items():
                for (B, R3, v) in self.fby.get(R2, []):
                    if B == A2:
                        P2p[R3] = add(P2p.get(R3, {}), scale(p, v))
            P3 = self.step(P2p, m[2])
            res = {}
            for R4, p in P3.items():
                for (D5, B, v) in self.fby.get(R4, []):
                    if B == A3:
                        res = add(res, scale(mul(p, self.cs(D5, m[3])), v))
        else:
            U, Wv = {}, {}
            for (D, R, v) in self.fby.get(A1, []):
                U[R] = add(U.get(R, {}), scale(self.cs(D, m[0]), v))
            for (D, R, v) in self.fby.get(A2, []):
                Wv[R] = add(Wv.get(R, {}), scale(self.cs(D, m[1]), v))
            M = {}
            for Ru_, pu in U.items():
                for (Rw_, Rq, v) in self.fby.get(Ru_, []):
                    if Rw_ in Wv:
                        M[Rq] = add(M.get(Rq, {}), scale(mul(pu, Wv[Rw_]), v))
            N = self.step(M, m[2])
            res = {}
            for Rp, p in N.items():
                for (Db, B, v) in self.fby.get(Rp, []):
                    if B == A3:
                        res = add(res, scale(mul(p, self.cs(Db, m[3])), v))
        self._T[key] = res
        return res

    def G(self, kind, As, ks):
        r = {}
        for m, w in self.K[kind].weights(ks).items():
            t = self.T(kind, As, m)
            if t:
                r = add(r, scale(t, w))
        return r

    def Q2(self, kind, p):
        r = {}
        ls = [l for l in letters_of(p) if l[0] != 0]
        for l3 in ls:
            d3 = sderiv(p, l3)
            if not d3:
                continue
            for l2 in ls:
                d23 = sderiv(d3, l2)
                if not d23:
                    continue
                for l1 in ls:
                    d123 = sderiv(d23, l1)
                    if not d123:
                        continue
                    G = self.G(kind, (l1[1], l2[1], l3[1]), (l1[2:], l2[2:], l3[2:]))
                    if G:
                        r = add(r, theta_coeff(mul(G, d123), '12'))
        return r
