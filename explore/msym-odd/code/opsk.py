"""Q_0 and the two-loop theta-graph operators on R x C^H with k odd coordinates, generic colour
tensors (arm convention).  Used for 7d (H = 3, k = 1), 3d (H = 1, k = 3), and for 5d (H = 2, k = 2)
as a cross-check of ops5.py.

Qhat_0 C^A = (1/2) f_{ABC} C^B C^C on the components:  Qhat_0 c_S = (-1)^{|S|} [theta_S-coefficient].

Theta graph with O = N (root), arms P = A, B, C from O to N' with interior vertices v^P_1..v^P_{m_P}
(m_P = n_P - 1); the contracted letter j sits on the spoke of arm j (z_j on the first edge of the
path, oriented from O).  Arm tensor:  arm_P[R] = f_{A_P D_1 R_1} C^{D_1}(v_1) f_{R_1 D_2 R_2} C^{D_2}(v_2) ...
with R = R_{m_P} (R = A_P if m_P = 0);
T^arm_Gamma = f_{R_A R_B R_C} arm_A[R_A] arm_B[R_B] arm_C[R_C],  outputs ordered v^A_1.., v^B_1.., v^C_1...
Qhat_Gamma = int d^k theta int prod_j dz_j [d_z^k K_Gamma / pi] T^arm_Gamma prod_j delta/delta C^{A_j}(x + z_j, theta).
"""
import os
from fractions import Fraction as Fr

import sympy as sp

MUTATE = os.environ.get('MUTATE', '')       # control runs only: 'scale' doubles the z-derivative weights

import htk_path  # noqa: F401
import htk
from salgk import SAlg
from theta_kernels import ThetaKernel


class TheoryK:
    def __init__(self, alg, n, f):
        self.S = alg
        self.n, self.f = n, f
        self.fby = {}
        for (A, B, C), v in f.items():
            self.fby.setdefault(A, []).append((B, C, v))
        self._q0 = {}

    def q0_letter(self, l):
        if l in self._q0:
            return self._q0[l]
        S = self.S
        t, A, j = l[0], l[1], l[2:]
        if any(j):
            a = next(i for i in range(3) if j[i])
            jj = list(j); jj[a] -= 1
            r = S.dpar(self.q0_letter((t, A) + tuple(jj)), a)
        else:
            R = {}
            for (B, C, v) in self.fby.get(A, []):
                R = S.add(R, S.scale(S.mul(S.csuper(B), S.csuper(C)), v / 2))
            mask = t - 1
            r = S.scale(S.theta_coeff(R, mask), (-1) ** bin(mask).count('1'))
        self._q0[l] = r
        return r

    def Q0(self, p):
        S = self.S
        r = {}
        for l in S.letters_of(p):
            d = S.lderiv(p, l)
            if d:
                r = S.add(r, S.mul(self.q0_letter(l), d))
        return r


class ThetaOp:
    """Kernel data of one theta graph (O = root = N)."""

    def __init__(self, G, paths, H, name):
        self.TK = ThetaKernel(G, paths, H, name)
        self.TK.check_pointwise(npts=3, seed=23)
        TK = self.TK
        self.H = H
        self.outv = [v for P in range(3) for v in TK.interior[P]]
        self.arm_sizes = [len(TK.interior[P]) for P in range(3)]
        zs = [[sp.Symbol('z%d_%d' % (j, i)) for i in range(H)] for j in range(3)]
        zedge = {}
        for P in range(3):
            e = TK.paths[P][0]
            sg = TK.signs[P][0]          # +1: spoke oriented from O
            zedge[e] = [sg * zs[P][i] for i in range(H)]
        X = TK.displacement(zedge)
        self.cz = []
        for v in self.outv:
            row = []
            for j in range(3):
                cs = [sp.expand(sp.diff(X[v][i], zs[j][i])) for i in range(H)]
                assert all(sp.expand(c - cs[0]) == 0 for c in cs)
                row.append(cs[0])
            self.cz.append(row)
        self.lam = [TK.lam[v] for v in self.outv]
        zero = {TK.lam[TK.Nprime][i]: 0 for i in range(H)} if TK.Nprime < G.nv - 1 else {}
        self.D = sp.expand(TK.Dtheta().subs(zero))
        self.s = TK.s
        self._W = {}

    def weights(self, ks):
        if ks in self._W:
            return self._W[ks]
        H, nout = self.H, len(self.outv)
        expr = self.s * self.D / 2
        for j in range(3):
            for a in range(3):
                if a >= H:
                    assert ks[j][a] == 0
                    continue
                if ks[j][a]:
                    expr *= sum(self.lam[v][a] * self.cz[v][j] for v in range(nout)) ** ks[j][a]
        res = self.TK.integrate_poly(expr)
        out = {}
        if res:
            assert set(res) == {1}, res
            allL = [self.lam[v][a] for v in range(nout) for a in range(H)]
            P = sp.Poly(sp.expand(res[1]), *allL)
            for mon, c in P.terms():
                m = tuple(tuple(mon[H * v + a] for a in range(H)) + (0,) * (3 - H) for v in range(nout))
                c = sp.Rational(c)
                out[m] = Fr(int(c.p), int(c.q))
        if MUTATE == 'scale' and any(any(k_) for k_ in ks):
            out = {m: 2 * c for m, c in out.items()}
        self._W[ks] = out
        return out


class TwoLoopK:
    def __init__(self, theory, ops):
        self.th = theory
        self.S = theory.S
        self.fby = theory.fby
        self.ops = ops
        self._T = {}
        self._cs = {}

    def cs(self, A, m):
        key = (A, m)
        if key not in self._cs:
            self._cs[key] = self.S.csuper(A, m)
        return self._cs[key]

    def arm(self, A0, ms):
        S = self.S
        vec = {A0: S.const(1)}
        for m in ms:
            out = {}
            for R, p in vec.items():
                for (D, Rn, v) in self.fby.get(R, []):
                    out[Rn] = S.add(out.get(Rn, {}), S.scale(S.mul(p, self.cs(D, m)), v))
            vec = out
        return vec

    def T(self, kind, As, m):
        key = (kind, As, m)
        if key in self._T:
            return self._T[key]
        S = self.S
        sizes = self.ops[kind].arm_sizes
        ms, k0 = [], 0
        for P in range(3):
            ms.append(m[k0:k0 + sizes[P]]); k0 += sizes[P]
        arms = [self.arm(As[P], ms[P]) for P in range(3)]
        res = {}
        for RA, pa in arms[0].items():
            for (RB, RC, v) in self.fby.get(RA, []):
                if RB in arms[1] and RC in arms[2]:
                    res = S.add(res, S.scale(S.mul(S.mul(pa, arms[1][RB]), arms[2][RC]), v))
        self._T[key] = res
        return res

    def G(self, kind, As, ks):
        S = self.S
        r = {}
        for m, w in self.ops[kind].weights(ks).items():
            t = self.T(kind, As, m)
            if t:
                r = S.add(r, S.scale(t, w))
        return r

    def Q2(self, kind, p):
        S = self.S
        r = {}
        ls = [l for l in S.letters_of(p) if l[0] != 0]
        for l3 in ls:
            d3 = S.sderiv(p, l3)
            if not d3:
                continue
            for l2 in ls:
                d23 = S.sderiv(d3, l2)
                if not d23:
                    continue
                for l1 in ls:
                    d123 = S.sderiv(d23, l1)
                    if not d123:
                        continue
                    G = self.G(kind, (l1[1], l2[1], l3[1]), (l1[2:], l2[2:], l3[2:]))
                    if G:
                        r = S.add(r, S.theta_coeff(S.mul(G, d123), S.full))
        return r
