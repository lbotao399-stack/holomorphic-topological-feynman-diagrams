"""Two-loop operators hat Q_A (double box) and hat Q_B (spider) of holomorphic CS on C^{3|2}.

hat Q_Gamma = int d^2theta int prod_j d^6z_j int dmu_Gamma  (colour tensor)  J_Gamma(d)
              [four output superfields at X_v(z)]  delta_{A1}(z_1) delta_{A2}(z_2) delta_{A3}(z_3)

(A) tree = path 1-2-3-4-5, spokes at 1,3,5 (letters 1,2,3); outputs at 1,2,4,5
    colour  f_{A1 D1 R1} f_{R1 D2 R2} f_{R2 A2 R3} f_{R3 D4 R4} f_{R4 D5 A3}
    J_A = rho^2 [q2 <124> + q3 <125>] [(1-s1) <145> + (1-s2) <245>]
(B) tree = centre q with arms u | w | a-b (a = v2 next to q, b = v1), spokes at u,w,b (letters 1,2,3)
    colour  f_{A1 Du Ru} f_{A2 Dw Rw} f_{Ru Rw Rq} f_{Rq Da R'} f_{R' Db A3}
    J_B = rho (1-rho) [q2 <uwa> + q3 <uwb>] [(1-ru) <uab> + (1-rw) <wab>]
<ijk> = eps^{abc} d_a^{(i)} d_b^{(j)} d_c^{(k)}.
"""
from fractions import Fraction as Fr
import sympy as sp

from salg import add, csuper, letters_of, mul, scale, sderiv, theta_coeff

rho, s1, s2, ru, rw = sp.symbols('rho s1 s2 ru rw')
q1, q2, q3 = sp.symbols('q1 q2 q3')
LAM = [[sp.Symbol(f'l{v}{a}') for a in range(3)] for v in range(4)]     # four outputs
ZW = lambda z: (1 - rho) * z[0] + rho * z[1] - z[2]


def br(i, j, k):
    return sp.Matrix([LAM[i], LAM[j], LAM[k]]).det()


# output order (A): (1,2,4,5) -> 0,1,2,3 ; (B): (u,w,a,b) -> 0,1,2,3
JA = rho ** 2 * (q2 * br(0, 1, 2) + q3 * br(0, 1, 3)) * ((1 - s1) * br(0, 2, 3) + (1 - s2) * br(1, 2, 3))
JB = rho * (1 - rho) * (q2 * br(0, 1, 2) + q3 * br(0, 1, 3)) * ((1 - ru) * br(0, 2, 3) + (1 - rw) * br(1, 2, 3))


def XA(z):
    W = ZW(z)
    return [q1 * (1 - s1) * W + s1 * z[0] + (1 - s1) * z[2],
            q1 * (1 - s2) * W + s2 * z[0] + (1 - s2) * z[2],
            q2 * W + z[2], q3 * W + z[2]]


def XB(z):
    W = ZW(z)
    return [q1 * (1 - ru) * W + ru * z[0] + (1 - ru) * z[2],
            q1 * (1 - rw) * W + rw * z[1] + (1 - rw) * z[2],
            q2 * W + z[2], q3 * W + z[2]]


def integrate_measure(expr, kind):
    """Exact integral of a polynomial over the ratio domain of (A) or (B)."""
    expr = sp.expand(expr)
    if kind == 'A':
        vs = (rho, s1, s2, q1, q2, q3)
    else:
        vs = (rho, ru, rw, q1, q2, q3)
    P = sp.Poly(expr, *vs)
    tot = sp.Integer(0)
    for mon, c in P.terms():
        r, a, b, e1, e2, e3 = mon
        w = sp.Rational(1, r + 1)
        if kind == 'A':
            w *= sp.Rational(1, (b + 1) * (a + b + 2))                 # 1 > s1 > s2 > 0
        else:
            w *= sp.Rational(1, (a + 1) * (b + 1))                     # ru, rw in (0,1)
        w *= sp.Rational(1, (e3 + 1) * (e2 + e3 + 2) * (e1 + e2 + e3 + 3))   # 1 > q1 > q2 > q3 > 0
        tot += c * w
    return sp.expand(tot)


_W = {}


def weights(kind, ks):
    """{(m_0,..,m_3): Fraction}: d_z^{k} at z = 0 of  int dmu J(lam) exp(sum_v lam_v . X_v(z))."""
    key = (kind, ks)
    if key in _W:
        return _W[key]
    J = JA if kind == 'A' else JB
    Xf = XA if kind == 'A' else XB
    expr = J
    for a in range(3):
        for j in range(3):
            if ks[j][a]:
                z = [sp.Integer(0)] * 3
                z[j] = sp.Integer(1)
                X = Xf(z)                                  # coefficient of z_j in X_v
                expr *= sum(LAM[v][a] * X[v] for v in range(4)) ** ks[j][a]
    allL = [LAM[v][a] for v in range(4) for a in range(3)]
    P = sp.Poly(sp.expand(expr), *allL)
    out = {}
    for mon, coeff in P.terms():
        val = integrate_measure(coeff, kind)
        if val != 0:
            val = sp.Rational(val)
            m = tuple(tuple(mon[3 * v + a] for a in range(3)) for v in range(4))
            out[m] = Fr(int(val.p), int(val.q))
    _W[key] = out
    return out


class TwoLoop:
    def __init__(self, theory):
        self.th = theory
        self.fby = theory.fby
        self._T = {}
        self._cs = {}

    def cs(self, A, m):
        key = (A, m)
        if key not in self._cs:
            self._cs[key] = csuper(A, m)
        return self._cs[key]

    def step(self, vec, m):
        """vec[R] -> out[R'] = sum_{R,D} vec[R] f_{R D R'} C^D_{,m}"""
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
        for m, w in weights(kind, ks).items():
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
