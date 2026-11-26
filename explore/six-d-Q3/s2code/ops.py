"""Q_0 and Q_1 of holomorphic CS on C^{3|2} as operators on polynomials in the letters.

Normalisation used here (hatted operators):
  Q_0 = i * Q0hat,   Q0hat C^A = (1/2) f_{ABC} C^B C^C,
  Q_1 = -(hbar/(16 pi^3)) * Q1hat,
  Q1hat = int d^2theta int d^6z_1 d^6z_2 int_{1>s1>s2>s3>0} d^3s  eps^{abc} f_{A1 D1 R1} f_{R1 D2 R2} f_{R2 D3 A2}
          d_a C^{D1}(x_1) d_b C^{D2}(x_2) d_c C^{D3}(x_3)  delta/delta C^{A1}(z_1) delta/delta C^{A2}(z_2),
  x_j = s_j z_1 + (1 - s_j) z_2.
"""
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import product

import sympy as sp

from salg import (EPS, add, const, csuper, dpar, gen, lderiv, letters_of, mul, scale, sderiv, theta_coeff)


class Theory:
    def __init__(self, n, f):
        self.n = n
        self.f = f                      # dict (A,B,C) -> Fraction, totally antisymmetric
        self.fby = {}
        for (A, B, C), v in f.items():
            self.fby.setdefault(A, []).append((B, C, v))
        self._q0 = {}
        self._T1 = {}
        self._W1 = {}

    # ---------------- Q0hat ----------------
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
            # Q0hat(c + th1 g1 + th2 g2 + th1 th2 b) = Q c - th1 Q g1 - th2 Q g2 + th1 th2 Q b
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

    # ---------------- Q1hat ----------------
    def W1(self, k1, k2):
        """Weights {(m1,m2,m3): Fraction}: d_{z1}^{k1} d_{z2}^{k2} at z=0 distributes jets m_j onto the outputs."""
        key = (k1, k2)
        if key in self._W1:
            return self._W1[key]
        s = sp.symbols('s1:4')
        u = [[sp.Symbol(f'u{j}{a}') for a in range(3)] for j in range(3)]
        expr = sp.Integer(1)
        for a in range(3):
            expr *= sum(s[j] * u[j][a] for j in range(3)) ** k1[a]
            expr *= sum((1 - s[j]) * u[j][a] for j in range(3)) ** k2[a]
        allu = [u[j][a] for j in range(3) for a in range(3)]
        poly = sp.Poly(sp.expand(expr), *allu)
        out = {}
        for mon, coeff in poly.terms():
            cp = sp.Poly(coeff, *s)
            val = sp.Integer(0)
            for (a, b, c), cc in cp.terms():
                val += cc * sp.Rational(1, (c + 1) * (b + c + 2) * (a + b + c + 3))
            if val != 0:
                m = tuple(tuple(mon[3 * j + a] for a in range(3)) for j in range(3))
                out[m] = Fr(int(val.p), int(val.q))
        self._W1[key] = out
        return out

    def T1(self, A1, A2, m):
        """eps^{abc} f_{A1 D1 R1} f_{R1 D2 R2} f_{R2 D3 A2} C^{D1}_{,a+m1} C^{D2}_{,b+m2} C^{D3}_{,c+m3} (superfields)."""
        key = (A1, A2, m)
        if key in self._T1:
            return self._T1[key]
        r = {}
        for (D1, R1, v1) in self.fby.get(A1, []):
            for (D2, R2, v2) in self.fby.get(R1, []):
                for (D3, B, v3) in self.fby.get(R2, []):
                    if B != A2:
                        continue
                    for (a, b, c), e in EPS.items():
                        ja = tuple(m[0][i] + (i == a) for i in range(3))
                        jb = tuple(m[1][i] + (i == b) for i in range(3))
                        jc = tuple(m[2][i] + (i == c) for i in range(3))
                        term = mul(mul(csuper(D1, ja), csuper(D2, jb)), csuper(D3, jc))
                        r = add(r, scale(term, v1 * v2 * v3 * e))
        self._T1[key] = r
        return r

    def G1(self, A1, A2, k1, k2):
        r = {}
        for m, w in self.W1(k1, k2).items():
            r = add(r, scale(self.T1(A1, A2, m), w))
        return r

    def Q1(self, p):
        r = {}
        ls = [l for l in letters_of(p) if l[0] != 0]
        for l2 in ls:
            d2 = sderiv(p, l2)
            if not d2:
                continue
            for l1 in ls:
                d12 = sderiv(d2, l1)
                if not d12:
                    continue
                G = self.G1(l1[1], l2[1], l1[2:], l2[2:])
                if G:
                    r = add(r, theta_coeff(mul(G, d12), '12'))
        return r


def L(t, A, k=(0, 0, 0)):
    return gen(({'c': 1, 'b': 2, 'g1': 3, 'g2': 4}[t], A) + tuple(k))
