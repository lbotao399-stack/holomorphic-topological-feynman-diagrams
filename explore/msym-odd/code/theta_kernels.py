"""Two-loop theta kernels at T = 1 (S3 Theorem thm-theta, eq. R-Ktheta, adapted to H holomorphic
directions and to graphs whose root O is the trivalent vertex N).

For a theta graph with trivalent vertices N = root (O) and N', paths P = A, B, C from N to N' with
n_P edges and interior vertices P_1..P_{n_P - 1}:

  R_P = sum_j t^P_j,  x^P_j = t^P_j / R_P,  w_P = (1/R_P) / sum_P' (1/R_P'),
  y^P_j = l_P + Q^P_j,  Q^P_j = sum_{k<j} lam^P_k,  mhat_P = sum_j x^P_j Q^P_j,
  l_P = w_P (lam_N + sum_P' mhat_P') - mhat_P       (S3 eq-series, eq-parallel),
  y_e = +- y^P_j (edge oriented along / against the path),

  K_Gamma(lam; z) = s_Gamma (D_Gamma / 2) int_{Delta^2} d^2w prod_P w_P^{H + 1/2 - n_P}
                    prod_P int_{Delta^{n_P - 1}} dx^P  exp(sum_v lam_v . X_v(x, w; z)),
  D_Gamma = det [[Lam_A, Lam_B, 0], [0, Lam_B, Lam_C]],  Lam_P = (lam^P_1 .. lam^P_{n_P-1}) (H x (n_P-1)),
  -sum_e y_e . z_e = sum_v lam_v . X_v.

s_Gamma = +-1 is determined by comparison with the brute-force density of htk.brute (Kirchhoff
solution at rational Schwinger times, no reduction):  g(xi) = D(t(xi)) |det[t ; dt/dxi]|.
"""
import random
import sys
from fractions import Fraction as Fr
from functools import lru_cache

import sympy as sp

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import brute
from htk.simplex import dirichlet, PiNumber


def to_fr(x):
    x = sp.Rational(x)
    return Fr(int(x.p), int(x.q))


class ThetaKernel:
    def __init__(self, G, paths, H, name=""):
        """G: htk Graph with root = N; paths: three lists of edge ids from N to N' (A, B, C)."""
        self.G, self.H, self.name = G, H, name
        self.N = G.root
        self.paths, self.signs, self.interior = [], [], []
        for pe in paths:
            v, sg, ints = self.N, [], []
            for pos, e in enumerate(pe):
                a, b = G.edges[e]
                if a == v:
                    sg.append(1); v = b
                elif b == v:
                    sg.append(-1); v = a
                else:
                    raise ValueError("path broken at edge %d" % e)
                if pos < len(pe) - 1:
                    ints.append(v)
            self.paths.append(list(pe)); self.signs.append(sg); self.interior.append(ints)
        self.Nprime = v
        self.nP = [len(p) for p in self.paths]
        assert sorted(e for p in self.paths for e in p) == list(range(G.ne))
        self.alpha = [Fr(2 * H + 1, 2) - n for n in self.nP]
        # symbols
        self.w = sp.symbols('wA wB wC', positive=True)
        self.x = [sp.symbols('x%s1:%d' % ("ABC"[i], n + 1), positive=True) for i, n in enumerate(self.nP)]
        g = G.nv - 1
        self.lam = [[sp.Symbol('l%d_%d' % (v, i)) for i in range(H)] for v in range(g)]
        self.s = None

    # -------------------------------------------------------------- currents
    def currents(self, lam=None, w=None, x=None):
        """y[e][i] as sympy expressions (lam: list over non-root vertices of H-vectors)."""
        H = self.H
        lam = self.lam if lam is None else lam
        w = self.w if w is None else w
        x = self.x if x is None else x
        g = self.G.nv - 1
        lamN = [-sum(lam[v][i] for v in range(g)) for i in range(H)]
        Q, mhat = [], []
        for P in range(3):
            Qs, acc = [[0] * H], [0] * H
            for u in self.interior[P]:
                acc = [acc[i] + lam[u][i] for i in range(H)]
                Qs.append(list(acc))
            Q.append(Qs)
            mhat.append([sum(x[P][j] * Qs[j][i] for j in range(self.nP[P])) for i in range(H)])
        tot = [lamN[i] + sum(mhat[P][i] for P in range(3)) for i in range(H)]
        y = [None] * self.G.ne
        for P in range(3):
            l = [w[P] * tot[i] - mhat[P][i] for i in range(H)]
            for j, e in enumerate(self.paths[P]):
                yy = [l[i] + Q[P][j][i] for i in range(H)]
                y[e] = yy if self.signs[P][j] > 0 else [-c for c in yy]
        return y

    def Dtheta(self, lam=None):
        lam = self.lam if lam is None else lam
        H = self.H
        A, B, C = self.interior
        rows = []
        for i in range(H):
            rows.append([lam[v][i] for v in A] + [lam[v][i] for v in B] + [0] * len(C))
        for i in range(H):
            rows.append([0] * len(A) + [lam[v][i] for v in B] + [lam[v][i] for v in C])
        return sp.Matrix(rows).det() if rows else sp.Integer(1)

    # -------------------------------------------------------------- map to Schwinger times
    def coords(self):
        """independent coordinates: wB, wC, then x^P_2..x^P_{n_P} for P = A, B, C."""
        out = [self.w[1], self.w[2]]
        for P in range(3):
            out += list(self.x[P][1:])
        return out

    def t_of(self):
        """t_e as sympy expressions in the independent coordinates (gauge sum_P w_P = 1, R_P = 1/w_P)."""
        wB, wC = self.w[1], self.w[2]
        wv = [1 - wB - wC, wB, wC]
        t = [None] * self.G.ne
        sub = {self.w[0]: wv[0]}
        for P in range(3):
            xs = list(self.x[P])
            x1 = 1 - sum(xs[1:])
            sub[xs[0]] = x1
            for j, e in enumerate(self.paths[P]):
                xj = x1 if j == 0 else xs[j]
                t[e] = xj / wv[P]
        return t, sub

    def random_point(self, rng):
        def simplex(k):
            ns = [rng.randint(1, 30) for _ in range(k)]
            return [Fr(a, sum(ns)) for a in ns]
        wv = simplex(3)
        xv = [simplex(n) for n in self.nP]
        return wv, xv

    # -------------------------------------------------------------- pointwise check against brute
    def check_pointwise(self, npts=3, seed=0, lam_num=None):
        rng = random.Random(seed)
        t_expr, _sub = self.t_of()
        cs = self.coords()
        Jm = sp.Matrix([t_expr] + [[sp.diff(te, c) for te in t_expr] for c in cs])
        signs = set()
        g = self.G.nv - 1
        for k in range(npts):
            lamn = lam_num or [[Fr(rng.randint(-6, 6)) for _ in range(self.H)] for _ in range(g)]
            wv, xv = self.random_point(rng)
            vals = {self.w[1]: wv[1], self.w[2]: wv[2]}
            for P in range(3):
                for j in range(1, self.nP[P]):
                    vals[self.x[P][j]] = xv[P][j]
            spv = {kk: sp.Rational(v.numerator, v.denominator) for kk, v in vals.items()}
            t = [to_fr(te.subs(spv)) for te in t_expr]
            J = to_fr(Jm.subs(spv).det())
            Dt = brute.density(self.G, t, lamn, T=1)
            lhs = Dt.scaled(abs(J))
            # claimed density s (D/2) prod w^alpha: coeff prod w^{H - n_P + 1}, radicand prod w
            Dn = to_fr(self.Dtheta([[sp.Rational(c.numerator, c.denominator) for c in v] for v in lamn]))
            coeff = Dn / 2
            rad = Fr(1)
            for P in range(3):
                coeff *= wv[P] ** (self.H - self.nP[P] + 1)
                rad *= wv[P]
            rhs = brute.Density(coeff, rad)
            if Dn == 0:
                assert lhs.coeff == 0
                continue
            if lhs.square() != rhs.square():
                raise AssertionError("density magnitude mismatch %s: %r vs %r" % (self.name, lhs, rhs))
            signs.add(lhs.sign() * rhs.sign())
            # currents
            lamsp = [[sp.Rational(c.numerator, c.denominator) for c in v] for v in lamn]
            wsp = [sp.Rational(c.numerator, c.denominator) for c in wv]
            xsp = [[sp.Rational(c.numerator, c.denominator) for c in xx] for xx in xv]
            yc = [[to_fr(c) for c in ye] for ye in self.currents(lamsp, wsp, xsp)]
            if yc != brute.currents(self.G, t, lamn):
                raise AssertionError("currents mismatch %s" % self.name)
        if len(signs) != 1:
            raise AssertionError("sign not constant: %s" % signs)
        self.s = signs.pop()
        return self.s

    # -------------------------------------------------------------- integration
    def displacement(self, zedge):
        """X_v (list over non-root vertices of H-vectors of sympy expressions) with
        -sum_e y_e . z_e = sum_v lam_v . X_v for z_e = zedge[e] (H-vectors of sympy expressions)."""
        y = self.currents()
        E = -sum(sum(y[e][i] * zedge[e][i] for i in range(self.H)) for e in zedge)
        E = sp.expand(E)
        g = self.G.nv - 1
        return [[sp.expand(sp.diff(E, self.lam[v][i])) for i in range(self.H)] for v in range(g)]

    def integrate_poly(self, poly_expr):
        """int prod w^alpha prod dx  poly_expr  (poly in w, x with coefficients in other symbols).
        Returns {pi power: sympy coefficient}."""
        gens = list(self.w) + [xx for P in range(3) for xx in self.x[P]]
        Pl = sp.Poly(sp.expand(poly_expr), *gens)
        acc = {}
        for mon, c in Pl.terms():
            q, kpi = Fr(1), 0
            dq, dk = dirichlet(tuple(self.alpha[P] + mon[P] for P in range(3)))
            q *= dq; kpi += dk
            off = 3
            for P in range(3):
                n = self.nP[P]
                if n >= 2:
                    dq, dk = dirichlet(tuple(Fr(mon[off + j]) for j in range(n)))
                    q *= dq; kpi += dk
                off += n
            acc[kpi] = acc.get(kpi, 0) + sp.Rational(q.numerator, q.denominator) * c
        return {k: sp.expand(v) for k, v in acc.items() if sp.expand(v) != 0}

    def K0(self):
        """K(lam; 0) symbolic: {pi power: coefficient}."""
        assert self.s is not None
        res = self.integrate_poly(sp.Integer(1))
        D = self.Dtheta()
        return {k: sp.expand(self.s * D / 2 * v) for k, v in res.items()}
