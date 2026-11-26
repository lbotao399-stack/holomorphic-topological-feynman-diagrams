"""Exact layer for the K4 subdivisions on C^3 (H=3, T=0).

Coordinates xi = (sigma_{E,k} (k=1..a_E-1 for the long chains E), gA, gB, rA, rB, rC):
  chain E = path p=u_0, u_1, ..., u_a = q;  x_{E,j} = sigma_{E,j-1} - sigma_{E,j} (sigma_{E,0}=1, sigma_{E,a}=0),
  t_{E,j} = R_E x_{E,j};
  K4 conductances (Y-Delta chart at the K4 vertex Dc, triangle (A,B,C) = other cubic vertices, G = 1):
     edge Dc-v: g_v (gC = 1-gA-gB);  edge w-x opposite v: h_v = g_w g_x r_v/(1-r_v);  R_E = 1/conductance.
Reference density (brute force, no reduction):  g_ref(xi) = D(t(xi)) |det[t; dt/dxi]|,  D = htk.brute.density.
Sources: lambda = e1, e2, e3, -(e1+e2+e3) at the degree-two vertices d1<d2<d3<d4;  Delta = det(l_d1,l_d2,l_d3) = 1.
"""
import sys, os, random
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk.graph import Graph
from htk import brute
from htk.linalg import det as fdet
from graphs import K4SUB, skeleton_chains, deg

class K4Case:
    def __init__(self, which, Dc_index=0, perm=(0, 1, 2)):
        self.name, self.edges = K4SUB[which]
        self.nv = 8
        self.G = Graph(self.nv, self.edges)
        d = deg(self.edges, self.nv)
        self.d2 = [v for v in range(self.nv) if d[v] == 2]
        self.lam = {self.d2[0]: [1, 0, 0], self.d2[1]: [0, 1, 0], self.d2[2]: [0, 0, 1], self.d2[3]: [-1, -1, -1]}
        self.cubic, self.chains = skeleton_chains(self.edges, self.nv)
        self.Dc = self.cubic[Dc_index]
        o = [v for v in self.cubic if v != self.Dc]
        self.tri = [o[i] for i in perm]          # (A, B, C)
        # symbols
        self.gA, self.gB, self.rA, self.rB, self.rC = sp.symbols('gA gB rA rB rC')
        self.sig = {}
        for ci, (path, es) in enumerate(self.chains):
            for k in range(1, len(es)):
                self.sig[(ci, k)] = sp.Symbol('s%d_%d' % (ci, k))
        self.xi = [self.sig[k] for k in sorted(self.sig)] + [self.gA, self.gB, self.rA, self.rB, self.rC]
        self._tmap()

    def kedge(self, ci):
        path, es = self.chains[ci]
        return path[0], path[-1]

    def _tmap(self):
        A, B, C = self.tri
        g = {A: self.gA, B: self.gB, C: 1 - self.gA - self.gB}
        r = {A: self.rA, B: self.rB, C: self.rC}
        cond = {}
        for v in self.tri:
            cond[frozenset((self.Dc, v))] = g[v]
            w, x = [u for u in self.tri if u != v]
            cond[frozenset((w, x))] = g[w] * g[x] * r[v] / (1 - r[v])
        t = [None] * 10
        for ci, (path, es) in enumerate(self.chains):
            a = len(es)
            R = 1 / cond[frozenset((path[0], path[-1]))]
            s = [sp.Integer(1)] + [self.sig[(ci, k)] for k in range(1, a)] + [sp.Integer(0)]
            for j, e in enumerate(es, start=1):
                t[e] = R * (s[j - 1] - s[j])
        self.t_expr = t
        Jm = sp.Matrix([[ti for ti in t]] + [[sp.diff(ti, x) for ti in t] for x in self.xi])
        self.jac_expr = Jm   # rows: t ; d t / d xi_k
        self.t_f = sp.lambdify(self.xi, t, 'sympy')
        self.J_f = sp.lambdify(self.xi, Jm, 'sympy')

    def random_point(self, rng, den=97):
        """random rational interior point"""
        vals = {}
        for ci, (path, es) in enumerate(self.chains):
            a = len(es)
            if a >= 2:
                xs = sorted(Fr(rng.randint(1, den - 1), den) for _ in range(a - 1))[::-1]
                while len(set(xs)) < len(xs):
                    xs = sorted(Fr(rng.randint(1, den - 1), den) for _ in range(a - 1))[::-1]
                for k in range(1, a):
                    vals[(ci, k)] = xs[k - 1]
        while True:
            a, b = Fr(rng.randint(1, den - 1), den), Fr(rng.randint(1, den - 1), den)
            if a + b < 1:
                break
        pt = [vals[k] for k in sorted(self.sig)] + [a, b] + [Fr(rng.randint(1, den - 1), den) for _ in range(3)]
        return pt

    def t_at(self, pt):
        return [Fr(str(sp.nsimplify(v))) if not isinstance(v, Fr) else v for v in
                [sp.Rational(x) for x in self.t_f(*[sp.Rational(p.numerator, p.denominator) for p in pt])]]

    def g_ref(self, pt):
        P = [sp.Rational(p.numerator, p.denominator) for p in pt]
        t = [Fr(int(sp.numer(x)), int(sp.denom(x))) for x in self.t_f(*P)]
        Jm = self.J_f(*P)
        Jf = [[Fr(int(sp.numer(Jm[i, j])), int(sp.denom(Jm[i, j]))) for j in range(10)] for i in range(10)]
        Dn = brute.density(self.G, t, self.lam).rational()
        return Dn * abs(fdet(Jf)), t

    def currents_brute(self, t):
        return brute.currents(self.G, t, self.lam)

if __name__ == "__main__":
    rng = random.Random(5)
    for w in range(6):
        c = K4Case(w)
        pt = c.random_point(rng)
        gr, t = c.g_ref(pt)
        print(c.name, "Dc", c.Dc, "tri", c.tri, "xi", c.xi, "\n   g_ref =", float(gr))
