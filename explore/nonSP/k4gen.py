"""K4 subdivision with arbitrary sources (H components) at any vertices: chain elimination + sourced Y-Delta chart
(same formulas as k4formula.py, with initial K4 sources mu_v = lam_v).  Used for the C^2 three-loop non-SP graph
[(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,4)] (K4 on {0,1,2,4} with the edge 01 subdivided by 3)."""
import sys, os, random, json
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk.graph import Graph
from htk import brute
from htk.linalg import det as fdet
from graphs import skeleton_chains

class K4Gen:
    def __init__(self, name, edges, nv, lam, H, Dc, tri):
        self.name, self.edges, self.nv, self.H = name, edges, nv, H
        self.G = Graph(nv, edges)
        self.lam = lam                      # dict v -> list of H Fractions (sum zero over all v)
        self.cubic, self.chains = skeleton_chains(edges, nv)
        self.Dc, self.tri = Dc, tri
        self.gA, self.gB, self.rA, self.rB, self.rC = sp.symbols('gA gB rA rB rC')
        self.sig = {}
        for ci, (path, es) in enumerate(self.chains):
            for k in range(1, len(es)):
                self.sig[(ci, k)] = sp.Symbol('s%d_%d' % (ci, k))
        self.xi = [self.sig[k] for k in sorted(self.sig)] + [self.gA, self.gB, self.rA, self.rB, self.rC]
        A, B, C = tri
        g = {A: self.gA, B: self.gB, C: 1 - self.gA - self.gB}
        r = {A: self.rA, B: self.rB, C: self.rC}
        cond = {}
        for v in tri:
            cond[frozenset((Dc, v))] = g[v]
            w, x = [u for u in tri if u != v]
            cond[frozenset((w, x))] = g[w] * g[x] * r[v] / (1 - r[v])
        t = [None] * len(edges)
        for ci, (path, es) in enumerate(self.chains):
            a = len(es); R = 1 / cond[frozenset((path[0], path[-1]))]
            s = [sp.Integer(1)] + [self.sig[(ci, k)] for k in range(1, a)] + [sp.Integer(0)]
            for j, e in enumerate(es, start=1):
                t[e] = R * (s[j - 1] - s[j])
        self.t_expr = t
        Jm = sp.Matrix([t] + [[sp.diff(ti, x) for ti in t] for x in self.xi])
        self.t_f = sp.lambdify(self.xi, t, 'sympy'); self.J_f = sp.lambdify(self.xi, Jm, 'sympy')
        # currents
        def _s(x):
            return sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.sympify(x)
        lamv = {v: sp.Matrix([_s(x) for x in lam.get(v, [0] * H)]) for v in range(nv)}
        mu = {v: lamv[v] for v in self.cubic}
        mhat, Qs = {}, {}
        for ci, (path, es) in enumerate(self.chains):
            p, q = path[0], path[-1]; a = len(es)
            mh = sp.zeros(H, 1)
            for k in range(1, a):
                sk = self.sig[(ci, k)]
                mu[p] = mu[p] + sk * lamv[path[k]]; mu[q] = mu[q] + (1 - sk) * lamv[path[k]]
                mh += sk * lamv[path[k]]
            mhat[ci] = mh
            Qs[ci] = [sum((lamv[path[k]] for k in range(1, j)), sp.zeros(H, 1)) for j in range(1, a + 1)]
        pv = {v: (1 - r[v]) * g[v] for v in tri}
        Q = pv[A] + pv[B] + pv[C]; self.Q = Q
        mup = {v: mu[v] + g[v] * mu[Dc] for v in tri}
        Ytri = lambda X, Y: (pv[Y] * mup[X] - pv[X] * mup[Y]) / Q
        third = lambda X, Y: [u for u in tri if u not in (X, Y)][0]
        def YK4(p, q):
            if p == Dc:
                return g[q] * mu[Dc] - sum(((1 - r[third(q, w)]) * Ytri(q, w) for w in tri if w != q), sp.zeros(H, 1))
            if q == Dc:
                return -YK4(q, p)
            return r[third(p, q)] * Ytri(p, q)
        y = [None] * len(edges)
        for ci, (path, es) in enumerate(self.chains):
            Y = YK4(path[0], path[-1])
            for j, e in enumerate(es, start=1):
                val = Y - mhat[ci] + Qs[ci][j - 1]
                y[e] = val if edges[e] == (path[j - 1], path[j]) else -val
        self.y = y
        self.extra = sorted(set().union(*[sp.sympify(x).free_symbols for v in lam.values() for x in v if not isinstance(x, Fr)]), key=str)
        self.y_f = sp.lambdify(self.xi + self.extra, [list(v) for v in y], 'sympy')

    def lam_list(self, subs=None):
        out = []
        for v in range(self.nv):
            row = []
            for x in self.lam.get(v, [Fr(0)] * self.H):
                if isinstance(x, Fr):
                    row.append(x)
                else:
                    val = sp.Rational(sp.sympify(x).subs(subs or {}))
                    row.append(Fr(int(sp.numer(val)), int(sp.denom(val))))
            out.append(row)
        return out

    def random_point(self, rng, den=97):
        vals = {}
        for ci, (path, es) in enumerate(self.chains):
            a = len(es)
            if a >= 2:
                xs = sorted(set(Fr(rng.randint(1, den - 1), den) for _ in range(a - 1)))[::-1]
                while len(xs) < a - 1:
                    xs = sorted(set(Fr(rng.randint(1, den - 1), den) for _ in range(a - 1)))[::-1]
                for k in range(1, a):
                    vals[(ci, k)] = xs[k - 1]
        while True:
            a, b = Fr(rng.randint(1, den - 1), den), Fr(rng.randint(1, den - 1), den)
            if a + b < 1: break
        return [vals[k] for k in sorted(self.sig)] + [a, b] + [Fr(rng.randint(1, den - 1), den) for _ in range(3)]

    @staticmethod
    def _fr(x):
        x = sp.Rational(x); return Fr(int(sp.numer(x)), int(sp.denom(x)))

    def check(self, pt, subs=None):
        P = [sp.Rational(p.numerator, p.denominator) for p in pt]
        t = [self._fr(x) for x in self.t_f(*P)]
        yb = brute.currents(self.G, t, self.lam_list(subs))
        yf = self.y_f(*(P + [(subs or {})[sy] for sy in self.extra]))
        okc = all(self._fr(yf[e][i]) == yb[e][i] for e in range(len(self.edges)) for i in range(self.H))
        Jm = self.J_f(*P); N = len(self.edges)
        Jf = [[self._fr(Jm[i, j]) for j in range(N)] for i in range(N)]
        g = brute.density(self.G, t, self.lam_list(subs)).rational() * abs(fdet(Jf))
        return okc, g

def write_wl(case, base, npts=3, seed=5, symsrc=None):
    ren = {s: sp.Symbol(str(s).replace('_', 'x')) for s in case.xi}
    cot = case.G.cotree()
    ys = [case.y[e][i].subs(ren) for e in cot for i in range(case.H)]
    xi = [ren[s] for s in case.xi]
    rng = random.Random(seed)
    pts, gref, oks = [], [], []
    for _ in range(npts):
        pt = case.random_point(rng)
        subs = {sy: sp.Rational(rng.randint(-9, 9), 7) for sy in (symsrc or [])}
        ok, g = case.check(pt, subs)
        pts.append((pt, subs)); gref.append(g); oks.append(ok)
    m = len(ys)
    L = ['ys = {%s};' % ', '.join(sp.mathematica_code(e) for e in ys),
         'Qp = %s;' % sp.mathematica_code(case.Q.subs(ren)),
         'vars = {%s};' % ', '.join(str(v) for v in xi),
         'Ns = Expand[Together[ys*Qp]];',
         'If[Union[Denominator /@ Ns] =!= {1}, Print["NOT POLYNOMIAL"]; Exit[1]];',
         'bord = Prepend[Table[Prepend[D[Ns[[j]], #] & /@ vars, Ns[[j]]], {j, %d}], Prepend[D[Qp, #] & /@ vars, Qp]];' % m,
         'detB = Expand[Det[bord]];', 'Print["detB terms ", Length[detB]];']
    for k, (pt, subs) in enumerate(pts):
        rule = ', '.join('%s -> %d/%d' % (v, p.numerator, p.denominator) for v, p in zip(xi, pt))
        for sy, val in subs.items():
            rule += ', %s -> %s' % (sy, sp.mathematica_code(val))
        L.append('Print["gsym%d ", InputForm[(detB/Qp^%d) /. {%s}]];' % (k, m + 1, rule))
    L.append('P = detB;')
    for ci, (path, es) in enumerate(case.chains):
        d = [str(ren[case.sig[(ci, k)]]) for k in range(1, len(es))]
        for k in range(len(d) - 1, -1, -1):
            L.append('P = Expand[Integrate[P, {%s, 0, %s}]];' % (d[k], '1' if k == 0 else d[k - 1]))
    L.append('Print["P terms ", Length[P]];')
    L.append('Put[P, "%s_P.m"]; Put[Qp, "%s_Q.m"];' % (base, base))
    open('wl/' + base.split('/')[-1] + '.wl', 'w').write('\n'.join(L))
    json.dump({'name': case.name, 'Dc': case.Dc, 'tri': case.tri, 'gref': [str(g) for g in gref], 'currents_ok': oks,
               'apow': m + 1}, open(base + '_ref.json', 'w'))
    print(case.name, 'Dc', case.Dc, 'currents == brute:', oks, 'apow', m + 1)

N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]
if __name__ == "__main__":
    # sources at 3 (deg 2), 2, 4 (root): lam3 = e1, lam2 = e2, lam4 = -(e1+e2);  K = c <lam3, lam2>^3
    lam = {3: [Fr(1), Fr(0)], 2: [Fr(0), Fr(1)], 4: [Fr(-1), Fr(-1)]}
    for Dc, tri, tag in [(4, [0, 1, 2], 'D4'), (0, [1, 2, 4], 'D0')]:
        c = K4Gen('C2n3[3;2,4]', N3, 5, lam, 2, Dc, tri)
        write_wl(c, 'out/c2n3_%s' % tag)
    # sources at 3, 0, 2: K4 vertex 4 source-free (rational expected)
    lam = {3: [Fr(1), Fr(0)], 0: [Fr(0), Fr(1)], 2: [Fr(-1), Fr(-1)]}
    c = K4Gen('C2n3[3;0,2]', N3, 5, lam, 2, 4, [0, 1, 2])
    write_wl(c, 'out/c2n3b_D4')
    # four sources: lam3 = e1, lam0 = e2, lam2 = (SA, SB), lam4 = -(sum): all K4 vertices sourced, no symmetry
    SA, SB = sp.symbols('SA SB')
    lam = {3: [Fr(1), Fr(0)], 0: [Fr(0), Fr(1)], 2: [SA, SB], 4: [-1 - SA, -1 - SB]}
    for Dc, tri, tag in [(4, [0, 1, 2], 'D4'), (1, [0, 2, 4], 'D1')]:
        c = K4Gen('C2n3[3;0,2,4]', N3, 5, lam, 2, Dc, tri)
        write_wl(c, 'out/c2n3c_%s' % tag, symsrc=[SA, SB])
