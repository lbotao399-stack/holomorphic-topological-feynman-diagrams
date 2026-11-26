"""C^2 (H=2, T=0), n=4: prism and K33 with sources at four vertices, two sourceless vertices.
Reduction (my construction; all steps are bijective changes of Schwinger variables):
  two source-free star-triangle moves at the sourceless vertices (merges -> fractions f, no merge -> new edge),
  then the sourced star-triangle chart on the remaining K4 (as in the C^3 K4 case, no chains):
     K4 vertex Dc, triangle (A,B,C); conductances g_v (Dc-v), h_v = g_w g_x r_v/(1-r_v);
     p_v = (1-r_v) g_v, Q = sum p_v, mu'_v = lam_v + g_v lam_Dc, Y_XY = (p_Y mu'_X - p_X mu'_Y)/Q,
     triangle edge X->Y: r_Z Y_XY;  star edge Dc->v: g_v lam_Dc - sum_w (1-r_{Z(vw)}) Y_vw.
  Undo of a source-free move at v with neighbours u_i (upper level: triangle edges T_ij, currents I_ij):
     merged edge: lower edge T/f with current f I, branch T/(1-f) with current (1-f) I; else branch T, current I;
     star a_{u_i} = tau_ij tau_ik / (tau_12+tau_13+tau_23),  y_{v->u_i} = sum_{j != i} J^br_{u_j -> u_i}.
Sources: lam_s0 = e1, lam_s1 = e2, lam_s2 = (SA, SB), lam_s3 = -(lam_s0+lam_s1+lam_s2), zero at the two sourceless vertices.
"""
import sys, os, random
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))
from htk.graph import Graph
from htk import brute
from htk.linalg import det as fdet
from graphs import PRISM, K33

def fe(a, b):
    return frozenset((a, b))

class C2Case:
    def __init__(self, name, edges, src, moves, Dc, tri, three=False):
        """src: 4 sourced vertices (s0,s1,s2,s3); moves: list of (v, (u1,u2,u3)) in the order applied
        to the original graph; Dc, tri: chart of the final K4."""
        self.name, self.edges = name, edges
        self.nv = 6
        self.G = Graph(6, edges)
        self.src = src
        self.SA, self.SB = sp.symbols('SA SB')
        lam = {src[0]: sp.Matrix([1, 0]), src[1]: sp.Matrix([0, 1]), src[2]: sp.Matrix([self.SA, self.SB])}
        if three:
            # three sources: lam_s0 = e1, lam_s1 = e2, lam_s2 = -(e1+e2) (SA, SB unused)
            lam = {src[0]: sp.Matrix([1, 0]), src[1]: sp.Matrix([0, 1]), src[2]: sp.Matrix([-1, -1])}
        else:
            lam[src[3]] = -(lam[src[0]] + lam[src[1]] + lam[src[2]])
        self.lamsym = {v: lam.get(v, sp.zeros(2, 1)) for v in range(6)}
        self.moves = moves
        self.Dc, self.tri = Dc, tri
        # replay the moves on edge sets to know which triangle edges merge
        E = set(fe(a, b) for a, b in edges)
        self.levels = []
        for v, us in moves:
            star = [fe(v, u) for u in us]
            assert all(s in E for s in star) and sum(1 for e in E if v in e) == 3
            pairs = [fe(us[0], us[1]), fe(us[0], us[2]), fe(us[1], us[2])]
            merged = {p: (p in E) for p in pairs}
            self.levels.append((v, us, merged))
            E = (E - set(star)) | set(pairs)
        assert E == set(fe(a, b) for a in [Dc] + tri for b in [Dc] + tri if a < b), E
        self.gA, self.gB, self.rA, self.rB, self.rC = sp.symbols('gA gB rA rB rC')
        self.fr = {}
        for li, (v, us, merged) in enumerate(self.levels):
            for p, m in merged.items():
                if m:
                    a, b = sorted(p)
                    self.fr[(li, a, b)] = sp.Symbol('f%d_%d%d' % (li, a, b))
        self.xi = [self.fr[k] for k in sorted(self.fr)] + [self.gA, self.gB, self.rA, self.rB, self.rC]
        self._build()

    def _build(self):
        A, B, C = self.tri
        Dc = self.Dc
        g = {A: self.gA, B: self.gB, C: 1 - self.gA - self.gB}
        r = {A: self.rA, B: self.rB, C: self.rC}
        lam = self.lamsym
        pv = {v: (1 - r[v]) * g[v] for v in self.tri}
        Q = pv[A] + pv[B] + pv[C]
        self.Q = Q
        mup = {v: lam[v] + g[v] * lam[Dc] for v in self.tri}
        def Ytri(X, Y):
            return (pv[Y] * mup[X] - pv[X] * mup[Y]) / Q
        def third(X, Y):
            return [u for u in self.tri if u not in (X, Y)][0]
        # top level: resistance and current (oriented min -> max) per edge
        net = {}
        for v in self.tri:
            a, b = sorted((Dc, v))
            Rv = 1 / g[v]
            yDv = g[v] * lam[Dc] - sum(((1 - r[third(v, w)]) * Ytri(v, w) for w in self.tri if w != v), sp.zeros(2, 1))
            net[fe(Dc, v)] = (Rv, yDv if a == Dc else -yDv)
            w, x = [u for u in self.tri if u != v]
            a, b = sorted((w, x))
            net[fe(w, x)] = ((1 - r[v]) / (r[v] * g[w] * g[x]), r[v] * Ytri(a, b))
        # undo the source-free moves in reverse order
        for li in reversed(range(len(self.levels))):
            v, us, merged = self.levels[li]
            tau, Jbr = {}, {}
            new = dict(net)
            for p, m in merged.items():
                a, b = sorted(p)
                T, I = net[p]
                del new[p]
                if m:
                    f = self.fr[(li, a, b)]
                    new[p] = (T / f, f * I)
                    tau[p], Jbr[p] = T / (1 - f), (1 - f) * I
                else:
                    tau[p], Jbr[p] = T, I
            S = sum(tau.values())
            for i, u in enumerate(us):
                o = [w for w in us if w != u]
                a_u = tau[fe(u, o[0])] * tau[fe(u, o[1])] / S
                yvu = sp.zeros(2, 1)
                for w in o:
                    J = Jbr[fe(u, w)]                 # oriented min -> max
                    yvu += J if w < u else -J         # current w -> u
                a, b = sorted((v, u))
                new[fe(v, u)] = (a_u, yvu if a == v else -yvu)
            net = new
        self.net = net
        assert set(net) == set(fe(a, b) for a, b in self.edges)
        self.t_expr = [net[fe(a, b)][0] for a, b in self.edges]
        self.y_expr = [net[fe(a, b)][1] if a < b else -net[fe(a, b)][1] for a, b in self.edges]
        self.allsyms = self.xi + [self.SA, self.SB]
        self.t_f = sp.lambdify(self.xi, self.t_expr, 'sympy')
        self.y_f = sp.lambdify(self.allsyms, [list(v) for v in self.y_expr], 'sympy')
        Jm = sp.Matrix([self.t_expr] + [[sp.diff(ti, x) for ti in self.t_expr] for x in self.xi])
        self.J_f = sp.lambdify(self.xi, Jm, 'sympy')

    def lam_num(self, SA, SB):
        out = []
        for v in range(6):
            col = self.lamsym[v].subs({self.SA: SA, self.SB: SB})
            out.append([Fr(int(sp.numer(c)), int(sp.denom(c))) for c in col])
        return out

    def random_point(self, rng, den=89):
        pt = [Fr(rng.randint(1, den - 1), den) for _ in self.fr]
        while True:
            a, b = Fr(rng.randint(1, den - 1), den), Fr(rng.randint(1, den - 1), den)
            if a + b < 1:
                break
        return pt + [a, b] + [Fr(rng.randint(1, den - 1), den) for _ in range(3)]

    @staticmethod
    def _fr(x):
        x = sp.Rational(x)
        return Fr(int(sp.numer(x)), int(sp.denom(x)))

    def check_currents(self, rng, SA, SB):
        pt = self.random_point(rng)
        P = [sp.Rational(p.numerator, p.denominator) for p in pt]
        t = [self._fr(x) for x in self.t_f(*P)]
        lam = self.lam_num(SA, SB)
        yb = brute.currents(self.G, t, lam)
        yf = self.y_f(*(P + [sp.Rational(SA), sp.Rational(SB)]))
        return all(self._fr(yf[e][i]) == yb[e][i] for e in range(9) for i in range(2))

    def g_ref(self, pt, SA, SB):
        P = [sp.Rational(p.numerator, p.denominator) for p in pt]
        t = [self._fr(x) for x in self.t_f(*P)]
        Jm = self.J_f(*P)
        Jf = [[self._fr(Jm[i, j]) for j in range(9)] for i in range(9)]
        return brute.density(self.G, t, self.lam_num(SA, SB)).rational() * abs(fdet(Jf))

def prism_case(Dc_choice=0):
    # prism: triangles (0,2,3), (1,4,5); sources at 0,1,2,3; moves at 4 (nbrs 1,3,5) then 5 (nbrs 1,2,3)
    K4v = [0, 1, 2, 3]
    Dc = K4v[Dc_choice]
    return C2Case("prism", PRISM, (0, 1, 2, 3), [(4, (1, 3, 5)), (5, (1, 2, 3))], Dc, [v for v in K4v if v != Dc])

def k33_case(Dc_choice=0):
    # K33: sides {0,2,4},{1,3,5}; sources at 0,1,3,5; moves at 2 (nbrs 1,3,5) then 4 (nbrs 1,3,5)
    K4v = [0, 1, 3, 5]
    Dc = K4v[Dc_choice]
    return C2Case("K33", K33, (0, 1, 3, 5), [(2, (1, 3, 5)), (4, (1, 3, 5))], Dc, [v for v in K4v if v != Dc])

def prism3_case(Dc_choice=3):
    # prism, sources at 0,1,2 (three); moves at 3 (nbrs 0,2,4) then 5 (nbrs 1,2,4); final K4 {0,1,2,4}
    K4v = [0, 1, 2, 4]
    Dc = K4v[Dc_choice]
    return C2Case("prism3", PRISM, (0, 1, 2), [(3, (0, 2, 4)), (5, (1, 2, 4))], Dc, [v for v in K4v if v != Dc], three=True)

def k33_3_case(Dc_choice=3):
    # K33, sources at 0,1,2 (three); moves at 3 (nbrs 0,2,4) then 5 (nbrs 0,2,4); final K4 {0,1,2,4}
    K4v = [0, 1, 2, 4]
    Dc = K4v[Dc_choice]
    return C2Case("K33_3", K33, (0, 1, 2), [(3, (0, 2, 4)), (5, (0, 2, 4))], Dc, [v for v in K4v if v != Dc], three=True)

if __name__ == "__main__":
    rng = random.Random(3)
    for mk in (prism3_case, k33_3_case, prism_case, k33_case):
        for dci in range(4):
            c = mk(dci)
            oks = [c.check_currents(rng, Fr(rng.randint(-9, 9), 7), Fr(rng.randint(-9, 9), 5)) for _ in range(3)]
            print(c.name, "Dc", c.Dc, "tri", c.tri, "moves", c.levels, "currents == brute:", oks)
