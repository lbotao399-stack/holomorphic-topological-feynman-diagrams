"""Three-loop graphs of Q_3 (hCS on C^{3|2}, S2): labelled graphs, reduction by the chain move
(S2 eq-chain), the bubble move (S2 eq-bubble) and the box (S2 eq-box-r, eq-box-currents), and an
independent brute-force layer (mesh currents, exact rationals).

Vertices 0..6 are the cubic vertices 1..7 of the text, vertex 7 is O (root).
Kernel (S2 eq-kernel-def):  I_Gamma = eps_T det(B_T) int_{omega_m} e^{-sum y.z}  wedge_{e not in T} d^3 y_e,
which equals the reduced kernel K_Gamma for n = 3 (CONVENTIONS.md Section 3).
"""
import random
from fractions import Fraction as Fr
from itertools import permutations

import sympy as sp

ROOT = 7
NV = 8

# ---------------------------------------------------------------------------------- graphs
# Each graph: edges (tail, head) in order; spokes [(edge index, sign)] for letters j = 1..4
# (z_e = sign * z_j); moves: list of ('chain', p, v, v', q) / ('bubble', p, v, v', q) / ('box', a, b, c),
# the chain/bubble edges being found from the current reduced graph.
# Text labels: vertex k -> k+1, ROOT -> O.

def _E(*pairs, root=ROOT):
    lab = lambda x: root if x == 'O' else x - 1
    return [(lab(a), lab(b)) for a, b in pairs]


GRAPHS = {
    # Gamma_1: triple box, path 1-2-3-4-5-6-7, spokes at 1, 3, 5, 7
    'G1': dict(edges=_E(('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), ('O', 5), (5, 6), (6, 7), (7, 'O')),
               spokes=[(0, 1), (3, 1), (6, 1), (9, -1)],
               moves=[('chain', 'O', 1, 2, 3), ('chain', 'O', 3, 4, 5), ('box', 5, 6, 7)]),
    # Gamma_2: box O-1-2-3 (spoke at 3), path 3-4-5, arms 5-6-O and 5-7-O
    'G2': dict(edges=_E(('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')),
               spokes=[(0, 1), (3, 1), (7, -1), (9, -1)],
               moves=[('chain', 'O', 1, 2, 3), ('bubble', 5, 6, 7, 'O'), ('box', 3, 4, 5)]),
    # Gamma_3: box O-1-2-3 (spoke at 3), edge 3-4, arms 4-5-O and 4-6-7-O
    'G3': dict(edges=_E(('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')),
               spokes=[(0, 1), (3, 1), (6, -1), (9, -1)],
               moves=[('chain', 'O', 1, 2, 3), ('bubble', 'O', 3, 5, 4), ('box', 4, 6, 7)]),
    # Gamma_4: arms O-1-3 and O-2-3, path 3-4-5, arms 5-6-O and 5-7-O
    'G4': dict(edges=_E(('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')),
               spokes=[(0, 1), (2, 1), (7, -1), (9, -1)],
               moves=[('bubble', 'O', 1, 2, 3), ('bubble', 5, 6, 7, 'O'), ('box', 3, 4, 5)]),
    # Gamma_5: arms O-1-3 and O-2-3, edge 3-4, arms 4-5-O and 4-6-7-O
    'G5': dict(edges=_E(('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')),
               spokes=[(0, 1), (2, 1), (6, -1), (9, -1)],
               moves=[('bubble', 'O', 1, 2, 3), ('bubble', 'O', 3, 5, 4), ('box', 4, 6, 7)]),
}

# two-loop graphs and the box of S2 (Sections 6, 7), for the validation of the conventions
GRAPHS_S2 = {
    'box': dict(edges=_E(('O', 1), (1, 2), (2, 3), (3, 'O'), root=3), spokes=[(0, 1), (3, -1)],
                moves=[('box', 1, 2, 3)], root=3),
    'A': dict(edges=_E(('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O'), root=5),
              spokes=[(0, 1), (3, 1), (6, -1)], moves=[('chain', 'O', 1, 2, 3), ('box', 3, 4, 5)], root=5),
    'B': dict(edges=_E(('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O'), root=5),
              spokes=[(0, 1), (2, 1), (6, -1)], moves=[('bubble', 'O', 1, 2, 3), ('box', 3, 4, 5)], root=5),
}
for g in GRAPHS.values():
    g['root'] = ROOT
for g in list(GRAPHS.values()) + list(GRAPHS_S2.values()):
    rt = g['root']
    g['moves'] = [(m[0],) + tuple(rt if x == 'O' else x - 1 for x in m[1:]) for m in g['moves']]
ALL = dict(GRAPHS, **GRAPHS_S2)


def degree(edges, v):
    return sum((a == v) + (b == v) for a, b in edges)


def outputs(name):
    """Vertices of degree two (one output letter each)."""
    E = ALL[name]['edges']
    return [v for v in range(ALL[name]['root']) if degree(E, v) == 2]


def tlabel(v):
    return 'O' if v == ROOT else str(v + 1)


# ---------------------------------------------------------------------------------- symbols
LAM = [[sp.Symbol('l%d_%d' % (v + 1, a + 1)) for a in range(3)] for v in range(ROOT)]
Z = [[sp.Symbol('z%d_%d' % (j + 1, a + 1)) for a in range(3)] for j in range(4)]
R = sp.symbols('r1:10')


def vec(x):
    return sp.Matrix(x)


def det3(a, b, c):
    return sp.Matrix.hstack(a, b, c).det()


# ---------------------------------------------------------------------------------- reduction
class Reduction:
    """Successive moves with the ratio conventions of S2:
       chain  (p,v,v',q):  T = t1+t2+t3, r1 = T/(T+t0), r2 = (t2+t3)/T, r3 = t3/T           (eq-chain)
       bubble (p,v,v',q):  T1 = t_pv+t_vq, T2 = t_pv'+t_v'q, r1 = T1/(T1+T2), r2 = t_vq/T1,
                           r3 = t_v'q/T2                                                     (eq-bubble)
       box    (a,b,c):     cycle O-a-b-c-O, r_j = (sum_{k>j} t_k)/T                          (eq-box-r)
    Ratios of move m are r_{3m-2}, r_{3m-1}, r_{3m}.  Gauge: T(box) = 1."""

    def __init__(self, name):
        g = ALL[name]
        self.name = name
        self.root = ROOT = g['root']
        self.edges = list(g['edges'])
        self.spokes = g['spokes']
        self.moves = g['moves']
        N = len(self.edges)
        cur = {e: ab for e, ab in enumerate(self.edges)}
        src = {v: vec(LAM[v]) for v in range(ROOT)}
        self.steps = []
        nid = N
        for mi, mv in enumerate(self.moves):
            r1, r2, r3 = R[3 * mi:3 * mi + 3]
            kind = mv[0]
            if kind in ('chain', 'bubble'):
                p, v, w, q = mv[1:]
                if kind == 'chain':
                    paths = [[p, v, w, q], [p, q]]
                else:
                    paths = [[p, v, q], [p, w, q]]
                chains = []
                for path in paths:
                    ch = []
                    for a, b in zip(path, path[1:]):
                        es = [e for e, ab in cur.items() if set(ab) == {a, b}]
                        assert len(es) == 1, (self.name, mi, a, b, es)
                        e = es[0]
                        ch.append((e, +1 if cur[e] == (a, b) else -1))
                    chains.append(ch)
                lv, lw = src[v], src[w]
                st = dict(kind=kind, p=p, v=v, w=w, q=q, chains=chains, new=nid, lv=lv, lw=lw,
                          r=(r1, r2, r3))
                for c in chains:
                    for e, _ in c:
                        del cur[e]
                cur[nid] = (p, q)
                if p != ROOT:
                    src[p] = src[p] + r2 * lv + r3 * lw
                if q != ROOT:
                    src[q] = src[q] + (1 - r2) * lv + (1 - r3) * lw
                del src[v], src[w]
                nid += 1
            else:
                a, b, c = mv[1:]
                path = [ROOT, a, b, c, ROOT]
                ch = []
                for x, y in zip(path, path[1:]):
                    es = [e for e, ab in cur.items() if set(ab) == {x, y}]
                    assert len(es) == 1
                    e = es[0]
                    ch.append((e, +1 if cur[e] == (x, y) else -1))
                st = dict(kind='box', verts=(a, b, c), chain=ch, lam=(src[a], src[b], src[c]), r=(r1, r2, r3))
                for e, _ in ch:
                    del cur[e]
            self.steps.append(st)
        assert not cur
        self._currents()
        self._tmap()

    def _currents(self):
        ROOT = self.root
        y = {}
        box = self.steps[-1]
        r1, r2, r3 = box['r']
        la, lb, lc = box['lam']
        y1 = -(r1 * la + r2 * lb + r3 * lc)
        along = [y1, y1 + la, y1 + la + lb, y1 + la + lb + lc]
        for (e, s), ya in zip(box['chain'], along):
            y[e] = s * ya
        box['y1'] = y1
        for st in reversed(self.steps[:-1]):
            r1, r2, r3 = st['r']
            Y = y[st['new']]                      # merged edge, oriented p -> q
            st['Y'] = Y
            lv, lw = st['lv'], st['lw']
            if st['kind'] == 'chain':
                ya = (1 - r1) * Y - r2 * lv - r3 * lw
                along = [[ya, ya + lv, ya + lv + lw], [r1 * Y]]
            else:
                ya = (1 - r1) * Y - r2 * lv
                yc = r1 * Y - r3 * lw
                along = [[ya, ya + lv], [yc, yc + lw]]
            for ch, al in zip(st['chains'], along):
                for (e, s), yy in zip(ch, al):
                    y[e] = s * yy
            st['ya'] = ya
        self.y = [sp.expand(y[e]) for e in range(len(self.edges))]

    def _tmap(self):
        """Schwinger times of the original edges as rational functions of the ratios (gauge T_box = 1)."""
        t = {}
        box = self.steps[-1]
        r1, r2, r3 = box['r']
        for (e, _), te in zip(box['chain'], [1 - r1, r1 - r2, r2 - r3, r3]):
            t[e] = te
        for st in reversed(self.steps[:-1]):
            r1, r2, r3 = st['r']
            tred = t[st['new']]
            if st['kind'] == 'chain':
                T = tred / (1 - r1)
                t0 = tred / r1
                (a, _), (b, _), (c, _) = st['chains'][0]
                (d, _), = st['chains'][1]
                t[a], t[b], t[c], t[d] = T * (1 - r2), T * (r2 - r3), T * r3, t0
            else:
                T1 = tred / (1 - r1)
                T2 = tred / r1
                (a, _), (b, _) = st['chains'][0]
                (c, _), (d, _) = st['chains'][1]
                t[a], t[b], t[c], t[d] = T1 * (1 - r2), T1 * r2, T2 * (1 - r3), T2 * r3
        self.t = [sp.simplify(t[e]) for e in range(len(self.edges))]

    # local factors (S2 eq-chain-jac, eq-bubble-jac, eq-box-jac) without signs
    def factors(self):
        out = []
        for st in self.steps[:-1]:
            r1 = st['r'][0]
            d = det3(st['Y'], st['lv'], st['lw'])
            pre = r1 ** 2 if st['kind'] == 'chain' else r1 * (1 - r1)
            out.append(sp.expand(pre * d))
        la, lb, lc = self.steps[-1]['lam']
        out.append(sp.expand(det3(la, lb, lc)))
        return out

    def domain(self):
        """[(kind, ratios)]: 'cube' r in (0,1); 'ord2' 1 > a > b > 0; 'ord3' 1 > a > b > c > 0."""
        dom = []
        for st in self.steps:
            r1, r2, r3 = st['r']
            if st['kind'] == 'chain':
                dom += [('cube', (r1,)), ('ord2', (r2, r3))]
            elif st['kind'] == 'bubble':
                dom += [('cube', (r1,)), ('cube', (r2,)), ('cube', (r3,))]
            else:
                dom += [('ord3', (r1, r2, r3))]
        return dom

    def cotree(self):
        """Cotree used by the moves: for each move the first original (not merged) edge of its first
        chain; for the box the first original edge of the cycle.  Modulo the differentials of the
        ratios of earlier moves, y of any chain edge equals y_a + (reduced sources)."""
        N = len(self.edges)
        ct = []
        for st in self.steps[:-1]:
            ct.append(next(e for e, _ in st['chains'][0] if e < N))
        ct.append(next(e for e, _ in self.steps[-1]['chain'] if e < N))
        return ct

    def displacement(self):
        """X_v (v over non-root vertices) with -sum_e y_e.z_e = sum_v lam_v . X_v; z_e = s z_j on spokes."""
        ROOT = self.root
        zs = {e: s * vec(Z[j]) for j, (e, s) in enumerate(self.spokes)}
        expo = 0
        for e, zz in zs.items():
            expo -= (self.y[e].T * zz)[0, 0]
        expo = sp.expand(expo)
        X = []
        for v in range(ROOT):
            Xv = []
            for a in range(3):
                Xv.append(sp.expand(sp.diff(expo, LAM[v][a])))
            X.append(Xv)
        # consistency: X_v^a must not depend on a beyond z_j^a
        return X, expo


def integrate_domain(poly, dom):
    """Exact integral of a polynomial in R over the product domain."""
    P = sp.Poly(sp.expand(poly), *R)
    idx = {r: i for i, r in enumerate(R)}
    tot = 0
    for mon, c in P.terms():
        w = sp.Integer(1)
        for kind, rs in dom:
            ks = [mon[idx[r]] for r in rs]
            if kind == 'cube':
                w *= sp.Rational(1, ks[0] + 1)
            elif kind == 'ord2':
                a, b = ks
                w *= sp.Rational(1, (b + 1) * (a + b + 2))
            else:
                a, b, c3 = ks
                w *= sp.Rational(1, (c3 + 1) * (b + c3 + 2) * (a + b + c3 + 3))
        tot += c * w
    return sp.expand(tot)


def random_point(dom, rng):
    pt = {}
    for kind, rs in dom:
        if kind == 'cube':
            pt[rs[0]] = Fr(rng.randint(1, 98), 99)
        else:
            vals = sorted(rng.sample(range(1, 997), len(rs)), reverse=True)
            for r, x in zip(rs, vals):
                pt[r] = Fr(x, 997)
    return pt


# ---------------------------------------------------------------------------------- brute force
def fdet(M):
    M = [list(map(Fr, row)) for row in M]
    n = len(M)
    d = Fr(1)
    for i in range(n):
        piv = next((k for k in range(i, n) if M[k][i] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            d = -d
        d *= M[i][i]
        for k in range(i + 1, n):
            f = M[k][i] / M[i][i]
            if f:
                for j in range(i, n):
                    M[k][j] -= f * M[i][j]
    return d


def finv(M):
    n = len(M)
    A = [list(map(Fr, row)) + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for i in range(n):
        piv = next(k for k in range(i, n) if A[k][i] != 0)
        A[i], A[piv] = A[piv], A[i]
        p = A[i][i]
        A[i] = [x / p for x in A[i]]
        for k in range(n):
            if k != i and A[k][i] != 0:
                f = A[k][i]
                A[k] = [x - f * y for x, y in zip(A[k], A[i])]
    return [row[n:] for row in A]


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        while seq[i] != i:
            j = seq[i]
            seq[i], seq[j] = seq[j], seq[i]
            s = -s
    return s


class Brute:
    """Kirchhoff network at exact rational Schwinger times, mesh formulation:
    y = y0 + C^T l, y0 = tree solution of the current law, (C t C^T) l = -C t y0."""

    def __init__(self, edges, cotree, root=ROOT):
        self.root = ROOT = root
        NV = root + 1
        self.edges = edges
        self.N = len(edges)
        self.cot = list(cotree)
        self.tree = [e for e in range(self.N) if e not in self.cot]
        g = ROOT
        # grounded incidence B_T (rows: vertices 0..g-1, columns: tree edges in increasing order)
        BT = [[0] * len(self.tree) for _ in range(g)]
        for j, e in enumerate(self.tree):
            a, b = edges[e]
            if a < g:
                BT[a][j] += 1
            if b < g:
                BT[b][j] -= 1
        self.detBT = fdet(BT)
        assert self.detBT in (1, -1), 'not a spanning tree'
        order = self.tree + self.cot
        self.epsT = perm_sign([order.index(i) for i in range(self.N)])
        self.BTinv = finv(BT)          # tree currents: y_T = BT^{-1} lam
        # fundamental cycles: C[r][c_r] = 1, cycle along c_r back through the tree
        adj = {v: [] for v in range(NV)}
        for e in self.tree:
            a, b = edges[e]
            adj[a].append((b, e, 1))
            adj[b].append((a, e, -1))
        C = []
        for c in self.cot:
            a, b = edges[c]
            prev = {b: None}
            stack = [b]
            while stack:
                x = stack.pop()
                for (u, e, s) in adj[x]:
                    if u not in prev:
                        prev[u] = (x, e, s)
                        stack.append(u)
            row = [0] * self.N
            row[c] = 1
            x = a
            while prev[x] is not None:
                u, e, s = prev[x]
                row[e] += s
                x = u
            C.append(row)
        self.C = C
        # check: C B = 0
        for row in C:
            for v in range(g):
                assert sum(row[e] * ((edges[e][0] == v) - (edges[e][1] == v)) for e in range(self.N)) == 0

    def currents(self, t, lam):
        """lam: list over the non-root vertices of 3-vectors (Fractions)."""
        ROOT = self.root
        N, n = self.N, len(self.cot)
        y0 = [[Fr(0)] * 3 for _ in range(N)]
        for j, e in enumerate(self.tree):
            for a in range(3):
                y0[e][a] = sum((self.BTinv[j][v] * lam[v][a] for v in range(ROOT)), Fr(0))
        M = [[sum((Fr(self.C[r][e] * self.C[s][e]) * t[e] for e in range(N)), Fr(0)) for s in range(n)] for r in range(n)]
        Minv = finv(M)
        rhs = [[-sum((self.C[r][e] * t[e] * y0[e][a] for e in range(N)), Fr(0)) for a in range(3)] for r in range(n)]
        l = [[sum((Minv[r][s] * rhs[s][a] for s in range(n)), Fr(0)) for a in range(3)] for r in range(n)]
        y = [[y0[e][a] + sum((self.C[r][e] * l[r][a] for r in range(n)), Fr(0)) for a in range(3)] for e in range(N)]
        return y, Minv

    def dy_dt(self, t, lam):
        """d y_e^a / d t_j = sum_r C_re dl_r^a/dt_j,  dl/dt_j = -M^{-1} C_{:,j} y_j."""
        y, Minv = self.currents(t, lam)
        N, n = self.N, len(self.cot)
        D = [[[Fr(0)] * N for _ in range(3)] for _ in range(N)]
        for j in range(N):
            for a in range(3):
                dl = [-sum((Minv[r][s] * self.C[s][j] for s in range(n)), Fr(0)) * y[j][a] for r in range(n)]
                for e in range(N):
                    D[e][a][j] = sum((self.C[r][e] * dl[r] for r in range(n)), Fr(0))
        return y, D

    def check_laws(self, t, lam, y):
        for v in range(self.root):
            for a in range(3):
                s = sum(((self.edges[e][0] == v) - (self.edges[e][1] == v)) * y[e][a] for e in range(self.N))
                assert s == lam[v][a]
        for row in self.C:
            for a in range(3):
                assert sum(row[e] * t[e] * y[e][a] for e in range(self.N)) == 0


def lam_sub(lam):
    return {LAM[v][a]: lam[v][a] for v in range(len(lam)) for a in range(3)}
