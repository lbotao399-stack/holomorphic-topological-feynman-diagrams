"""Independent exact layer for the Q_3 graphs (own code; nothing imported from the explorer).

Kernel of S2 eq-kernel-def:
  I_Gamma = eps_T det(B_T) int_{omega_m} e^{-sum y.z} wedge_{e not in T} d^3 y_e,
  omega_m = iota_E dt_1...dt_m,  (B_T)_{ve} = +1 (v = tail e), -1 (v = head e), v != O, e in T,
  eps_T = sign of (tree edges increasing, cotree edges increasing).
Currents: nodal analysis, L = B^T diag(1/t) B, phi = L^{-1} lambda, y = diag(1/t) B phi
(so that sum_{tail=v} y - sum_{head=v} y = lambda_v and t_e y_e = phi_tail - phi_head, S2 eq-KK).
Density in ratio coordinates r (any chart t(r)):
  g(r) = eps_T det B_T * det[d y_{c,a} / d r_k]_{(c,a),k} * sgn det[t ; dt/dr],
rows (c,a): cotree edges increasing, a = 1,2,3.
Chart t(r): backward from the final box (T_box = 1) through the moves of S2 eq-chain, eq-bubble,
eq-box-r, with forward-mode exact differentiation (Dual numbers over Fractions).
"""
import random
from fractions import Fraction as Fr

GRAPHS = {
    'G1': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), ('O', 5), (5, 6), (6, 7), (7, 'O')],
               moves=[('chain', 'O', 1, 2, 3), ('chain', 'O', 3, 4, 5), ('box', 'O', 5, 6, 7)]),
    'G2': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
               moves=[('chain', 'O', 1, 2, 3), ('bubble', 5, 6, 7, 'O'), ('box', 'O', 3, 4, 5)]),
    'G3': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
               moves=[('chain', 'O', 1, 2, 3), ('bubble', 'O', 3, 5, 4), ('box', 'O', 4, 6, 7)]),
    'G4': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 'O'), (5, 7), (7, 'O')],
               moves=[('bubble', 'O', 1, 2, 3), ('bubble', 5, 6, 7, 'O'), ('box', 'O', 3, 4, 5)]),
    'G5': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O'), (4, 6), (6, 7), (7, 'O')],
               moves=[('bubble', 'O', 1, 2, 3), ('bubble', 'O', 3, 5, 4), ('box', 'O', 4, 6, 7)]),
    # S2 graphs (Sections 6, 7) for validation of the layer
    'box': dict(edges=[('O', 1), (1, 2), (2, 3), (3, 'O')], moves=[('box', 'O', 1, 2, 3)]),
    'A': dict(edges=[('O', 1), (1, 2), (2, 3), ('O', 3), (3, 4), (4, 5), (5, 'O')],
              moves=[('chain', 'O', 1, 2, 3), ('box', 'O', 3, 4, 5)]),
    'B': dict(edges=[('O', 1), (1, 3), ('O', 2), (2, 3), (3, 4), (4, 5), (5, 'O')],
              moves=[('bubble', 'O', 1, 2, 3), ('box', 'O', 3, 4, 5)]),
}
for _g in GRAPHS.values():
    _g['verts'] = sorted({x for e in _g['edges'] for x in e if x != 'O'})
    # spokes in edge order: letter j on the j-th edge incident to O; z_e = +z_j (O->v), -z_j (v->O)
    _g['spokes'] = [(i, 1 if e[0] == 'O' else -1) for i, e in enumerate(_g['edges']) if 'O' in e]
    _g['outputs'] = [v for v in _g['verts'] if sum(v in e for e in _g['edges']) == 2]


# ----------------------------------------------------------------------------- linear algebra
def det(M):
    M = [[Fr(x) for x in r] for r in M]
    n = len(M)
    d = Fr(1)
    for i in range(n):
        p = next((k for k in range(i, n) if M[k][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        inv = 1 / M[i][i]
        for k in range(i + 1, n):
            f = M[k][i] * inv
            if f:
                for j in range(i, n):
                    M[k][j] -= f * M[i][j]
    return d


def solve(A, Bcols):
    """A X = B; A n x n, B n x k (list of rows)."""
    n = len(A)
    k = len(Bcols[0])
    M = [list(A[i]) + list(Bcols[i]) for i in range(n)]
    for i in range(n):
        p = next(kk for kk in range(i, n) if M[kk][i] != 0)
        M[i], M[p] = M[p], M[i]
        inv = 1 / M[i][i]
        M[i] = [x * inv for x in M[i]]
        for kk in range(n):
            if kk != i and M[kk][i] != 0:
                f = M[kk][i]
                M[kk] = [a - f * b for a, b in zip(M[kk], M[i])]
    return [row[n:] for row in M]


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        while seq[i] != i:
            j = seq[i]
            seq[i], seq[j] = seq[j], seq[i]
            s = -s
    return s


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


# ----------------------------------------------------------------------------- network
class Network:
    def __init__(self, edges, verts):
        self.edges = edges
        self.verts = verts
        self.vi = {v: i for i, v in enumerate(verts)}
        self.B = [[(1 if e[0] == v else 0) - (1 if e[1] == v else 0) for v in verts] for e in edges]

    def solve(self, t, lam):
        """lam: dict v -> 3-vector. Returns y[e][a], phi[v][a], X = L^{-1} B^T."""
        n, m = len(self.verts), len(self.edges)
        L = [[sum(Fr(self.B[e][i] * self.B[e][j]) / t[e] for e in range(m)) for j in range(n)] for i in range(n)]
        rhs = [[Fr(lam.get(v, (0, 0, 0))[a]) for a in range(3)] for v in self.verts]
        BT = [[Fr(self.B[e][i]) for e in range(m)] for i in range(n)]
        sol = solve(L, [rhs[i] + BT[i] for i in range(n)])
        phi = [row[:3] for row in sol]
        X = [row[3:] for row in sol]
        y = [[sum(self.B[e][i] * phi[i][a] for i in range(n)) / t[e] for a in range(3)] for e in range(m)]
        return y, phi, X

    def dy_dt(self, t, lam):
        n, m = len(self.verts), len(self.edges)
        y, phi, X = self.solve(t, lam)
        # d phi / d t_j = (1/t_j^2) X[:, j] (b_j . phi)
        D = [[[Fr(0)] * m for _ in range(3)] for _ in range(m)]   # D[e][a][j]
        for j in range(m):
            bphi = [sum(self.B[j][i] * phi[i][a] for i in range(n)) for a in range(3)]
            dphi = [[X[i][j] * bphi[a] / t[j] ** 2 for a in range(3)] for i in range(n)]
            for e in range(m):
                for a in range(3):
                    val = sum(self.B[e][i] * dphi[i][a] for i in range(n)) / t[e]
                    if e == j:
                        val -= y[e][a] / t[e]
                    D[e][a][j] = val
        return y, D

    def tree_sign(self, tree):
        m = len(self.edges)
        tree = sorted(tree)
        cot = [e for e in range(m) if e not in tree]
        eps = perm_sign(tree + cot)
        BT = [[self.B[e][i] for e in tree] for i in range(len(self.verts))]
        return eps * det(BT), cot

    def is_tree(self, tree):
        par = {v: v for v in self.verts + ['O']}

        def f(x):
            while par[x] != x:
                x = par[x]
            return x
        for e in tree:
            a, b = self.edges[e]
            ra, rb = f(a), f(b)
            if ra == rb:
                return False
            par[ra] = rb
        return len(tree) == len(self.verts)

    def random_tree(self, rng):
        order = list(range(len(self.edges)))
        rng.shuffle(order)
        tree = []
        for e in order:
            if self._acyclic(tree + [e]):
                tree.append(e)
        assert self.is_tree(tree)
        return sorted(tree)

    def _acyclic(self, es):
        par = {v: v for v in self.verts + ['O']}

        def f(x):
            while par[x] != x:
                x = par[x]
            return x
        for e in es:
            a, b = self.edges[e]
            ra, rb = f(a), f(b)
            if ra == rb:
                return False
            par[ra] = rb
        return True


# ----------------------------------------------------------------------------- dual numbers
class Dual:
    __slots__ = ('v', 'g')

    def __init__(self, v, g):
        self.v = v
        self.g = g

    def __add__(self, o):
        o = lift(o, len(self.g))
        return Dual(self.v + o.v, [a + b for a, b in zip(self.g, o.g)])
    __radd__ = __add__

    def __sub__(self, o):
        o = lift(o, len(self.g))
        return Dual(self.v - o.v, [a - b for a, b in zip(self.g, o.g)])

    def __rsub__(self, o):
        return lift(o, len(self.g)) - self

    def __mul__(self, o):
        o = lift(o, len(self.g))
        return Dual(self.v * o.v, [self.v * b + o.v * a for a, b in zip(self.g, o.g)])
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = lift(o, len(self.g))
        return Dual(self.v / o.v, [(a * o.v - self.v * b) / o.v ** 2 for a, b in zip(self.g, o.g)])


def lift(x, k):
    return x if isinstance(x, Dual) else Dual(Fr(x), [Fr(0)] * k)


# ----------------------------------------------------------------------------- moves and chart
def find_edge(cur, a, b):
    ks = [k for k, (u, w, _) in cur.items() if {u, w} == {a, b}]
    assert len(ks) == 1, (a, b, ks)
    return ks[0]


def build_moves(name):
    """Forward pass: returns list of move records (kind, roles -> edge keys, merged key) and final box."""
    g = GRAPHS[name]
    cur = {('e', i): (e[0], e[1], None) for i, e in enumerate(g['edges'])}
    recs = []
    for mi, mv in enumerate(g['moves']):
        kind = mv[0]
        if kind == 'chain':
            p, v, v2, q = mv[1:]
            roles = dict(t1=find_edge(cur, p, v), t2=find_edge(cur, v, v2), t3=find_edge(cur, v2, q),
                         t0=find_edge(cur, p, q))
        elif kind == 'bubble':
            p, v, v2, q = mv[1:]
            roles = dict(pv=find_edge(cur, p, v), vq=find_edge(cur, v, q), pv2=find_edge(cur, p, v2),
                         v2q=find_edge(cur, v2, q))
        else:
            o, a, b, c = mv[1:]
            roles = dict(t1=find_edge(cur, o, a), t2=find_edge(cur, a, b), t3=find_edge(cur, b, c),
                         t4=find_edge(cur, c, o))
            recs.append((kind, mv, roles, None))
            continue
        for k in roles.values():
            del cur[k]
        mk = ('m', mi)
        cur[mk] = (p, q, None)
        recs.append((kind, mv, roles, mk))
    return recs


def chart(name, r):
    """t(r) as Duals (gradient w.r.t. r_1..r_k), r list of Fractions (k = 3 * #moves)."""
    k = len(r)
    R = [Dual(Fr(r[i]), [Fr(int(j == i)) for j in range(k)]) for i in range(k)]
    recs = build_moves(name)
    tval = {}
    # final box
    kind, mv, roles, _ = recs[-1]
    b = 3 * (len(recs) - 1)
    r1, r2, r3 = R[b], R[b + 1], R[b + 2]
    tval[roles['t1']] = 1 - r1
    tval[roles['t2']] = r1 - r2
    tval[roles['t3']] = r2 - r3
    tval[roles['t4']] = r3
    for mi in range(len(recs) - 2, -1, -1):
        kind, mv, roles, mk = recs[mi]
        r1, r2, r3 = R[3 * mi], R[3 * mi + 1], R[3 * mi + 2]
        tr = tval.pop(mk)
        if kind == 'chain':
            t0 = tr / r1
            T = tr / (1 - r1)
            tval[roles['t0']] = t0
            tval[roles['t1']] = T * (1 - r2)
            tval[roles['t2']] = T * (r2 - r3)
            tval[roles['t3']] = T * r3
        else:
            T2 = tr / r1
            T1 = tr / (1 - r1)
            tval[roles['vq']] = r2 * T1
            tval[roles['pv']] = (1 - r2) * T1
            tval[roles['v2q']] = r3 * T2
            tval[roles['pv2']] = (1 - r3) * T2
    m = len(GRAPHS[name]['edges'])
    ts = [tval[('e', i)] for i in range(m)]
    return [x.v for x in ts], [x.g for x in ts]


def density(name, r, lam, tree):
    g = GRAPHS[name]
    net = Network(g['edges'], g['verts'])
    t, dt = chart(name, r)
    y, D = net.dy_dt(t, lam)
    sT, cot = net.tree_sign(tree)
    k = len(r)
    m = len(t)
    rows = []
    for c in cot:
        for a in range(3):
            rows.append([sum(D[c][a][j] * dt[j][kk] for j in range(m)) for kk in range(k)])
    jac = det(rows)
    orient = det([t] + [[dt[j][kk] for j in range(m)] for kk in range(k)])
    return sT * jac * (1 if orient > 0 else -1), y, t


def exponent(name, y, z):
    """-sum_e y_e . z_e with z_e = sign * z_j on spoke j."""
    g = GRAPHS[name]
    tot = Fr(0)
    for j, (e, s) in enumerate(g['spokes']):
        tot -= s * sum(y[e][a] * z[j][a] for a in range(3))
    return tot


def rand_frac(rng, den=97):
    return Fr(rng.randint(1, den - 1), den)


def random_point(name, rng):
    """Random rational point of D_Gamma (the domain depends on move kinds)."""
    r = []
    for mv in GRAPHS[name]['moves']:
        if mv[0] == 'chain':
            a = rand_frac(rng)
            b, c = sorted([rand_frac(rng, 101), rand_frac(rng, 103)], reverse=True)
            r += [a, b, c]
        elif mv[0] == 'bubble':
            r += [rand_frac(rng), rand_frac(rng, 101), rand_frac(rng, 103)]
        else:
            r += sorted([rand_frac(rng), rand_frac(rng, 101), rand_frac(rng, 103)], reverse=True)
    return r


def reduction_factors(name, r, lam):
    """Sequential reduction with own solver on the reduced networks.  Returns the ratios recomputed
    from t(r), the product of local factors (chain r1^2 det(y_red, l_v, l_v'), bubble r1(1-r1)
    det(y_red, l_v, l_v'), box det(l_a, l_b, l_c)), and the reduced currents on retained edges
    compared with the full currents."""
    g = GRAPHS[name]
    t, _ = chart(name, r)
    edges = {('e', i): (e[0], e[1], t[i]) for i, e in enumerate(g['edges'])}
    src = {v: [Fr(x) for x in lam.get(v, (0, 0, 0))] for v in g['verts']}
    full_net = Network(g['edges'], g['verts'])
    yfull, _, _ = full_net.solve(t, {v: tuple(src[v]) for v in src})
    verts = list(g['verts'])
    prod = Fr(1)
    rr = []
    ok_currents = True
    recs = build_moves(name)
    for mi, (kind, mv, roles, mk) in enumerate(recs):
        if kind == 'box':
            o, a, b, c = mv[1:]
            T = sum(edges[roles[x]][2] for x in ('t1', 't2', 't3', 't4'))
            r1 = (edges[roles['t2']][2] + edges[roles['t3']][2] + edges[roles['t4']][2]) / T
            r2 = (edges[roles['t3']][2] + edges[roles['t4']][2]) / T
            r3 = edges[roles['t4']][2] / T
            rr += [r1, r2, r3]
            prod *= det3(src[a], src[b], src[c])
            break
        p, v, v2, q = mv[1:]
        if kind == 'chain':
            t1, t2, t3, t0 = (edges[roles[x]][2] for x in ('t1', 't2', 't3', 't0'))
            T = t1 + t2 + t3
            r1, r2, r3 = T / (T + t0), (t2 + t3) / T, t3 / T
            tred = T * t0 / (T + t0)
            w = r1 ** 2
        else:
            tpv, tvq, tpv2, tv2q = (edges[roles[x]][2] for x in ('pv', 'vq', 'pv2', 'v2q'))
            T1, T2 = tpv + tvq, tpv2 + tv2q
            r1, r2, r3 = T1 / (T1 + T2), tvq / T1, tv2q / T2
            tred = T1 * T2 / (T1 + T2)
            w = r1 * (1 - r1)
        rr += [r1, r2, r3]
        lv, lv2 = src[v], src[v2]
        for k in roles.values():
            del edges[k]
        edges[mk] = (p, q, tred)
        if p != 'O':
            src[p] = [src[p][a] + r2 * lv[a] + r3 * lv2[a] for a in range(3)]
        if q != 'O':
            src[q] = [src[q][a] + (1 - r2) * lv[a] + (1 - r3) * lv2[a] for a in range(3)]
        verts = [x for x in verts if x not in (v, v2)]
        del src[v], src[v2]
        keys = list(edges)
        net = Network([edges[k][:2] for k in keys], verts)
        yred, _, _ = net.solve([edges[k][2] for k in keys], {x: tuple(src[x]) for x in verts})
        for i, kk in enumerate(keys):
            if kk[0] == 'e' and yred[i] != yfull[kk[1]]:
                ok_currents = False
        ymerged = yred[keys.index(mk)]
        prod *= w * det3(ymerged, lv, lv2)
    return rr, prod, ok_currents


# ----------------------------------------------------------------------------- explorer formulas (C3, C4)
def vadd(*terms):
    out = [Fr(0)] * 3
    for c, v in terms:
        for a in range(3):
            out[a] += c * v[a]
    return out


def J_formula(name, r, l):
    r1, r2, r3, r4, r5, r6, r7, r8, r9 = r
    mu = vadd((1 - r2, l[1]), (1 - r3, l[2]))
    if name == 'G1':
        return (r1 ** 2 * r4 ** 2
                * det3(l[1], l[2], vadd((r6 + (1 - r4) * (1 - r6) * r7, l[4]), ((1 - r4) * r8, l[6]), ((1 - r4) * r9, l[7])))
                * det3(mu, l[4], vadd((r8, l[6]), (r9, l[7])))
                * det3(vadd((1 - r5, mu), (1 - r6, l[4])), l[6], l[7]))
    def d2():
        return (det3(l[1], l[2], vadd((r8, l[4]), (r9 * r5, l[6]), (r9 * r6, l[7])))
                * det3(vadd((1 - r7, mu), (1 - r8, l[4])), l[6], l[7])
                * det3(mu, l[4], vadd((r5, l[6]), (r6, l[7]))))

    def d3():
        return (det3(l[1], l[2], vadd((r7 * (1 - r6), l[5]), (r8, l[6]), (r9, l[7])))
                * det3(mu, l[5], vadd((r8, l[6]), (r9, l[7])))
                * det3(vadd((1 - r5, mu), (1 - r6, l[5])), l[6], l[7]))
    if name == 'G2':
        return r1 ** 2 * r4 * (1 - r4) * d2()
    if name == 'G3':
        return r1 ** 2 * r4 * (1 - r4) ** 2 * d3()
    if name == 'G4':
        return r1 * (1 - r1) * r4 * (1 - r4) * d2()
    if name == 'G5':
        return r1 * (1 - r1) * r4 * (1 - r4) ** 2 * d3()
    raise KeyError(name)


I_SIGN = {'G1': 1, 'G2': 1, 'G3': 1, 'G4': -1, 'G5': -1}
S_SIGN = {'G1': 1, 'G2': -1, 'G3': 1, 'G4': 1, 'G5': -1}


def X_formula(name, r, z):
    r1, r2, r3, r4, r5, r6, r7, r8, r9 = r
    z1, z2, z3, z4 = z
    W1 = vadd((1 - r1, z1), (r1, z2))
    if name == 'G1':
        W2 = vadd((1 - r4, W1), (r4, z3))
        P5 = vadd((1, z4), (r7, W2), (-r7, z4))
        P3 = vadd((r5, W1), (1 - r5, P5))
        return {1: vadd((r2, z1), (1 - r2, P3)), 2: vadd((r3, z1), (1 - r3, P3)),
                4: vadd((r6, W1), (1 - r6, P5)), 6: vadd((1, z4), (r8, W2), (-r8, z4)),
                7: vadd((1, z4), (r9, W2), (-r9, z4))}
    if name in ('G2', 'G4'):
        W2 = vadd((1 - r4, z3), (r4, z4))
        P3 = vadd((r7, W1), (1 - r7, W2))
        P5 = vadd((r9, W1), (1 - r9, W2))
        X2 = vadd((r3, z1), (1 - r3, P3)) if name == 'G2' else vadd((r3, z2), (1 - r3, P3))
        return {1: vadd((r2, z1), (1 - r2, P3)), 2: X2, 4: vadd((r8, W1), (1 - r8, W2)),
                6: vadd((r5, P5), (1 - r5, z3)), 7: vadd((r6, P5), (1 - r6, z4))}
    W2 = vadd((1 - r4, W1), (r4, z3))
    P4 = vadd((1, z4), (r7, W2), (-r7, z4))
    P3 = vadd((r5, W1), (1 - r5, P4))
    X2 = vadd((r3, z1), (1 - r3, P3)) if name == 'G3' else vadd((r3, z2), (1 - r3, P3))
    return {1: vadd((r2, z1), (1 - r2, P3)), 2: X2, 5: vadd((r6, z3), (1 - r6, P4)),
            6: vadd((1, z4), (r8, W2), (-r8, z4)), 7: vadd((1, z4), (r9, W2), (-r9, z4))}
