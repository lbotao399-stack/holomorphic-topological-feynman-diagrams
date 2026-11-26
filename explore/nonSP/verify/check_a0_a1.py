"""Exact pointwise checks of claims a0 (chain + sourced star-triangle currents) and a1 (density formula).
The formulas are re-implemented here from the claim text; the reference is vk.currents_exact / vk.density_exact
(fresh Kirchhoff solver on the full graph) and the sympy Jacobian of the explicit map xi -> t.
Coordinates xi = (s_{E,k} (chains in order, k = 1..a-1), gA, gB, rA, rB, rC); gC = 1 - gA - gB;
star edge Dc-v: R = 1/g_v; triangle edge w-x: R = (1 - r_v)/(r_v g_w g_x), v the third triangle vertex;
t_{E,j} = R_E (s_{j-1} - s_j), s_0 = 1, s_a = 0."""
import sys, random
from fractions import Fraction as Fr
import sympy as sp
from vk import G, currents_exact, density_exact, fdet
from hepp import chains_of
sys.path.insert(0, '..')
from graphs import K4SUB

C2N3 = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 4)]


def build(nv, E, lam, Dc):
    """lam: dict v -> list (H); returns symbols, t(xi), y(xi) (formula), Q, N, chains."""
    H = len(next(iter(lam.values())))
    cubic, chains = chains_of(nv, E)
    tri = [v for v in cubic if v != Dc]
    A, B, C = tri
    gA, gB, rA, rB, rC = sp.symbols('gA gB rA rB rC')
    g = {A: gA, B: gB, C: 1 - gA - gB}
    r = {A: rA, B: rB, C: rC}
    sv = {}
    for ci, (path, es) in enumerate(chains):
        for k in range(1, len(es)):
            sv[(ci, k)] = sp.Symbol('s%d_%d' % (ci, k))
    xi = [sv[k] for k in sorted(sv)] + [gA, gB, rA, rB, rC]
    third = lambda X, Y: [u for u in tri if u not in (X, Y)][0]
    def Rof(p, q):
        if Dc in (p, q):
            v = q if p == Dc else p
            return 1 / g[v]
        v = third(p, q)
        return (1 - r[v]) / (r[v] * g[p] * g[q])
    t = [None] * len(E)
    for ci, (path, es) in enumerate(chains):
        a = len(es)
        s = [sp.Integer(1)] + [sv[(ci, k)] for k in range(1, a)] + [sp.Integer(0)]
        R = Rof(path[0], path[-1])
        for j, e in enumerate(es, start=1):
            t[e] = sp.together(R * (s[j - 1] - s[j]))
    L = {v: sp.Matrix([sp.Rational(x) for x in lam.get(v, [0] * H)]) for v in range(nv)}
    mu = {v: L[v] for v in cubic}
    mhat = {}; Qs = {}
    for ci, (path, es) in enumerate(chains):
        p, q = path[0], path[-1]; a = len(es)
        mh = sp.zeros(H, 1)
        for k in range(1, a):
            s = sv[(ci, k)]
            mu[p] = mu[p] + s * L[path[k]]; mu[q] = mu[q] + (1 - s) * L[path[k]]
            mh = mh + s * L[path[k]]
        mhat[ci] = mh
        Qs[ci] = [sum((L[path[k]] for k in range(1, j)), sp.zeros(H, 1)) for j in range(1, a + 1)]
    p_ = {v: (1 - r[v]) * g[v] for v in tri}
    Q = p_[A] + p_[B] + p_[C]
    mup = {v: mu[v] + g[v] * mu[Dc] for v in tri}
    QY = lambda X, Y: p_[Y] * mup[X] - p_[X] * mup[Y]          # Q * Y_XY
    def QYK4(p, q):                                              # Q * (K4 current p -> q)
        if p == Dc:
            v = q
            return Q * g[v] * mu[Dc] - sum(((1 - r[third(v, w)]) * QY(v, w) for w in tri if w != v), sp.zeros(H, 1))
        if q == Dc:
            return -QYK4(q, p)
        return r[third(p, q)] * QY(p, q)
    Qy = [None] * len(E)
    N = []
    for ci, (path, es) in enumerate(chains):
        QYE = QYK4(path[0], path[-1])
        for j, e in enumerate(es, start=1):
            val = QYE - Q * mhat[ci] + Q * Qs[ci][j - 1]
            sgn = 1 if E[e] == (path[j - 1], path[j]) else -1
            Qy[e] = (sgn * val).expand()
        if path[0] in tri and path[-1] in tri:
            N.append((QYE - Q * mhat[ci]).expand())             # first edge, chain direction
    Nflat = [N[c][i] for i in range(H) for c in range(3)]
    return xi, t, Qy, sp.expand(Q), Nflat, chains, tri


def rand_point(rng, xi):
    vals = {}
    nchainvars = len(xi) - 5
    # s_{E,k}: decreasing in k inside each chain
    names = [str(x) for x in xi[:nchainvars]]
    groups = {}
    for x in xi[:nchainvars]:
        ci = str(x).split('_')[0]
        groups.setdefault(ci, []).append(x)
    for ci, xs in groups.items():
        ss = sorted([Fr(k, 98) for k in rng.sample(range(1, 98), len(xs))], reverse=True)
        for x, v in zip(xs, ss):
            vals[x] = v
    while True:
        a, b = Fr(rng.randint(1, 40), 61), Fr(rng.randint(1, 40), 61)
        if a + b < 1:
            break
    vals[xi[-5]] = a; vals[xi[-4]] = b
    for x in xi[-3:]:
        vals[x] = Fr(rng.randint(1, 30), 31)
    return vals


def tofr(v):
    v = sp.nsimplify(v) if not isinstance(v, (sp.Rational, sp.Integer)) else v
    return Fr(int(sp.numer(v)), int(sp.denom(v)))


def check_case(nv, E, lam, Dc, rng, npts=3):
    xi, t, Qy, Q, N, chains, tri = build(nv, E, lam, Dc)
    g = G(nv, E)
    H = len(next(iter(lam.values())))
    lamlist = [[Fr(x) for x in lam.get(v, [0] * H)] for v in range(nv - 1)]
    Jt = sp.Matrix([t] + [[sp.diff(ti, x) for ti in t] for x in xi])
    dQ = [sp.diff(Q, x) for x in xi]
    dN = [[sp.diff(nn, x) for x in xi] for nn in N]
    cur_ok = True; ratios = []
    for _ in range(npts):
        pt = rand_point(rng, xi)
        sub = {k: sp.Rational(v.numerator, v.denominator) for k, v in pt.items()}
        tv = [tofr(ti.subs(sub)) for ti in t]
        assert all(x > 0 for x in tv)
        yb = currents_exact(g, tv, lamlist)
        Qv = tofr(Q.subs(sub))
        for e in range(len(E)):
            for i in range(H):
                if tofr(Qy[e][i].subs(sub)) / Qv != yb[e][i]:
                    cur_ok = False
        if H == 3:
            Jv = [[tofr(x.subs(sub)) for x in row] for row in Jt.tolist()]
            gref = density_exact(g, tv, lamlist) * abs(fdet(Jv))
            Mx = [[Qv] + [tofr(x.subs(sub)) for x in dQ]] + \
                 [[tofr(N[k].subs(sub))] + [tofr(x.subs(sub)) for x in dN[k]] for k in range(9)]
            gf = fdet(Mx) / Qv ** 10
            ratios.append(gref / gf if gf != 0 else None)
    return cur_ok, ratios


if __name__ == "__main__":
    rng = random.Random(2026)
    print("C^3 K4 subdivisions, sources e1,e2,e3,-(sum) at the degree-two vertices:")
    for name, E in K4SUB:
        deg = [0] * 8
        for a, b in E:
            deg[a] += 1; deg[b] += 1
        d2 = [v for v in range(8) if deg[v] == 2]
        lam = {d2[0]: [1, 0, 0], d2[1]: [0, 1, 0], d2[2]: [0, 0, 1], d2[3]: [-1, -1, -1]}
        cubic = [v for v in range(8) if deg[v] == 3]
        for Dc in cubic:
            ok, rat = check_case(8, E, lam, Dc, rng)
            print("  %s Dc=%d  currents(formula) == Kirchhoff: %s   g_ref / [det/Q^10]: %s" % (name, Dc, ok, [str(x) for x in rat]), flush=True)
    print("C^3 K4 subdivisions, generic sources at all vertices (currents only):")
    for name, E in K4SUB[:2]:
        lam = {v: [rng.randint(-4, 4) for _ in range(3)] for v in range(7)}
        lam[7] = [-sum(lam[v][i] for v in range(7)) for i in range(3)]
        for Dc in [v for v in range(8) if sum(1 for e in E if v in e) == 3][:2]:
            ok, _ = check_case(8, E, lam, Dc, rng, npts=2)
            print("  %s Dc=%d currents ok: %s" % (name, Dc, ok), flush=True)
    print("C^2 n=3 graph, sources lam3, lam0, lam2 generic, lam1 = 0, lam4 = -(sum):")
    lam = {3: [1, 0], 0: [0, 1], 2: [Fr(1, 3), Fr(2, 5)], 1: [0, 0]}
    lam[4] = [-(1 + Fr(1, 3)), -(1 + Fr(2, 5))]
    for Dc in (0, 1, 2, 4):
        ok, _ = check_case(5, C2N3, lam, Dc, rng)
        print("  Dc=%d currents ok: %s" % (Dc, ok), flush=True)
