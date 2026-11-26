"""Task 2, dressed tight series-parallel element (bubble) + 3-loop core Gamma'.

Gamma = S cup_{p,q} O, S = bubble (two chains p -> q, 4 edges, 1 loop, two sourced bivalent vertices
u, w; tight at H = 3).  Covariance K(g lambda) = det(g)^n K(lambda): g with g lambda_u = e1,
g lambda_w = e2.  Then D_S(t_S; sigma, e1, e2) = sigma^3 c(t_S) (rows of component 3 contain only
sigma), sigma = J - s, s = h_u e1 + h_w e2, s^3 = 0, hence

    K_Gamma(g lambda) = kappa sum_{a,b <= 3} m[a,b] [alpha^a beta^b] K_Gamma'[y^3_{e*}](lambda'(alpha, beta)),
    m[a,b] = int_{RP^3} c(t_S) h_u^a h_w^b,   lambda'_p = lambda_p + alpha e1 + beta e2,
    lambda'_q = lambda_q + (1 - alpha) e1 + (1 - beta) e2,

K_Gamma'[y^3_{e*}] of degree <= 3 in alpha and in beta (alpha enters only component-1 rows, 3 loops).
m[a,b]: exact (sympy, ratio coordinates of the bubble, polynomial integrand).  K_Gamma': cores.kernel
(exact, insertion y^3_{e*}) at the 16 points (alpha, beta) in {0..3}^2.  kappa: chart-level check.
usage: python3 routeB_bubble.py IDX [IDX ...]   (twistor4 indices; bulk3 as b7 etc.)"""
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

import locality as L
from htk import brute
from htk.graph import Graph, from_networkx
from htk.linalg import det, inverse

sys.path.insert(0, os.path.join(L.HERE, '../atlas'))
sys.path.insert(0, os.path.join(L.HERE, '../cores'))
H = 3


def e3_sources(Gh, lab, S0, root, seed):
    """identical to atlas/e3_open_kernels.py"""
    rng = random.Random(seed)
    lam = {lab[v]: [rng.randint(-3, 3) for _ in range(H)] for v in S0 if v != root}
    full = [lam.get(v, [0] * H) for v in range(Gh.nv - 1)]
    full.append([-sum(r[i] for r in full) for i in range(H)])
    return full


def bubble_moments(GS, u_idx, w_idx, amax=3):
    """m[a][b] = int c(t) h_u^a h_w^b over RP^3 for the bubble standalone graph GS (q = root),
    c = D_S at sources u: e1, w: e2, p: e3.  Coordinates: chains from p to q, t = R_i x_{i,j},
    R = (w0, 1 - w0); measure |det[t; dt/dxi]| (positive), polynomial integrand."""
    nv = GS.nv
    p, q = nv - 2, nv - 1
    # chains from p to q
    adj = {}
    for e, (a, b) in enumerate(GS.edges):
        adj.setdefault(a, []).append((b, e))
        adj.setdefault(b, []).append((a, e))
    chains = []
    for (v, e) in adj[p]:
        ch = [e]
        prev, cur = p, v
        verts = []
        while cur != q:
            verts.append(cur)
            nxt = [(x, f) for (x, f) in adj[cur] if f != ch[-1]][0]
            ch.append(nxt[1])
            prev, cur = cur, nxt[0]
        chains.append((ch, verts))
    assert len(chains) == 2 and sum(len(c) for c, _ in chains) == GS.ne
    w0 = sp.Symbol('w0')
    xs, tsub, hv = [], {}, {}
    Rs = [w0, 1 - w0]
    for ci, (ch, verts) in enumerate(chains):
        k = len(ch)
        xv = [sp.Symbol('x%d_%d' % (ci, j)) for j in range(k - 1)]
        xs.append(xv)
        fr = xv + [1 - sum(xv)]
        for j, e in enumerate(ch):
            tsub[e] = Rs[ci] * fr[j]
        for j, v in enumerate(verts):
            hv[v] = sum(fr[j + 1:])               # Dirichlet harmonic measure of p at v
    tsym = [tsub[e] for e in range(GS.ne)]
    xi = [w0] + [x for xv in xs for x in xv]
    Jm = sp.Matrix([tsym] + [[sp.diff(te, z) for te in tsym] for z in xi])
    jac = sp.factor(Jm.det())
    # D_S symbolic (brute formula with symbolic t)
    T = sp.symbols('T0:%d' % GS.ne, positive=True)
    g = nv - 1
    Lm = sp.zeros(g, g)
    for (a, b), te in zip(GS.edges, T):
        c = 1 / te
        if a < g:
            Lm[a, a] += c
        if b < g:
            Lm[b, b] += c
        if a < g and b < g:
            Lm[a, b] -= c
            Lm[b, a] -= c
    Li = Lm.inv()
    src = [[0] * 3 for _ in range(g)]
    src[u_idx] = [1, 0, 0]
    src[w_idx] = [0, 1, 0]
    src[p] = [0, 0, 1]
    phi = Li * sp.Matrix(src)
    y = []
    for e, (a, b) in enumerate(GS.edges):
        pa = phi.row(a) if a < g else sp.zeros(1, 3)
        pb = phi.row(b) if b < g else sp.zeros(1, 3)
        y.append([(pa[i] - pb[i]) / T[e] for i in range(3)])
    C = GS.cycle_matrix()
    N = GS.ne
    e0 = N - 1
    rows = [[C[0][e] * y[e][i] for e in range(N) if e != e0] for i in range(3)]
    U = sum(C[0][e] ** 2 * T[e] for e in range(N))
    D = sp.Matrix(rows).det() / (T[e0] * U ** 3)
    if e0 % 2:
        D = -D
    Dxi = sp.simplify(D.subs({T[e]: tsym[e] for e in range(N)}) * sp.Abs(jac) if False else
                      D.subs({T[e]: tsym[e] for e in range(N)}) * jac)
    Dxi = sp.factor(sp.cancel(Dxi))
    # sign of jac on the domain: evaluate at an interior point
    pt = {w0: sp.Rational(1, 3)}
    for xv in xs:
        for j, x in enumerate(xv):
            pt[x] = sp.Rational(1, len(xv) + 2)
    sgn = 1 if jac.subs(pt) > 0 else -1
    dens = sp.expand(sp.cancel(Dxi * sgn))
    assert dens.is_polynomial(*xi), dens
    # exact integration over [0,1] (w0) x simplices
    def integ(expr):
        expr = sp.expand(expr)
        res = sp.integrate(expr, (w0, 0, 1))
        for xv in xs:
            # simplex {x >= 0, sum <= 1}
            for j in range(len(xv) - 1, -1, -1):
                upper = 1 - sum(xv[:j])
                res = sp.integrate(res, (xv[j], 0, upper))
        return sp.nsimplify(sp.expand(res))
    huv = hv[u_idx], hv[w_idx]
    m = [[Fr(str(integ(dens * huv[0] ** a * huv[1] ** b))) for b in range(amax + 1)] for a in range(amax + 1)]
    return m, chains


def kernel_dict_add(acc, val, c):
    for k, v in val.items():
        acc[k] = acc.get(k, Fr(0)) + c * Fr(v)
    return acc


def run(setname, idx):
    import networkx as nx
    from atlas import set_bulk3, set_twistor
    import cores
    if setname == 'twistor4':
        tw = set_twistor(4)
        items = [it for k in ('tree', 'unicyclic') for it in tw[k]]
        G0, O, S0 = items[idx]
    else:
        G0 = set_bulk3(3, 4)[idx][0]
        S0 = [v for v in G0 if G0.degree(v) == 2]
        O = S0[-1]
    Gh, lab = from_networkx(G0, root=O)
    full = e3_sources(Gh, lab, S0, O, 1000 * idx + 7)
    pcs = json.load(open(os.path.join(L.HERE, 'pieces_d3n4.json')))
    rec = [s for s in pcs if s['set'] == setname and s['idx'] == idx][0]
    pc = rec['tightS_SP_coreO'][0]
    p, q, C = lab[pc['p']], lab[pc['q']], [lab[v] for v in pc['C']]
    Sed = [e for e, (a, b) in enumerate(Gh.edges) if a in C or b in C]
    t0 = time.time()
    # normalization g: lambda_u -> e1, lambda_w -> e2
    u, w = sorted(C)
    lu, lw = [Fr(x) for x in full[u]], [Fr(x) for x in full[w]]
    for v3 in ([0, 0, 1], [0, 1, 0], [1, 0, 0]):
        Mcol = [[lu[i], lw[i], Fr(v3[i])] for i in range(3)]
        if det(Mcol) != 0:
            break
    g = inverse(Mcol)
    lam_t = [[sum((g[i][j] * Fr(x[j]) for j in range(3)), Fr(0)) for i in range(3)] for x in full]
    n = Gh.loops
    # two-terminal split (root must not be internal to S: here the root is O or a bivalent vertex)
    tt = L.TwoTerminal(Gh, Sed, p=p, q=q)
    # kappa
    rng = random.Random(5)
    kappa = None
    for _ in range(6):
        tS = L.rand_times(len(Sed), rng)
        tp = L.rand_times(len(tt.sp.O) + 1, rng)
        DG, JPhi, DS, Dp, sigma = tt.check_point(tS, tp, lam_t)
        if DS.coeff * Dp.coeff != 0:
            kappa = DG.coeff * abs(JPhi) / (DS.coeff * Dp.coeff)    # integration sign incl. orientation of Phi
            assert kappa in (1, -1), kappa
            break
        assert DG.coeff == 0
    if kappa is None:
        kappa = 1
    # S standalone indices of u, w
    m, chains = bubble_moments(tt.GS, tt.Sidx[u], tt.Sidx[w])
    t1 = time.time()
    # Gamma' kernels with insertion y^3_{e*} at alpha, beta in 0..3
    Gp = tt.Gp
    estar = Gp.ne - 1
    vals = {}
    for al in range(4):
        for be in range(4):
            lamp = []
            for v in tt.vp:
                x = list(lam_t[v])
                if v == p:
                    x = [x[0] + al, x[1] + be, x[2]]
                if v == q:
                    x = [x[0] + 1 - al, x[1] + 1 - be, x[2]]
                lamp.append(x)
            vals[(al, be)] = cores.kernel(Gp.nv, Gp.edges, lamp, [(estar, 2)])
    t2 = time.time()
    # tensor interpolation in alpha, beta (nodes 0..3): coefficient extraction
    V = [[Fr(a) ** k for k in range(4)] for a in range(4)]
    Vi = inverse(V)
    labels = sorted({k for v in vals.values() for k in v})
    K = {}
    for lb in labels:
        F = [[Fr(vals[(a, b)].get(lb, 0)) for b in range(4)] for a in range(4)]
        coef = [[sum((Vi[i][a] * Vi[j][b] * F[a][b] for a in range(4) for b in range(4)), Fr(0))
                 for j in range(4)] for i in range(4)]
        K[lb] = kappa * sum((m[i][j] * coef[i][j] for i in range(4) for j in range(4)), Fr(0)) / det(g) ** n
    K = {k: v for k, v in K.items() if v}
    t3 = time.time()
    return dict(K=K, kappa=kappa, chains=[len(c) for c, _ in chains], NS=len(Sed), Gp_ne=Gp.ne,
                t_moments=t1 - t0, t_kernels=t2 - t1, t_total=t3 - t0)


if __name__ == '__main__':
    logs = {}
    for fn in ('e3_twistor4_s0.log', 'e3_twistor4_s1.log', 'e3_bulk3_s0.log', 'e3_bulk3_s1.log'):
        for line in open(os.path.join(L.HERE, '../atlas', fn)):
            if 'fastker K' in line:
                head = line.split(' SP=')[0]
                s = line.split('fastker K = ')[1].split(' (')[0]
                logs[head] = eval(s)
    for arg in sys.argv[1:]:
        setname, idx = ('bulk3', int(arg[1:])) if arg.startswith('b') else ('twistor4', int(arg))
        res = run(setname, idx)
        key = ('bulk3(3,4) #%d' % idx) if setname == 'bulk3' else ('twistor4 #%d' % idx)
        ref = {k: Fr(v) for k, v in logs[key].items()}
        ok = ref == res['K']
        print('%s: locality K = %s | log K = %s | equal: %s | kappa=%d chains=%s N_S=%d N\'=%d | '
              'moments %.1fs, 16 Gamma\' kernels %.1fs' % (
                  key, {k: str(v) for k, v in res['K'].items()}, {k: str(v) for k, v in ref.items()}, ok,
                  res['kappa'], res['chains'], res['NS'], res['Gp_ne'], res['t_moments'], res['t_kernels']),
              flush=True)
