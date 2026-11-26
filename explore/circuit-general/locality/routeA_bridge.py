"""Task 2, dressed non-series-parallel two-terminal element: the Wheatstone bridge S = K4 - pq
(terminals p, q; internal a, b; tight at H = 2).  Dressed-edge data = moments of S computed with the
exact core evaluator (closure S^ = S + e* = K4) through

    int_{RP^4} P(t) / (U_S^A U_{S/pq}^B) = Gamma(A+B)/(Gamma(A) Gamma(B)) I_{A+B}(K4; alpha_S, alpha_* = A),
    U_{S^} = t_* U_S + U_{S/pq},  R_eff = U_{S/pq} / U_S.

K_Gamma(g lambda) = kappa sum_beta sum_{i,j} M[beta; i, j] [alpha^i beta^j] K_Gamma'[(y^1_*)^b1 (y^2_*)^b2](lambda'(alpha, beta)),
D_S = sum_beta c_beta(t) (J^1)^b1 (J^2)^b2 (J = y_{e*}), M[beta; i, j] = int c_beta h_a^i h_b^j,
g lambda_a = e1, g lambda_b = e2, s = h_a e1 + h_b e2, lambda'_p = lambda_p + s, lambda'_q = lambda_q + e1 + e2 - s.
Gamma' = O + e* (series-parallel): htk.kernel_exact with symbolic (alpha, beta).
usage: python3 routeA_bridge.py"""
import json
import os
import random
import sys
import time
from fractions import Fraction as Fr
from math import factorial

import sympy as sp

import locality as L
from htk import brute
import htk
from htk.graph import Graph
from htk.linalg import det, inverse

sys.path.insert(0, os.path.join(L.HERE, '../cores'))
import cores                                             # noqa: E402

K4E = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]   # Shat vertices: p=0, q=1, a=2, b=3


class Bridge:
    """symbolic data of the bare bridge S (standalone GS from TwoTerminal: vertex order a, b, p, q)"""

    def __init__(self, tt, H):
        GS = tt.GS
        assert GS.nv == 4 and GS.ne == 5
        self.GS, self.H = GS, H
        T = sp.symbols('t0:5', positive=True)
        self.T = T
        g = 3
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
        C = GS.cycle_matrix()
        self.C = C
        M = sp.Matrix([[sum(C[r][e] * C[s][e] * T[e] for e in range(5)) for s in range(2)] for r in range(2)])
        self.U = sp.expand(sp.factor(M.det()))
        R = sp.factor(Li[2, 2])                                  # q grounded, p index 2
        self.V = sp.expand(sp.cancel(R * self.U))
        # unit flow p -> q
        phi = [Li[x, 2] for x in range(3)] + [0]
        self.a = [sp.cancel((phi[GS.edges[e][0]] - phi[GS.edges[e][1]]) / T[e]) for e in range(5)]
        # Dirichlet: p, q grounded, unit sources at a (0), b (1)
        L00 = Lm[:2, :2]
        L00i = L00.inv()
        self.g = []
        for v in range(2):
            ph = [L00i[x, v] for x in range(2)] + [0, 0]
            self.g.append([sp.cancel((ph[GS.edges[e][0]] - ph[GS.edges[e][1]]) / T[e]) for e in range(5)])
        # harmonic measure of p at a, b: E = -L00^{-1} L0p
        L0p = Lm[:2, 2]
        E = -(L00i * L0p)
        self.h = [sp.cancel(E[0]), sp.cancel(E[1])]

    def numerators(self):
        """polynomial numerators: a_e = A_e / U, g_{v,e} = G_{v,e} / V, h_v = Hh_v / V"""
        if not hasattr(self, '_num'):
            A = [sp.expand(sp.cancel(x * self.U)) for x in self.a]
            Gn = [[sp.expand(sp.cancel(x * self.V)) for x in row] for row in self.g]
            Hh = [sp.expand(sp.cancel(x * self.V)) for x in self.h]
            self._num = (A, Gn, Hh)
        return self._num

    def DS(self, J):
        """D_S = det(Nmat) / (t_e0 U^H (U V)^(2H)) with polynomial Nmat; internal sources e1 at a,
        e2 at b (H = 2).  Returns (numerator polynomial in t and J, denominator)."""
        H, C, T = self.H, self.C, self.T
        A, Gn, Hh = self.numerators()
        U, V = self.U, self.V
        e0 = 4
        rows = [[sp.expand(C[r][e] * (A[e] * V * J[i] + (Gn[i][e] * U if i < 2 else 0)))
                 for e in range(5) if e != e0] for i in range(H) for r in range(2)]
        num = sp.expand(sp.Matrix(rows).det(method='berkowitz'))
        den = T[e0] * U ** H * (U * V) ** (2 * H)
        return (-num if e0 % 2 else num), den


def split_UV(expr, U, V, T):
    """expr = P / (U^A V^B), P polynomial: returns (Poly P, A, B)"""
    num, den = sp.fraction(sp.cancel(sp.together(expr)))
    den = sp.Poly(den, *T)
    Up, Vp = sp.Poly(U, *T), sp.Poly(V, *T)
    A = B = 0
    while den.total_degree() > 0:
        q2, r2 = sp.div(den, Vp)
        if r2.is_zero:
            den, B = q2, B + 1
            continue
        q1, r1 = sp.div(den, Up)
        if r1.is_zero:
            den, A = q1, A + 1
            continue
        raise ArithmeticError('denominator not of the form U^A V^B: %s' % den)
    return sp.Poly(sp.expand(num / den.as_expr()), *T), A, B


_pcache = {}


def period_K4(alpha, Q):
    key = (tuple(alpha), Q)
    if key not in _pcache:
        _pcache[key] = cores.period(K4E, list(alpha), Q)
    return _pcache[key]


def moment(br, num, A, B, divide_t=None):
    """int_{RP^4} num / (t_{divide_t} U^A V^B)  via K4 periods (num polynomial in T)"""
    T = br.T
    P = sp.Poly(num, *T)
    if divide_t is not None:
        q_, r_ = sp.div(P, sp.Poly(T[divide_t], *T))
        assert r_.is_zero
        P = q_
    Up, Vp = sp.Poly(br.U, *T), sp.Poly(br.V, *T)
    while A > 1:
        q_, r_ = sp.div(P, Up)
        if not r_.is_zero:
            break
        P, A = q_, A - 1
    while B > 1:
        q_, r_ = sp.div(P, Vp)
        if not r_.is_zero:
            break
        P, B = q_, B - 1
    if B == 0:
        P, B = P * Vp, 1
    if A == 0:
        P, A = P * Up, 1
    Q = A + B
    pref = Fr(factorial(Q - 1), factorial(A - 1) * factorial(B - 1))
    vmap = {0: 2, 1: 3, 2: 0, 3: 1}
    eidx = []
    for (u, w) in br.GS.edges:
        x, y = sorted((vmap[u], vmap[w]))
        eidx.append(K4E.index((x, y)))
    out = {}
    for mon, c in P.terms():
        alpha = [0] * 6
        for e, k in enumerate(mon):
            alpha[eidx[e]] = k + 1
        alpha[K4E.index((0, 1))] = A
        assert sum(alpha) == 3 * Q, (alpha, Q)
        val = period_K4(alpha, Q)
        cf = Fr(int(sp.numer(c)), int(sp.denom(c))) * pref
        for lb, v in val.items():
            out[lb] = out.get(lb, Fr(0)) + cf * v
    return out, (A, B, len(P.terms()))


def run(G, S, lam, name, ref=None):
    H = 2
    t0 = time.time()
    tt = L.TwoTerminal(G, S)
    a_v, b_v = tt.sp.S0
    lf, _ = L.lam_full(G, lam)
    la, lb = lf[a_v], lf[b_v]
    g = inverse([[la[0], lb[0]], [la[1], lb[1]]])
    lam_t = [[g[i][0] * x[0] + g[i][1] * x[1] for i in range(2)] for x in lf]
    n = G.loops
    rng = random.Random(3)
    tS = L.rand_times(5, rng)
    tp = L.rand_times(len(tt.sp.O) + 1, rng)
    DG, JPhi, DSv, Dp, sigma = tt.check_point(tS, tp, lam_t)
    kappa = DG.coeff * abs(JPhi) / (DSv.coeff * Dp.coeff)
    assert kappa in (1, -1)
    br = Bridge(tt, H)
    # check symbolic D_S against brute at the point
    Jv = [sigma[i] + sum((brh * lam_t[v][i] for brh, v in zip(
        [sp.Rational(x) for x in tt.reff(tS)[1]], tt.sp.S0)), 0) for i in range(H)]
    sub = {br.T[e]: sp.Rational(str(tS[e])) for e in range(5)}
    J1, J2 = sp.symbols('J1 J2')
    numJ, den = br.DS([J1, J2])
    chk = (numJ / den).subs(sub).subs({J1: sp.Rational(str(Jv[0])), J2: sp.Rational(str(Jv[1]))})
    assert sp.nsimplify(chk) == sp.Rational(str(DSv.coeff)), (chk, DSv)
    DSJ = sp.Poly(numJ, J1, J2)
    print('  D_S numerator: %d J-monomials, timing %.1fs' % (len(DSJ.terms()), time.time() - t0), flush=True)
    # Gamma' kernels with symbolic shift and insertion
    al, be = sp.symbols('alpha beta')
    estar = tt.Gp.ne - 1
    tot = {}
    nmom = 0
    for (b1, b2), cpoly in DSJ.terms():
        lamp = []
        for v in tt.vp:
            x = list(lam_t[v])
            if v == tt.p:
                x = [x[0] + al, x[1] + be]
            if v == tt.q:
                x = [x[0] + 1 - al, x[1] + 1 - be]
            lamp.append(x)
        powers = {}
        if b1:
            powers[(estar, 0)] = b1
        if b2:
            powers[(estar, 1)] = b2
        Kp = htk.kernel_exact(tt.Gp, lamp, T=0, powers=powers or None, check=1)
        Kp = sp.Poly(sp.expand(Kp), al, be)
        A_, Gn_, Hh = br.numerators()
        for (i, j), cf in Kp.terms():
            mom, info = moment(br, cpoly * Hh[0] ** i * Hh[1] ** j, 3 * H, 2 * H + i + j, divide_t=4)
            nmom += info[2]
            for lbl, v in mom.items():
                tot[lbl] = tot.get(lbl, Fr(0)) + Fr(int(sp.numer(cf)), int(sp.denom(cf))) * v
    detg = det(g)
    K = {k: kappa * v / detg ** n for k, v in tot.items() if v}
    t1 = time.time()
    print('%s: locality K = %s | reference %s | equal %s | kappa=%d | %d moment monomials, %d K4 periods, %.1fs'
          % (name, {k: str(v) for k, v in K.items()}, {k: str(v) for k, v in (ref or {}).items()},
             (ref == K) if ref is not None else '-', kappa, nmom, len(_pcache), t1 - t0), flush=True)
    return K


if __name__ == '__main__' and len(sys.argv) == 1:
    rng = random.Random(11)
    # D=2, n=3: K4 with p-q subdivided: p=0, q=1, a=2, b=3, v=4 (root)
    G = Graph(5, [(0, 2), (0, 3), (2, 1), (3, 1), (2, 3), (0, 4), (4, 1)])
    S = [0, 1, 2, 3, 4]
    for trial in range(3):
        lam = [[rng.randint(-4, 4), rng.randint(-4, 4)] for _ in range(4)]
        t = time.time()
        ref = cores.kernel(G.nv, G.edges, lam)
        print('  reference cores.kernel (whole graph) %.1fs' % (time.time() - t))
        run(G, S, lam, 'D2 n3 K4[2,1,1,1,1,1] sources %s' % lam, ref={k: Fr(v) for k, v in ref.items()})


def find_bridges(G):
    """two-terminal pieces S = K4 - pq (bare Wheatstone bridge) with internal a, b of degree 3 in G,
    root not internal, Gamma' = Gamma - {a, b} + e* series-parallel"""
    import networkx as nx
    out = []
    nb = {v: [] for v in range(G.nv)}
    for e, (x, y) in enumerate(G.edges):
        nb[x].append((y, e))
        nb[y].append((x, e))
    for a in range(G.nv - 1):
        for b in range(a + 1, G.nv - 1):
            if G.degree(a) != 3 or G.degree(b) != 3:
                continue
            na = sorted(x for x, _ in nb[a])
            nbb = sorted(x for x, _ in nb[b])
            if b not in na:
                continue
            pa = sorted(set(na) - {b})
            pb = sorted(set(nbb) - {a})
            if pa != pb or len(pa) != 2:
                continue
            S = sorted(set(e for _, e in nb[a]) | set(e for _, e in nb[b]))
            if len(S) != 5:
                continue
            Gp = nx.MultiGraph()
            Gp.add_edges_from([G.edges[e] for e in range(G.ne) if e not in S])
            Gp.add_edge(*pa)
            if htk.catalog.is_sp(Gp):
                out.append(S)
    return out


def main_n4(maxgraphs=6):
    sys.path.insert(0, os.path.join(L.HERE, '../cores'))
    import fastker
    rng = random.Random(29)
    done = 0
    for G in htk.catalog.tight_graphs(2, 4):
        brs = find_bridges(G)
        if not brs:
            continue
        S = brs[0]
        lam = [[rng.randint(-4, 4), rng.randint(-4, 4)] for _ in range(G.nv - 1)]
        lam.append([-sum(x[0] for x in lam), -sum(x[1] for x in lam)])
        t = time.time()
        ref = fastker.kernel(G.nv, G.edges, lam, verbose=False)
        tref = time.time() - t
        print('  D2 n4 graph %s: reference fastker (whole graph) %.1fs' % (G.edges, tref), flush=True)
        try:
            run(G, S, lam[:-1], 'D2 n4 bridge S=%s' % S, ref={k: Fr(v) for k, v in ref.items()})
        except Exception as ex:
            print('  FAILED', repr(ex), flush=True)
        done += 1
        if done >= maxgraphs:
            break


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'n4':
    main_n4(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
