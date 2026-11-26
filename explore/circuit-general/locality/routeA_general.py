"""Task 2, route A for any tight two-terminal S (non-series-parallel allowed), any H, T = 0:

  K_Gamma = kappa int_{t_S} sum_beta c_beta(t_S) K_Gamma'[(J - s(t_S))^beta](lambda'_p = lambda_p + s(t_S)),
  D_S(t_S; sigma, lambda_S0) = sum_beta c_beta(t_S) sigma^beta,  c_beta = P_beta / U_S^(n_S + H),
  s = sum_v h_v lambda_v,  h_v = Hn_v / V,  V = U_{S/pq},
  int_{RP^{N_S-1}} P / (U_S^A V^B) = Gamma(A+B)/(Gamma(A) Gamma(B)) I_{A+B}(S + e*; alpha_S, alpha_* = A).

P_beta, Hn_v: exact interpolation (homogeneous polynomials) from htk.brute on S standalone;
K_Gamma'[J^b](h): htk.kernel_exact with symbolic h_v (Gamma' series-parallel);
periods: cores.period on S^ = S + e* (series reduction to W_n / prism / K33)."""
import itertools
import os
import random
import sys
import time
from fractions import Fraction as Fr
from math import comb, factorial

import sympy as sp

import locality as L
import htk
from htk import brute
from htk.graph import Graph
from htk.linalg import solve

sys.path.insert(0, os.path.join(L.HERE, '../cores'))
import cores                                                   # noqa: E402


def monomials(nvar, deg):
    if nvar == 1:
        yield (deg,)
        return
    for i in range(deg + 1):
        for rest in monomials(nvar - 1, deg - i):
            yield (i,) + rest


def interp_homog(f, nvar, deg, rng):
    """exact homogeneous interpolation: f(t) (Fraction) = sum_m c_m t^m, deg fixed"""
    mons = list(monomials(nvar, deg))
    rows, rhs = [], []
    while len(rows) < len(mons):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 7)) for _ in range(nvar)]
        rows.append([prod_pow(t, m) for m in mons])
        rhs.append([f(t)])
    c = solve(rows, rhs)
    out = {m: v[0] for m, v in zip(mons, c) if v[0] != 0}
    for _ in range(2):                       # verification at independent points
        t = [Fr(rng.randint(1, 60), rng.randint(1, 11)) for _ in range(nvar)]
        assert sum((cc * prod_pow(t, m) for m, cc in out.items()), Fr(0)) == f(t), 'interpolation failed'
    return out


def prod_pow(t, m):
    out = Fr(1)
    for x, k in zip(t, m):
        if k:
            out *= x ** k
    return out


class PolyD(dict):
    """sparse polynomial {exponent tuple: Fraction}"""

    def mul(self, other):
        out = PolyD()
        for m1, c1 in self.items():
            for m2, c2 in other.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                v = out.get(m, 0) + c1 * c2
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
        return out

    def pw(self, k, cache):
        key = (id(self), k)
        if key not in cache:
            if k == 0:
                nv = len(next(iter(self)))
                cache[key] = PolyD({(0,) * nv: Fr(1)})
            else:
                cache[key] = self.pw(k - 1, cache).mul(self)
        return cache[key]


def to_sympy(P, T):
    return sp.Poly.from_dict({m: sp.Rational(c.numerator, c.denominator) for m, c in P.items()}, *T)


class Dressed:
    def __init__(self, tt, lam_full, H, seed=1):
        self.tt, self.H = tt, H
        GS = tt.GS
        self.GS = GS
        N, nS = GS.ne, GS.loops
        self.N, self.nS = N, nS
        rng = random.Random(seed)
        sp_ = tt.sp
        self.lam0 = [lam_full[v] for v in sp_.S0]
        self.src_int = [a for a, v in enumerate(sp_.S0) if any(x != 0 for x in lam_full[v])]
        # U_S, V = U_{S/pq} (= R_eff U_S)
        C = GS.cycle_matrix()

        def U_of(t):
            return brute.symanzik_U(GS, t, C)

        def V_of(t):
            R, h, Li = tt.reff(t)
            return R * U_of(t)
        self.U = PolyD(interp_homog(U_of, N, nS, rng))
        self.V = PolyD(interp_homog(V_of, N, nS + 1, rng))
        # harmonic measures: Hn_v = h_v V
        self.Hn = []
        for a in self.src_int:
            self.Hn.append(PolyD(interp_homog(lambda t, a=a: tt.reff(t)[1][a] * V_of(t), N, nS + 1, rng)))
        # D_S as polynomial in sigma: coefficients via a tensor grid in sigma
        degs = [nS] * H
        grid = list(itertools.product(*[range(d + 1) for d in degs]))
        self.betas = grid
        dP = (nS + H) * nS - N

        def coeffs_at(t):
            vals = []
            for sg in grid:
                lamS = [list(x) for x in self.lam0] + [list(sg)]
                vals.append(brute.density(GS, t, lamS).coeff)
            # tensor Vandermonde on nodes 0..nS
            V1 = [[Fr(x) ** k for k in range(nS + 1)] for x in range(nS + 1)]
            Vi = L.inverse(V1)
            out = {}
            for b in grid:
                acc = Fr(0)
                for gi, sg in enumerate(grid):
                    w = Fr(1)
                    for i in range(H):
                        w *= Vi[b[i]][sg[i]]
                    acc += w * vals[gi]
                out[b] = acc
            return out
        # interpolate P_beta = U^(nS+H) c_beta (degree dP) simultaneously
        mons = list(monomials(N, dP))
        rows, rhs = [], []
        while len(rows) < len(mons):
            t = [Fr(rng.randint(1, 40), rng.randint(1, 7)) for _ in range(N)]
            u = U_of(t) ** (nS + H)
            cf = coeffs_at(t)
            rows.append([prod_pow(t, m) for m in mons])
            rhs.append([cf[b] * u for b in grid])
        sol = solve(rows, rhs)
        for _ in range(2):                   # verification at independent points
            t = [Fr(rng.randint(1, 60), rng.randint(1, 11)) for _ in range(N)]
            u = U_of(t) ** (nS + H)
            cf = coeffs_at(t)
            for j, b in enumerate(grid):
                assert sum((sol[i][j] * prod_pow(t, m) for i, m in enumerate(mons)), Fr(0)) == cf[b] * u, \
                    'P_beta interpolation failed'
        self.P = {}
        for j, b in enumerate(grid):
            P = PolyD({m: sol[i][j] for i, m in enumerate(mons) if sol[i][j] != 0})
            if P:
                self.P[b] = P
        self.A0 = nS + H

    def moment(self, P, A, B, Sedges, pq):
        """int P / (U^A V^B) over RP^{N-1} via periods of S^ = S + e*"""
        T = sp.symbols('t0:%d' % self.N)
        Pp = to_sympy(P, T)
        Up, Vp = to_sympy(self.U, T), to_sympy(self.V, T)
        while A > 1:
            q_, r_ = sp.div(Pp, Up)
            if not r_.is_zero:
                break
            Pp, A = q_, A - 1
        while B > 1:
            q_, r_ = sp.div(Pp, Vp)
            if not r_.is_zero:
                break
            Pp, B = q_, B - 1
        if B == 0:
            Pp, B = Pp * Vp, 1
        if A == 0:
            Pp, A = Pp * Up, 1
        Q = A + B
        pref = Fr(factorial(Q - 1), factorial(A - 1) * factorial(B - 1))
        edges = list(Sedges) + [pq]
        out = {}
        for mon, c in Pp.terms():
            alpha = [k + 1 for k in mon] + [A]
            key = (tuple(map(tuple, edges)), tuple(alpha), Q)
            if key not in _pc:
                _pc[key] = cores.period(edges, alpha, Q)
            cf = Fr(int(sp.numer(c)), int(sp.denom(c))) * pref
            for lb, v in _pc[key].items():
                out[lb] = out.get(lb, Fr(0)) + cf * v
        return out, len(Pp.terms())


_pc = {}


def run(G, S, lam, name, ref=None, H=None):
    t0 = time.time()
    tt = L.TwoTerminal(G, S)
    lf, H = L.lam_full(G, lam)
    rng = random.Random(3)
    kappa = None
    for _ in range(6):
        tS = L.rand_times(len(S), rng)
        tp = L.rand_times(len(tt.sp.O) + 1, rng)
        DG, JPhi, DSv, Dp, sigma = tt.check_point(tS, tp, lam)
        if DSv.coeff * Dp.coeff != 0:
            kappa = DG.coeff * abs(JPhi) / (DSv.coeff * Dp.coeff)
            assert kappa in (1, -1), kappa
            break
        assert DG.coeff == 0
    if kappa is None:                     # density vanishes at 6 points: kappa irrelevant if K = 0
        kappa = 1
    dr = Dressed(tt, lf, H)
    t1 = time.time()
    # check the interpolated D_S at the chart point
    U0 = sum((c * prod_pow(tS, m) for m, c in dr.U.items()), Fr(0))
    val = sum((sum((c * prod_pow(tS, m) for m, c in P.items()), Fr(0)) * prod_pow(sigma, b)
               for b, P in dr.P.items()), Fr(0)) / U0 ** dr.A0
    assert val == DSv.coeff, (val, DSv.coeff)     # interpolated D_S = brute D_S at the last chart point
    # Gamma' kernels: symbolic h_v
    hs = sp.symbols('h0:%d' % len(dr.src_int))
    s = [sum((hs[j] * dr.lam0[a][i] for j, a in enumerate(dr.src_int)), 0) for i in range(H)]
    tot0 = [sum((x[i] for x in dr.lam0), Fr(0)) for i in range(H)]
    lamp = []
    for v in tt.vp:
        x = [sp.Rational(str(c)) for c in lf[v]]
        if v == tt.p:
            x = [x[i] + s[i] for i in range(H)]
        if v == tt.q:
            x = [x[i] + sp.Rational(str(tot0[i])) - s[i] for i in range(H)]
        lamp.append(x)
    estar = tt.Gp.ne - 1
    Kcache = {}

    def Kp(bp):
        if bp not in Kcache:
            powers = {(estar, i): k for i, k in enumerate(bp) if k}
            Kcache[bp] = sp.expand(htk.kernel_exact(tt.Gp, lamp, T=0, powers=powers or None, check=0))
        return Kcache[bp]
    # sum over beta: K'[(J - s)^beta] = sum_{b' <= b} prod binom(b_i, b'_i) (-s_i)^(b_i - b'_i) K'[J^b']
    total = sp.Integer(0)
    poly_by_beta = {}
    for b in dr.P:
        expr = sp.Integer(0)
        for bp in itertools.product(*[range(k + 1) for k in b]):
            cf = sp.Integer(1)
            for i in range(H):
                cf *= comb(b[i], bp[i]) * (-s[i]) ** (b[i] - bp[i])
            expr += cf * Kp(bp)
        poly_by_beta[b] = sp.Poly(sp.expand(expr), *hs) if hs else sp.Poly(expr, sp.Symbol('dummy'))
    t2 = time.time()
    res = {}
    nmon = 0
    pcache = {}
    Sedges = tt.GS.edges
    pq = (tt.Sidx[tt.p], tt.Sidx[tt.q])
    for b, PB in poly_by_beta.items():
        for mon, cf in PB.terms():
            if not hs:
                mon = ()
            P = dr.P[b]
            for j, k in enumerate(mon):
                if k:
                    P = P.mul(dr.Hn[j].pw(k, pcache))
            mom, nm = dr.moment(P, dr.A0, sum(mon), Sedges, pq)
            nmon += nm
            c = Fr(int(sp.numer(cf)), int(sp.denom(cf)))
            for lb, v in mom.items():
                res[lb] = res.get(lb, Fr(0)) + kappa * c * v
    res = {k: v for k, v in res.items() if v}
    t3 = time.time()
    print('%s: locality K = %s | reference %s | equal %s | kappa=%d | P_beta terms %d, moment monomials %d, '
          'S^ periods %d | S data %.1fs, Gamma\' %.1fs, moments %.1fs' % (
              name, {k: str(v) for k, v in res.items()}, {k: str(v) for k, v in (ref or {}).items()},
              (ref == res) if ref is not None else '-', kappa, sum(len(P) for P in dr.P.values()), nmon,
              len(_pc), t1 - t0, t2 - t1, t3 - t2), flush=True)
    return res


if __name__ == '__main__' and (len(sys.argv) < 2 or sys.argv[1] == 'n3'):
    import fastker
    from check_locality import k4sub, relabel_root
    which = sys.argv[1] if len(sys.argv) > 1 else 'n3'
    if which == 'n3':
        rng = random.Random(41)
        for chains in ([3, 2, 1, 1, 2, 1], [3, 1, 2, 1, 1, 2]):
            nv, ed, ch = k4sub(chains, 3)
            G, new = relabel_root(nv, ed, 4)
            S = [e for c in ch[1:] for e in c]
            for _ in range(2):
                lam = [[rng.randint(-3, 3) for _ in range(3)] for _ in range(G.nv - 1)]
                full = lam + [[-sum(x[i] for x in lam) for i in range(3)]]
                t = time.time()
                ref = fastker.kernel(G.nv, G.edges, full, verbose=False)
                print('  reference fastker (whole graph, N=%d) %.1fs' % (G.ne, time.time() - t), flush=True)
                run(G, S, lam, 'D3 n3 K4%s bridge->box' % chains, ref={k: Fr(v) for k, v in ref.items()})


def run_logged(setname, idx):
    """bridge graphs of the D = 3, n = 4 sets with the sources of atlas/e3_open_kernels.py; if the root
    is internal to S, relabel (K -> sgn(pi)^D K, CONVENTIONS 2.2; checked with brute at one point)."""
    import json
    sys.path.insert(0, os.path.join(L.HERE, '../atlas'))
    from atlas import set_bulk3, set_twistor
    from htk.graph import from_networkx
    from routeB_bubble import e3_sources
    from htk.linalg import perm_sign
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
    pc = [x for x in rec['tightS_nonSP'] if x['NS'] == 7][0]
    C = [lab[v] for v in pc['C']]
    sign = 1
    G, lam_full = Gh, full
    if Gh.root in C:
        newroot = [v for v in range(Gh.nv) if v not in C and v not in (lab[pc['p']], lab[pc['q']])][0]
        perm = [v for v in range(Gh.nv) if v != newroot] + [newroot]       # new label i <- old perm[i]
        new = {v: i for i, v in enumerate(perm)}
        G = Graph(Gh.nv, [(new[a], new[b]) for a, b in Gh.edges])
        lam_full = [full[perm[i]] for i in range(Gh.nv)]
        C = [new[v] for v in C]
        # K(G) = sgn^D K(Gh): pointwise density check
        rng = random.Random(1)
        t = L.rand_times(Gh.ne, rng)
        d0 = brute.density(Gh, t, full).coeff
        d1 = brute.density(G, t, lam_full).coeff
        rule = perm_sign([perm[i] for i in range(Gh.nv)]) ** 3        # CONVENTIONS 2.2, D = 3
        if d0 != 0:
            sign = d1 / d0
            assert sign == rule
        else:
            assert d1 == 0
            sign = rule
    S = [e for e, (a, b) in enumerate(G.edges) if a in C or b in C]
    logs = {}
    for fn in ('e3_twistor4_s0.log', 'e3_twistor4_s1.log', 'e3_bulk3_s0.log', 'e3_bulk3_s1.log'):
        for line in open(os.path.join(L.HERE, '../atlas', fn)):
            if 'fastker K' in line:
                logs[line.split(' SP=')[0]] = eval(line.split('fastker K = ')[1].split(' (')[0])
    key = ('bulk3(3,4) #%d' % idx) if setname == 'bulk3' else ('twistor4 #%d' % idx)
    ref = {k: Fr(v) * sign for k, v in logs[key].items()}
    print('%s: root relabel sign %d' % (key, sign), flush=True)
    return run(G, S, lam_full[:-1], key + ' bridge (N_S=7) -> Gamma\' (2 loops, SP)', ref=ref)


if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[1] == 'logged':
    for arg in sys.argv[2:]:
        st, ix = ('bulk3', int(arg[1:])) if arg.startswith('b') else ('twistor4', int(arg))
        run_logged(st, ix)
