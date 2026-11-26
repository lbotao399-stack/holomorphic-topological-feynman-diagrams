"""Task 1: exact pointwise checks of the core factorization against htk.brute.

T = 0 (any H):  for a tight polynomial partial reduction Gamma -> core (tight chain moves; D = 2 also
LLC (b), (c) source-free stars), coordinates xi of the moves and affine core times T,
    D_Gamma(t(xi,T)) J_aff (-1)^k = sigma_G^H eps(n,H) sigma_c^H eps(n_c,H) J_poly(xi, ell_core) D_core(T; lam_core(xi)),
J_aff = det d t / d(xi, T),  J_poly = det d ell_Gamma / d(xi, ell_core)  (polynomial in xi, ell_core, lambda).
LLC local factors: single (b): J_poly = +-r1 r2 det(y_ux, y_wx); paired (c): +-(r1-r2)(1-r3-r4)<y_uA,y_AB><y_uB,y_AB> g_m.

T = 1 (any H), k = 2 tight chain moves whose merged edges e*_m survive in the core:
    D_Gamma J_aff (-1)^k = s prod_m (w_m1 w_m2)^{-1/2} J_x(x, ell_core) D_core^{[F]},  F = prod_m (-sqrt(T_e*) w_e*),
D^{[F]} = U^{-H-1/2} <F prod_r psi^tau_r> ^ prod psi^i_r / iota_E vol (Gaussian w, covariance (delta/t - C^t M^-1 C)/2),
J_x = det d ell_Gamma / d(x, ell_core), x = series coordinates; s = +-1 constant (checked).
"""
import random
import sys
import time
from fractions import Fraction as Fr

import networkx as nx

from graphs_enum import admissible
from reduce_search import tight_search, Network, search
from polyreduce import PartialReduction, FullReduction, Dual, lam_dict
from htk.graph import from_networkx, Graph
from htk import brute, catalog
from htk.linalg import (det as fdet, inverse, one_form, wedge, wedge_all, form_add, form_scale,
                        basic_coefficient)
from htk.conventions import sign_dy


def nxg(G):
    H = nx.Graph()
    H.add_nodes_from(range(G.nv))
    H.add_edges_from(G.edges)
    return H


def rand_sources(G, S0, H, rng):
    lam = [[0] * H for _ in range(G.nv - 1)]
    for v in S0:
        if v != G.root:
            lam[v] = [rng.randint(-6, 6) for _ in range(H)]
    return lam


def rand_point(red, rng):
    xi = {}
    for blk in red.el.blocks:
        ns = [rng.randint(1, 30) for _ in blk]
        s = sum(ns)
        for b, kk in zip(blk, ns):
            xi[b] = Fr(kk, s)
    return xi


# ---------------------------------------------------------------------------- T = 0

def t0_checks(cases, H_of, npts=3, seed=1, log=print):
    rng = random.Random(seed)
    tot = ok = nonzero = 0
    llc_single = []
    llc_pair = []
    for D, n, idx, G, S0, plan, kind in cases:
        H = H_of(D)
        red = PartialReduction(G, plan)
        lam = rand_sources(G, S0, H, rng)
        good = True
        for _ in range(npts):
            xi = rand_point(red, rng)
            Tc = [Fr(rng.randint(1, 50), rng.randint(1, 50)) for _ in range(red.core.ne)]
            lhs, rhs, Jp = red.check(lam, xi, Tc)
            tot += 1
            if lhs != 0:
                nonzero += 1
            if lhs == rhs:
                ok += 1
            else:
                good = False
            if kind == "single" and H == 2:
                llc_single.append(_llc_single_ratio(red, lam, xi, Tc, Jp))
            if kind == "paired" and H == 2:
                llc_pair.append(_llc_pair_ratio(red, lam, xi, Tc, Jp))
        log("T=0 D=%d n=%d #%d %-7s plan %s core |V|=%d |E|=%d n_c=%d : %s" % (
            D, n, idx, kind, plan, red.core.nv, red.core.ne, red.core.loops, "exact" if good else "FAIL"))
    return tot, ok, [x for x in llc_single if x is not None], [x for x in llc_pair if x is not None], nonzero


def _core_y(red, lam, xi, Tc):
    G, el = red.G, red.el
    lam0, H = lam_dict(G, lam)
    bvn = el.bary_values({c: Fr(xi[c]) for c in el.coords}, Fr(1))
    _, talln = el.times(bvn, {e: Fr(Tc[i]) for i, e in enumerate(red.core_edges)})
    lamc = el.sources(lam0, bvn, talln, None)
    lam_core = [lamc[red.inv[i]] for i in range(red.core.nv - 1)]
    yc = brute.currents(red.core, [Fr(x) for x in Tc], lam_core)
    return {red.core_edges[e]: yc[e] for e in range(red.core.ne)}, bvn


def _dir(red, eid, a, b, y):
    """current of core edge eid in the direction a -> b"""
    p, q = red.el.final[eid]
    return y if (p, q) == (a, b) else [-x for x in y]


def _det2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _llc_single_ratio(red, lam, xi, Tc, Jp):
    st = red.el.steps[0]
    assert st.kind == "Y"
    yc, bv = _core_y(red, lam, xi, Tc)
    Ps = [ps for ps in red.el.steps[1:] if ps.kind == "P"]
    merged = {}
    for ps in Ps:
        orig = [e for e in ps.edges if e not in st.new][0]
        merged[frozenset((ps.p, ps.q))] = (ps, bv[ps.block[ps.edges.index(orig)]])
    nb = set(st.nbrs)
    pairs = [frozenset(p) for p in ((st.nbrs[0], st.nbrs[1]), (st.nbrs[0], st.nbrs[2]), (st.nbrs[1], st.nbrs[2]))]
    absent = [p for p in pairs if p not in merged]
    assert len(absent) == 1 and len(merged) == 2
    u, w = sorted(absent[0])
    x = (nb - {u, w}).pop()
    ps_ux, r1 = merged[frozenset((u, x))]
    ps_wx, r2 = merged[frozenset((w, x))]
    y_ux = _dir(red, ps_ux.new[0], u, x, yc[ps_ux.new[0]])
    y_wx = _dir(red, ps_wx.new[0], w, x, yc[ps_wx.new[0]])
    llc = r1 * r2 * _det2(y_ux, y_wx)
    if llc == 0:
        assert Jp == 0
        return None
    return Jp / llc


def _llc_pair_ratio(red, lam, xi, Tc, Jp):
    s1, s2 = red.el.steps[0], [s for s in red.el.steps if s.kind == "Y"][1]
    yc, bv = _core_y(red, lam, xi, Tc)
    N1, N2 = set(s1.nbrs), set(s2.nbrs)
    assert N1 == N2
    Ps = [s for s in red.el.steps if s.kind == "P"]
    # identify u, A, B: AB is merged after the first star
    first_P = [s for s in red.el.steps[1:red.el.steps.index(s2)] if s.kind == "P"]
    assert len(first_P) == 1
    A, B = first_P[0].p, first_P[0].q
    u = (N1 - {A, B}).pop()
    fP = first_P[0]
    orig_AB = [e for e in fP.edges if e < red.G.ne][0]
    f_orig = bv[fP.block[fP.edges.index(orig_AB)]]
    f_1 = 1 - f_orig
    later = {frozenset((s.p, s.q)): s for s in Ps if s is not fP}
    sAB = later[frozenset((A, B))]
    g_m = bv[sAB.block[sAB.edges.index(fP.new[0])]]
    g_2 = 1 - g_m
    def frac_v1(pair):
        s = later[frozenset(pair)]
        e1 = [e for e in s.edges if e in s1.new][0]
        return bv[s.block[s.edges.index(e1)]], s
    r1, suA = frac_v1((u, A))
    r2, suB = frac_v1((u, B))
    r3 = g_m * f_1
    r4 = g_2
    yuA = _dir(red, suA.new[0], u, A, yc[suA.new[0]])
    yuB = _dir(red, suB.new[0], u, B, yc[suB.new[0]])
    yAB = _dir(red, sAB.new[0], A, B, yc[sAB.new[0]])
    llc = (r1 - r2) * (1 - r3 - r4) * _det2(yuA, yAB) * _det2(yuB, yAB) * g_m
    if llc == 0:
        assert Jp == 0
        return None
    return Jp / llc


# ---------------------------------------------------------------------------- T = 1

def wick(items, cov):
    """E[prod items], items = (dict edge -> form, degree), linear in the centred Gaussian w."""
    if not items:
        return {0: Fr(1)}
    if len(items) % 2:
        return {}
    a, rest = items[0], items[1:]
    out = {}
    acc = 0
    for j in range(len(rest)):
        b = rest[j]
        sign = -1 if (b[1] * acc) % 2 else 1
        pf = {}
        for e, fa in a[0].items():
            for f, fb in b[0].items():
                c = cov[e][f]
                if c != 0:
                    pf = form_add(pf, form_scale(wedge(fa, fb), c))
        if pf:
            sub = wick(rest[:j] + rest[j + 1:], cov)
            if sub:
                out = form_add(out, form_scale(wedge(pf, sub), sign))
        acc += b[1]
    return out


def density_T1_insert(G, t, lam, ins):
    """(coeff, U): D = coeff U^{-1/2} for U^{-H-1/2} <prod_{e in ins} w_e prod psi^tau> ^ prod psi^i / iota_E vol."""
    t = [Fr(x) for x in t]
    C = G.cycle_matrix()
    n, N = len(C), G.ne
    M = [[sum((Fr(C[a][e] * C[b][e]) * t[e] for e in range(N)), Fr(0)) for b in range(n)] for a in range(n)]
    Minv = inverse(M) if n else []
    K = [[sum((C[a][e] * Minv[a][b] * C[b][f] for a in range(n) for b in range(n)), Fr(0)) for f in range(N)]
         for e in range(N)]
    cov = [[(((1 / t[e]) if e == f else Fr(0)) - K[e][f]) / 2 for f in range(N)] for e in range(N)]
    items = [({e: {0: Fr(1)}}, 0) for e in ins]
    items += [({e: {1 << e: Fr(C[r][e])} for e in range(N) if C[r][e]}, 1) for r in range(n)]
    top = wick(items, cov)
    y = brute.currents(G, t, lam)
    form = wedge(top, wedge_all(brute.psi_forms(G, t, y, C)))
    U = brute.symanzik_U(G, t, C)
    H = len(y[0]) if y else 0
    coeff = basic_coefficient(form, t, None)
    return coeff / U ** H, U


def t1_checks(cases, npts=3, seed=2, log=print):
    rng = random.Random(seed)
    tot = ok = 0
    signs_all = set()
    for D, n, idx, G, S0, plan in cases:
        H = D - 1
        red = PartialReduction(G, plan)
        el = red.el
        Psteps = [s for s in el.steps if s.kind == "P"]
        estar = [s.new[0] for s in Psteps]
        assert all(e in red.core_edges for e in estar) and all(len(s.edges) == 2 for s in Psteps)
        lam = rand_sources(G, S0, H, rng)
        lam0, _ = lam_dict(G, lam)
        xcoords = [b for s in el.steps if s.kind == "S" for b in s.block[1:]]
        signs = set()
        good = True
        for _ in range(npts):
            xi = rand_point(red, rng)
            Tc = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(red.core.ne)]
            k, Nc = el.k, red.core.ne
            nvar = k + Nc
            X = {c: Dual.var(xi[c], j, nvar) for j, c in enumerate(el.coords)}
            one = Dual.const(1, nvar)
            bv = el.bary_values(X, one)
            t, _ = el.times(bv, {e: Dual.var(Tc[i], k + i, nvar) for i, e in enumerate(red.core_edges)})
            J_aff = fdet([[x.g[j] for j in range(nvar)] for x in t])
            tv = [x.v for x in t]
            Dg = brute.density(G, tv, lam, T=1)
            lhs_c, lhs_U = Dg.coeff * J_aff * (-1) ** k, Dg.radicand
            # core with insertion
            bvn = el.bary_values({c: Fr(xi[c]) for c in el.coords}, Fr(1))
            _, talln = el.times(bvn, {e: Fr(Tc[i]) for i, e in enumerate(red.core_edges)})
            lamc = el.sources(lam0, bvn, talln, None)
            lam_core = [lamc[red.inv[i]] for i in range(red.core.nv - 1)]
            ins = [red.core_edges.index(e) for e in estar]
            cc, cU = density_T1_insert(red.core, [Fr(x) for x in Tc], lam_core, ins)
            # F = prod (-sqrt(T_e*) w_e*): rational part (-1)^m, radicand prod T_e*
            rad = cU
            rat = cc * (-1) ** len(estar)
            for e in ins:
                rad = rad / Fr(Tc[e])
            # prod_m (w_m1 w_m2)^{-1/2}
            for s in Psteps:
                rad = rad * bvn[s.block[0]] * bvn[s.block[1]]
            # J_x by duals over (x, ell_core)
            Jx = _Jx(red, lam, xi, Tc, xcoords, H)
            rhs_c = rat * Jx
            # compare sign and square:  lhs_c / sqrt(lhs_U)  vs  s * rhs_c / sqrt(rad)
            if lhs_c == 0 and rhs_c == 0:
                continue
            tot += 1
            sq_ok = lhs_c ** 2 / lhs_U == rhs_c ** 2 / rad
            sg = (1 if lhs_c > 0 else -1) * (1 if rhs_c > 0 else -1)
            signs.add(sg)
            if sq_ok:
                ok += 1
            else:
                good = False
        signs_all |= {(len(signs) == 1, s) for s in signs}
        log("T=1 D=%d (H=%d) n=%d #%d plan %s core n_c=%d insert %d edge(s): squares %s, relative sign %s" % (
            D, H, n, idx, plan, red.core.loops, len(estar), "exact" if good else "FAIL", sorted(signs)))
    return tot, ok, signs_all


def _Jx(red, lam, xi, Tc, xcoords, H):
    """det d ell_Gamma / d(x, ell_core) with the parallel fractions held fixed."""
    G, el, core = red.G, red.el, red.core
    lam0, _ = lam_dict(G, lam)
    nx_ = len(xcoords)
    nc_ = core.loops
    nv2 = nx_ + nc_ * H
    X2 = {}
    for c in el.coords:
        if c in xcoords:
            X2[c] = Dual.var(xi[c], xcoords.index(c), nv2)
        else:
            X2[c] = Dual.const(xi[c], nv2)
    one2 = Dual.const(1, nv2)
    bv2 = el.bary_values(X2, one2)
    _, tall2 = el.times(bv2, {e: one2 * Fr(Tc[i]) for i, e in enumerate(red.core_edges)})
    lam2 = el.sources({v: [one2 * x for x in vec] for v, vec in lam0.items()}, bv2, tall2, None)
    lamcore2 = [lam2[red.inv[i]] for i in range(core.nv - 1)]
    tree = core.spanning_tree()
    cot = core.cotree(tree)
    B = core.incidence()
    BTinv = inverse([[Fr(B[e][v]) for e in tree] for v in range(core.nv - 1)])
    C = core.fundamental_cycles(tree)
    lam_num = [[x.v for x in vec] for vec in lamcore2]
    yc_num = brute.currents(core, [Fr(x) for x in Tc], lam_num)
    ycore = {}
    for pos, e in enumerate(tree):
        ycore[e] = [sum((BTinv[pos][v] * lamcore2[v][i] for v in range(core.nv - 1) if BTinv[pos][v] != 0),
                        Dual.const(0, nv2)) for i in range(H)]
    for e in cot:
        ycore[e] = [Dual.const(0, nv2) for _ in range(H)]
    for r, c in enumerate(cot):
        for i in range(H):
            lv = Dual.var(yc_num[c][i], nx_ + r * H + i, nv2)
            for e in range(core.ne):
                if C[r][e]:
                    ycore[e][i] = ycore[e][i] + C[r][e] * lv
    yG = el.currents({red.core_edges[e]: ycore[e] for e in range(core.ne)}, bv2, tall2)
    cotG = G.cotree()
    return fdet([[yG[c][i].g[j] for j in range(nv2)] for c in cotG for i in range(H)])


# ---------------------------------------------------------------------------- case selection

def enlarge_sources(Hx, S0, plan):
    """all vertices except the star centres of the plan, if every centre stays source-free at its
    elimination (Kron source transfer); else S0."""
    centres = [v for k, v in plan if k == "Y"]
    S1 = [v for v in Hx if v not in centres]
    net = Network(Hx, S1)
    X = 0
    for k, v in plan:
        nb, src = net.state(X)
        if k == "Y" and src >> v & 1:
            return S0
        X |= 1 << v
    return S1


def tight_cases(D, n, maxcases, want_kinds=("chain", "single", "paired")):
    """non-SP graphs with a nonempty tight partial reduction to a core with >= 1 loop."""
    out = []
    for idx, Gx in enumerate(admissible(D, n)):
        if catalog.is_sp(Gx):
            continue
        G, _ = from_networkx(Gx)
        Hx = nxg(G)
        S0 = sorted(set(v for v in Hx if Hx.degree(v) == 2) | {G.root})
        ts = tight_search(Hx, S0, D)
        path = ts["path"] if ts["reducible"] else ts["stuck_path"]
        if not path:
            continue
        # truncate to a partial plan leaving at least one loop
        steps = []
        kinds = set()
        i = 0
        while i < len(path):
            kind, v = path[i]
            if kind == "Y":
                nxt = path[i + 1] if i + 1 < len(path) else None
                kinds.add("Y")
            steps.append(path[i])
            i += 1
        red_full_loops = None
        # find the longest prefix leaving a core with >= 1 loop (prefixes at move boundaries)
        for cut in range(len(steps), 0, -1):
            try:
                red = PartialReduction(G, steps[:cut])
            except ValueError:
                continue
            if red.core.loops >= 1 and red.core.ne == D * red.core.loops + 1:
                kind = "mixed" if any(k == "Y" for k, _ in steps[:cut]) else "chain"
                out.append((D, n, idx, G, enlarge_sources(Hx, S0, steps[:cut]), steps[:cut], kind))
                break
        if len(out) >= maxcases:
            break
    return out


if __name__ == "__main__":
    logf = open("check_factorization.log", "w")

    def log(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        logf.write(s + "\n")
        logf.flush()

    t0 = time.time()
    log("== T = 0: tight partial reductions, H = D ==")
    cases = []
    for D, n, m in [(2, 4, 8), (2, 5, 12), (3, 4, 10), (4, 3, 0)]:
        cases += tight_cases(D, n, m)
    # dedicated LLC single / paired moves at D = 2: first move only
    singles, pairs = [], []
    for n in (4, 5):
        for idx, Gx in enumerate(admissible(2, n)):
            G, _ = from_networkx(Gx)
            Hx = nxg(G)
            S0 = sorted(set(v for v in Hx if Hx.degree(v) == 2) | {G.root})
            for v in range(G.nv):
                nb = set(Hx.neighbors(v))
                if v in S0 or len(nb) != 3:
                    continue
                for u in nb:
                    for w in nb:
                        if u < w:
                            x = (nb - {u, w}).pop()
                            if not Hx.has_edge(u, w) and Hx.has_edge(u, x) and Hx.has_edge(w, x):
                                if len(singles) < 8:
                                    singles.append((2, n, idx, G, [a for a in range(G.nv) if a != v], [("Y", v)], "single"))
            for v1 in range(G.nv):
                for v2 in range(v1 + 1, G.nv):
                    if v1 in S0 or v2 in S0:
                        continue
                    N1, N2 = set(Hx.neighbors(v1)), set(Hx.neighbors(v2))
                    if len(N1) == 3 and N1 == N2:
                        for u in N1:
                            A, B = sorted(N1 - {u})
                            if Hx.has_edge(A, B) and not Hx.has_edge(u, A) and not Hx.has_edge(u, B):
                                if len(pairs) < 8:
                                    pairs.append((2, n, idx, G, [a for a in range(G.nv) if a not in (v1, v2)], [("Y", v1), ("Y", v2)], "paired"))
    tot, ok, ls, lp, nz = t0_checks(cases + singles + pairs, lambda D: D, log=log)
    log("T=0 summary: %d/%d points exact (%d with nonzero density) over %d reductions" % (ok, tot, nz, len(cases + singles + pairs)))
    log("LLC b-01159 single-star factor: J_poly / (r1 r2 det(y_ux,y_wx)) values %s" % sorted(set(ls)))
    log("LLC b-01173 paired factor: J_poly / ((r1-r2)(1-r3-r4)<y_uA,y_AB><y_uB,y_AB> g_m) values %s" % sorted(set(lp)))

    log("== T = 1: insertion-free Wick density = htk.brute.density(T=1) ==")
    rng = random.Random(9)
    nchk = 0
    for D, n in [(2, 2), (2, 4), (3, 2), (3, 3)]:
        for Gx in admissible(D, n)[:6]:
            G, _ = from_networkx(Gx)
            lam = [[rng.randint(-5, 5) for _ in range(D - 1)] for _ in range(G.nv - 1)]
            t = [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(G.ne)]
            c, U = density_T1_insert(G, t, lam, [])
            Db = brute.density(G, t, lam, T=1)
            assert Db.coeff == c and (c == 0 or Db.radicand == U), (c, U, Db)
            nchk += 1
    log("wick (no insertion) = brute T=1 density on %d graphs: exact" % nchk)
    log("== T = 1: k = 2 tight chain moves, merged edges in the core, H = D - 1 ==")
    t1cases = []
    for D, n in [(2, 4), (3, 4), (2, 6)]:
        cnt = 0
        for idx, Gx in enumerate(admissible(D, n)):
            if catalog.is_sp(Gx) or cnt >= (4 if n < 6 else 2):
                continue
            G, _ = from_networkx(Gx)
            Hx = nxg(G)
            S0 = sorted(set(v for v in Hx if Hx.degree(v) == 2) | {G.root})
            # one k = 2 tight chain move: chains p..q with a1 + a2 = D + 1
            found = None
            for p in range(G.nv):
                for q in range(p + 1, G.nv):
                    from reduce_search import _chains_between
                    adjm = {v: set(Hx.neighbors(v)) for v in Hx}
                    ch = _chains_between(adjm, p, q)
                    for c1 in range(len(ch)):
                        for c2 in range(c1 + 1, len(ch)):
                            if len(ch[c1]) + len(ch[c2]) + 2 == D + 1:
                                found = [("S", v) for v in ch[c1] + ch[c2]]
                                break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                red = PartialReduction(G, found)
                if not catalog.is_sp(nxg(red.core)):
                    t1cases.append((D, n, idx, G, list(range(G.nv)), found))
                    cnt += 1
    tot1, ok1, sg = t1_checks(t1cases, log=log)
    log("T=1 summary: %d/%d points exact (sign and square) over %d reductions; per-reduction sign constant: %s" % (
        ok1, tot1, len(t1cases), sorted(sg)))
    log("time %.1fs" % (time.time() - t0))
