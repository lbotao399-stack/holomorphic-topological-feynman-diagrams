"""Verification of a1, a2, a5 (rho part), a7, a8 with the own exterior-algebra layer (vcore), exact.
 A1  Berezin form: [(-1)^n/n! sum_{r_1..r_n} E[a_r1..a_rn] Z_r1...Z_rn]_{theta_1..theta_n} == Pf(-1/2 dM M^-1 dM),
     Z_r = sum_s dM_rs theta_s, E = Gaussian moments with covariance M^{-1}/2 (Wick), i.e. the a-integral of
     exp(-a^T M a - sum_e A_e dt_e Theta_e) done by moments, not by completing the square.
 A2  a plain chain {e,f}: Pf depends on (t_e, t_f) only through t_e + t_f and on (dt_e, dt_f) only through
     dt_e + dt_f (components compared exactly); in own ratio coordinates no dx-component for plain chains.
 A5  plain bananas, C_rec rows e_a - e_{n+1}, t_i = 1/w_i: Pf = ((n-1)!!/2^{n/2}) / prod w  (i.e.
     rho = ((n-1)!!/2^{n/2}) prod w^{-1/2} dw_1..dw_n), n = 2, 4, 6; non-plain banana (2,2,3,3,3) in
     (w, x) coordinates; int prod w^{-1/2} * (n-1)!!/2^{n/2} over the simplex = pi^{n/2} (own Gamma).
 A7  G2 small plan, C_rec built from the stated rule: Pf == -prod_{m=1}^3 [t_m dw_m1/(w_m1 w_m2)] ^ (1/4T) dlog(t9 t10 t11/T^3)
     as 4-forms in the 9 dt_e, at random rational t (no gauge); U == prod_m t_m/(w_m1 w_m2); second form at t_4 = 1 (float).
 A8  rho ^ rho = 0: Sym_4 (10 generators), Sym_6 restricted to a random 12-dim family, n=6 SP graph of d=2
     pulled back (12 generators, top form), nested (2,1) graph n = 4 (12 generators); tr (dM M^-1)^{2j} = 0.
"""
import sys, random, time, math, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/ht-families/verify")
import vcore as vc
import vsp
import vdir

rng = random.Random(31)


def slice_dts(N):
    return [vc.oneform([1 if j == e else 0 for j in range(N - 1)]) if e < N - 1 else {} for e in range(N)]


def full_dts(N):
    return [vc.oneform([1 if j == e else 0 for j in range(N)]) for e in range(N)]


def wick(idx, S):
    if not idx:
        return Fr(1)
    if len(idx) % 2:
        return Fr(0)
    a = idx[0]
    tot = Fr(0)
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        if S[a][idx[j]] != 0:
            tot += S[a][idx[j]] * wick(rest, S)
    return tot


def A1(C, t, dts, G):
    n, N = len(C), len(t)
    M = [[sum(C[r][e] * C[s][e] * Fr(t[e]) for e in range(N)) for s in range(n)] for r in range(n)]
    Mi = vc.inv_fr(M)
    S = [[x / 2 for x in row] for row in Mi]
    dM = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            f = {}
            for e in range(N):
                if C[r][e] * C[s][e]:
                    f = vc.fadd(f, dts[e], C[r][e] * C[s][e])
            dM[r][s] = f
    Z = []
    for r in range(n):
        z = {}
        for s in range(n):
            z = vc.fadd(z, vc.wedge(dM[r][s], {1 << (G + s): Fr(1)}))
        Z.append(z)
    tot = {}
    for idx in itertools.product(range(n), repeat=n):
        w = wick(list(idx), S)
        if w == 0:
            continue
        p = {0: Fr(1)}
        for r in idx:
            p = vc.wedge(p, Z[r])
        tot = vc.fadd(tot, p, w)
    tot = vc.fscale(tot, Fr((-1) ** n, math.factorial(n)))
    thmask = ((1 << n) - 1) << G
    out = {}
    for k, v in tot.items():
        if k & thmask == thmask:
            out[k & ~thmask] = out.get(k & ~thmask, 0) + v
        elif k & thmask:
            pass
    pf, U, _ = vc.rho_forms(C, t, dts)
    return {k: v for k, v in out.items() if v != 0} == pf


G2E = [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), (3, 4), (4, 5)]
G4E = [(0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (2, 3), (3, 4), (4, 5)]
GALL = {
    0: [(0, 4), (0, 5), (1, 4), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)],
    1: [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (2, 5), (3, 4)],
    2: G2E, 3: [(0, 1), (0, 2), (0, 5), (1, 2), (1, 3), (1, 4), (1, 5), (3, 4), (4, 5)], 4: G4E}


def rand_t(N):
    return [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(N)]


if __name__ == "__main__":
    # ---------------- A1
    from v_b2 import banana, nested
    t0 = time.time()
    for name, nv, edges in [("theta(1,2,2)", 4, [(0, 3), (0, 1), (1, 3), (0, 2), (2, 3)]), ("G2", 6, G2E), ("G4", 6, G4E),
                            ("banana(2,2,3,3,3)",) + banana((2, 2, 3, 3, 3))[:2]]:
        net = vc.Net(nv, edges)
        C = net.cycle_basis()
        res = [A1(C, rand_t(net.N), slice_dts(net.N), net.N - 1) for _ in range(2)]
        print("A1 %-18s n=%d: Berezin moments == Pf at 2 rational points: %s" % (name, net.n, res))
    print("   (%.1fs)" % (time.time() - t0))

    # ---------------- A2
    for k, edges in GALL.items():
        net = vc.Net(6, edges)
        C = net.cycle_basis()
        deg = {v: [e for e, ab in enumerate(edges) if v in ab] for v in range(6)}
        chains = [tuple(deg[v]) for v in range(6) if len(deg[v]) == 2]
        ok = True
        for (e, f) in chains:
            same_or = (edges[e][1] == edges[f][0]) or (edges[e][0] == edges[f][1])
            # columns of C equal up to the orientation sign
            sg = 1 if all(C[r][e] == C[r][f] for r in range(len(C))) else (-1 if all(C[r][e] == -C[r][f] for r in range(len(C))) else 0)
            ok &= sg != 0
            t = rand_t(net.N)
            d = Fr(rng.randint(1, 9), 97) * min(t[e], t[f])
            t2 = list(t)
            t2[e] += d
            t2[f] -= d
            dts = full_dts(net.N)
            p1, U1, _ = vc.rho_forms(C, t, dts)
            p2, U2, _ = vc.rho_forms(C, t2, dts)
            ok &= p1 == p2 and U1 == U2
            # dependence on (dt_e, dt_f) only through dt_e + dt_f (C_re C_se = C_rf C_sf):
            # (dt_e, dt_f) -> (0, dt_e + dt_f) reproduces Pf
            dts3b = list(dts)
            dts3b[e] = {}
            dts3b[f] = vc.fadd(dts[f], dts[e])
            p3, _, _ = vc.rho_forms(C, t, dts3b)
            ok &= p3 == p1
        # own ratio coordinates: no dx component for plain chains, coefficients independent of x
        rc = vsp.RatioCoords(vsp.sp_decompose(6, edges), net.N)
        plainx = [b for nd, nm in rc.blocks if nd.kind == 'S' and all(c.kind == 'E' for c in nd.children) for b in nm[1:]]
        xi = rc.random_point(rng)
        full = rc.complete(xi)
        tt = rc.t_of(full)
        rows = rc.jac_rows(full, tt)
        dts = [vc.oneform([rows[c][e] for c in range(net.N - 1)]) for e in range(net.N)]
        pf, U, _ = vc.rho_forms(C, tt, dts)
        bad = [kk for kk in pf if any(kk >> rc.coords.index(b) & 1 for b in plainx)]
        # vary the x of plain chains
        xi2 = dict(xi)
        for nd, nm in rc.blocks:
            if nd.kind == 'S' and all(c.kind == 'E' for c in nd.children):
                parts = [Fr(rng.randint(1, 20)) for _ in nm]
                s = sum(parts)
                for b, p in zip(nm[1:], parts[1:]):
                    xi2[b] = p / s
        full2 = rc.complete(xi2)
        tt2 = rc.t_of(full2)
        rows2 = rc.jac_rows(full2, tt2)
        dts2 = [vc.oneform([rows2[c][e] for c in range(net.N - 1)]) for e in range(net.N)]
        pf2, U2, _ = vc.rho_forms(C, tt2, dts2)
        print("A2 G%d: plain chains %s: Pf(t) depends on chain only via (R, dR): %s ; ratio coords: #components with dx(plain) = %d ; "
              "Pf, U unchanged under x(plain) change: %s" % (k, chains, ok, len(bad), pf == pf2 and U == U2))

    # ---------------- A5
    for n in (2, 4, 6):
        k = n + 1
        edges = [(0, 1)] * k
        C = [[1 if e == a else (-1 if e == k - 1 else 0) for e in range(k)] for a in range(n)]
        okall = True
        for _ in range(2):
            parts = [Fr(rng.randint(1, 30)) for _ in range(k)]
            s = sum(parts)
            w = [p / s for p in parts]
            t = [1 / x for x in w]
            dts = []
            for i in range(k):
                if i < n:
                    dts.append(vc.oneform([-1 / w[i] ** 2 if j == i else 0 for j in range(n)]))
                else:
                    dts.append(vc.oneform([1 / w[k - 1] ** 2] * n))
            pf, U, _ = vc.rho_forms(C, t, dts)
            c = Fr(math.prod(range(n - 1, 0, -2)), 2 ** (n // 2))
            pw = math.prod(w)
            okall &= pf == {(1 << n) - 1: c / pw} and U == 1 / pw
        q, kk = vdir.dirichlet([Fr(-1, 2)] * k)
        val = c * q
        print("A5 plain banana n=%d: Pf == ((n-1)!!/2^{n/2})/prod w, U == 1/prod w (2 pts): %s ; ((n-1)!!/2^{n/2}) int prod w^{-1/2} = %s * pi^%s" % (
            n, okall, val, Fr(kk, 2)))
    # non-plain banana (2,2,3,3,3), coordinates (w_1..w_4; x of each chain), C_rec rows chain a - chain 5
    lens = (2, 2, 3, 3, 3)
    nv, edges, ints = banana(lens)
    k = 5
    chains, e0 = [], 0
    for a in lens:
        chains.append(list(range(e0, e0 + a)))
        e0 += a
    C = [[(1 if e in chains[a] else 0) - (1 if e in chains[k - 1] else 0) for e in range(13)] for a in range(4)]
    names = ["w%d" % i for i in range(1, 6)] + ["x%d_%d" % (i, j) for i in range(5) for j in range(lens[i])]
    coords = ["w%d" % i for i in range(1, 5)] + ["x%d_%d" % (i, j) for i in range(5) for j in range(1, lens[i])]
    elim = {c: ("w5" if c.startswith("w") else c.split("_")[0] + "_0") for c in coords}
    texp = []
    for i, ch in enumerate(chains):
        for j, e in enumerate(ch):
            texp.append({"w%d" % (i + 1): -1, "x%d_%d" % (i, j): 1})
    ok = True
    for _ in range(2):
        val = {}
        parts = [Fr(rng.randint(1, 30)) for _ in range(5)]
        for i in range(5):
            val["w%d" % (i + 1)] = parts[i] / sum(parts)
        for i in range(5):
            pp = [Fr(rng.randint(1, 30)) for _ in range(lens[i])]
            for j in range(lens[i]):
                val["x%d_%d" % (i, j)] = pp[j] / sum(pp)
        t = [math.prod(val[b] ** ex for b, ex in te.items()) for te in texp]
        dts = []
        for e in range(13):
            vec = [t[e] * (Fr(texp[e].get(c, 0)) / val[c] - Fr(texp[e].get(elim[c], 0)) / val[elim[c]]) for c in coords]
            dts.append(vc.oneform(vec))
        pf, U, _ = vc.rho_forms(C, t, dts)
        pw = math.prod(val["w%d" % i] for i in range(1, 6))
        ok &= pf == {(1 << 4) - 1: Fr(3, 4) / pw} and U == 1 / pw
    print("A5 banana(2,2,3,3,3) in (w, x) coordinates: Pf == (3/4)/prod w dw_1..dw_4 (no dx), U == 1/prod w (2 pts): %s" % ok)

    # ---------------- A7
    plan = [(0, 1, [[0], [1, 4]], [[1], [1, -1]]), (0, 4, [[2], [3, 8]], [[1], [1, -1]]),
            (1, 4, [[5, 7], [6]], [[1, 1], [1]]), (0, 1, [[9], [10, 11]], [[1], [1, -1]])]
    N0 = 9
    flow = {e: {e: 1} for e in range(N0)}
    rows_by = []
    for mi, (p, q, chs, sgs) in enumerate(plan):
        def cf(ch, sg):
            f = {}
            for e, s in zip(ch, sg):
                for ee, c in flow[e].items():
                    f[ee] = f.get(ee, 0) + s * c
            return f
        fk = cf(chs[-1], sgs[-1])
        rr = []
        for a in range(len(chs) - 1):
            fa = cf(chs[a], sgs[a])
            rr.append([fa.get(e, 0) - fk.get(e, 0) for e in range(N0)])
        rows_by.append(rr)
        flow[N0 + mi] = fk
    Crec = [r for rr in reversed(rows_by) for r in rr]
    net = vc.Net(6, G2E)
    B = net.Bmat()
    for r in Crec:
        assert all(sum(B[e][v] * r[e] for e in range(9)) == 0 for v in range(5))
    print("A7 C_rec rows:", Crec)
    okA, okU = True, True
    for _ in range(3):
        t = rand_t(9)
        dts = full_dts(9)
        pf, U, _ = vc.rho_forms(Crec, t, dts)
        T = {e: Fr(t[e]) for e in range(9)}
        dT = {e: dts[e] for e in range(9)}
        tm, wm, dwm = {}, {}, {}
        Uf = Fr(1)
        for mi, (p, q, chs, sgs) in enumerate(plan):
            R = [sum(T[e] for e in ch) for ch in chs]
            dR = []
            for ch in chs:
                f = {}
                for e in ch:
                    f = vc.fadd(f, dT[e])
                dR.append(f)
            tr = 1 / sum(1 / x for x in R)
            w = [tr / x for x in R]
            dtr = {}
            for i in range(len(chs)):
                dtr = vc.fadd(dtr, dR[i], w[i] ** 2)
            dw = [vc.fadd(vc.fscale(dtr, 1 / R[i]), dR[i], -tr / R[i] ** 2) for i in range(len(chs))]
            T[N0 + mi], dT[N0 + mi] = tr, dtr
            tm[mi], wm[mi], dwm[mi] = tr, w, dw
            Uf *= tr ** (len(chs) - 1) / math.prod(w)
        rhs = {0: Fr(1)}
        for mi in range(3):
            rhs = vc.wedge(rhs, vc.fscale(dwm[mi][0], tm[mi] / (wm[mi][0] * wm[mi][1])))
        TT = T[9] + T[10] + T[11]
        dlog = {}
        for e in (9, 10, 11):
            dlog = vc.fadd(dlog, dT[e], 1 / T[e])
        dlog = vc.fadd(dlog, vc.fadd(vc.fadd(dT[9], dT[10]), dT[11]), -3 / TT)
        rhs = vc.fscale(vc.wedge(rhs, dlog), Fr(-1) / (4 * TT))
        okA &= pf == rhs and len(pf) > 0
        okU &= U == Uf
    print("A7 G2 small plan: Pf(C_rec) == -prod_m[t_m dw_m1/(w_m1 w_m2)] ^ (1/4T) dlog(t9t10t11/T^3) at 3 points: %s ; U == prod t_m/(w w'): %s" % (okA, okU))
    # second form in the gauge t_4 = 1 (float): rho = -(1/2) prod (w w')^{-1/2} dw_m1 ^ d sqrt(Pi)
    import numpy as np
    t = rand_t(9)
    # rescale t so that t_4 = 1
    T = {e: Fr(t[e]) for e in range(9)}
    for mi, (p, q, chs, sgs) in enumerate(plan):
        R = [sum(T[e] for e in ch) for ch in chs]
        T[N0 + mi] = 1 / sum(1 / x for x in R)
    sc = T[12]
    t = [x / sc for x in t]
    dts = full_dts(9)
    pf, U, _ = vc.rho_forms(Crec, t, dts)
    rho = {kk: float(v) / math.sqrt(float(U)) for kk, v in pf.items()}
    T = {e: Fr(t[e]) for e in range(9)}
    dT = {e: dts[e] for e in range(9)}
    lhs = {0: 1.0}
    for mi, (p, q, chs, sgs) in enumerate(plan):
        R = [sum(T[e] for e in ch) for ch in chs]
        dR = []
        for ch in chs:
            f = {}
            for e in ch:
                f = vc.fadd(f, dT[e])
            dR.append(f)
        tr = 1 / sum(1 / x for x in R)
        w = [tr / x for x in R]
        dtr = {}
        for i in range(len(chs)):
            dtr = vc.fadd(dtr, dR[i], w[i] ** 2)
        dw = [vc.fadd(vc.fscale(dtr, 1 / R[i]), dR[i], -tr / R[i] ** 2) for i in range(len(chs))]
        T[N0 + mi], dT[N0 + mi] = tr, dtr
        if mi < 3:
            lhs = vc.wedge(lhs, {kk: float(v) / math.sqrt(float(w[0] * w[1])) for kk, v in dw[0].items()})
    TT = T[9] + T[10] + T[11]
    Pi = T[9] * T[10] * T[11] / TT ** 3
    dlog = {}
    for e in (9, 10, 11):
        dlog = vc.fadd(dlog, dT[e], 1 / T[e])
    dlog = vc.fadd(dlog, vc.fadd(vc.fadd(dT[9], dT[10]), dT[11]), -3 / TT)
    dsq = {kk: 0.5 * math.sqrt(float(Pi)) * float(v) for kk, v in dlog.items()}
    lhs = vc.wedge(lhs, dsq)
    lhs = {kk: -0.5 * v for kk, v in lhs.items()}
    keys = set(lhs) | set(rho)
    err = max(abs(lhs.get(kk, 0) - rho.get(kk, 0)) for kk in keys) / max(abs(v) for v in rho.values())
    print("A7 second form -(1/2) prod (w w')^{-1/2} dw_m1 ^ d sqrt(Pi) vs rho at t_4 = 1 (float): max rel. deviation %.1e" % err)

    # ---------------- A8
    t0 = time.time()
    # Sym_4: M = M0 + sum_k s_k E_k
    n = 4
    while True:
        Lm = [[Fr(rng.randint(-5, 5)) for _ in range(n)] for _ in range(n)]
        M0 = [[sum(Lm[i][k] * Lm[j][k] for k in range(n)) + (1 if i == j else 0) for j in range(n)] for i in range(n)]
        if vc.det_fr(M0) > 0:
            break
    pairs = [(i, j) for i in range(n) for j in range(i, n)]

    def rho_family(M0, dM):
        n = len(M0)
        Mi = vc.inv_fr(M0)
        A = [[{} for _ in range(n)] for _ in range(n)]
        for r in range(n):
            for s in range(n):
                f = {}
                for p in range(n):
                    for q in range(n):
                        if Mi[p][q] != 0:
                            f = vc.fadd(f, vc.wedge(dM[r][p], dM[q][s]), Mi[p][q])
                A[r][s] = vc.fscale(f, Fr(-1, 2))
        return vc.pfaffian_forms(A)

    dM = [[{} for _ in range(n)] for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        dM[i][j] = vc.fadd(dM[i][j], {1 << k: Fr(1)})
        if i != j:
            dM[j][i] = vc.fadd(dM[j][i], {1 << k: Fr(1)})
    pf = rho_family(M0, dM)
    rr = vc.wedge(pf, pf)
    print("A8 Sym_4 (10 generators): #components of Pf = %d, of Pf^Pf = %d" % (len(pf), len(rr)))
    # Sym_6 on a random 12-dim linear family dM = sum_k E^(k) ds_k with random symmetric integer E^(k)
    n = 6
    while True:
        Lm = [[Fr(rng.randint(-3, 3)) for _ in range(n)] for _ in range(n)]
        M0 = [[sum(Lm[i][k] * Lm[j][k] for k in range(n)) + (1 if i == j else 0) for j in range(n)] for i in range(n)]
        if vc.det_fr(M0) > 0:
            break
    dM = [[{} for _ in range(n)] for _ in range(n)]
    for k in range(12):
        for i in range(n):
            for j in range(i, n):
                c = rng.randint(-3, 3)
                if c:
                    dM[i][j] = vc.fadd(dM[i][j], {1 << k: Fr(c)})
                    if i != j:
                        dM[j][i] = vc.fadd(dM[j][i], {1 << k: Fr(c)})
    pf = rho_family(M0, dM)
    rr = vc.wedge(pf, pf)
    print("A8 Sym_6 random 12-dim family (top form): #components of Pf = %d, of Pf^Pf = %d  (%.1fs)" % (len(pf), len(rr), time.time() - t0))
    # n=6 SP graph at d=2 (pullback to the slice t_12 = 1, 12 generators)
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
    import enumerate_graphs as eg
    import networkx as nx
    gs = [Gx for Gx in eg.enumerate_tight(2, 6) if eg.is_sp(Gx)]
    for gi in (0, 6):
        Gx = gs[gi]
        nodes = sorted(Gx.nodes())
        edges = sorted(tuple(sorted((nodes.index(a), nodes.index(b)))) for a, b in Gx.edges())
        net = vc.Net(len(nodes), edges)
        C = net.cycle_basis()
        pf, U, _ = vc.rho_forms(C, rand_t(net.N), slice_dts(net.N))
        rr = vc.wedge(pf, pf)
        print("A8 d=2 n=6 SP graph %d (N=%d, 12 generators): #components of Pf = %d, of Pf^Pf = %d" % (gi, net.N, len(pf), len(rr)))
    nv, edges = nested((1, 3, 3))
    net = vc.Net(nv, edges)
    C = net.cycle_basis()
    pf, U, _ = vc.rho_forms(C, rand_t(net.N), slice_dts(net.N))
    print("A8 nested theta(1,3,3) (2,1) n=4, 12 generators: #components of Pf = %d, of Pf^Pf = %d" % (len(pf), len(vc.wedge(pf, pf))))
    # traces
    net = vc.Net(6, G2E)
    C = net.cycle_basis()
    t = rand_t(9)
    dts = full_dts(9)
    n = 4
    M = [[sum(C[r][e] * C[s][e] * t[e] for e in range(9)) for s in range(n)] for r in range(n)]
    Mi = vc.inv_fr(M)
    dM = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            f = {}
            for e in range(9):
                if C[r][e] * C[s][e]:
                    f = vc.fadd(f, dts[e], C[r][e] * C[s][e])
            dM[r][s] = f
    Bf = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(n):
            f = {}
            for q in range(n):
                f = vc.fadd(f, vc.fscale(dM[r][q], Mi[q][s]))
            Bf[r][s] = f
    P = [[{0: Fr(int(i == j))} for j in range(n)] for i in range(n)]
    trs = []
    for j in range(1, 5):
        Q = [[{} for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for k2 in range(n):
                f = {}
                for l in range(n):
                    f = vc.fadd(f, vc.wedge(P[i][l], Bf[l][k2]))
                Q[i][k2] = f
        P = Q
        if j % 2 == 0:
            tr = {}
            for i in range(n):
                tr = vc.fadd(tr, P[i][i])
            trs.append(len(tr))
    print("A8 G2: #components of tr B^2, tr B^4 (B = dM M^-1, 9 generators): %s" % trs)
