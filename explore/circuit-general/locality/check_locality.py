"""Pointwise exact checks of the locality factorization against htk.brute (T = 0, 1).
(A) form level, mode R: Omega_Gamma = eps * Omega_S ^ pi^* Omega_out   (any S, any k)
(B) form level, mode P: same without the dh terms (k = 2, tight S)
(C) chart level, k = 2: D_Gamma J_Phi = kappa D_S(t_S; sigma, lambda_S0) D_Gamma'(t'; lambda'(t_S))
eps, kappa must be constant (+-1) over points and sources."""
import random
from fractions import Fraction as Fr
import sys
from fractions import Fraction as Fr

import locality as L
from htk.graph import Graph


def k4sub(chains, H, order=None):
    """K4 on 0..3 with skeleton edges 01,02,03,12,13,23 subdivided into chains of the given
    lengths; returns (Graph, chain edge lists).  Subdivision vertices numbered from 4."""
    sk = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    nv = 4
    edges, ch = [], []
    for (a, b), ln in zip(sk, chains):
        path = [a] + list(range(nv, nv + ln - 1)) + [b]
        nv += ln - 1
        ids = []
        for u, v in zip(path, path[1:]):
            ids.append(len(edges))
            edges.append((u, v))
        ch.append(ids)
    return nv, edges, ch


def relabel_root(nv, edges, root):
    """move vertex `root` to the last label"""
    perm = [v for v in range(nv) if v != root] + [root]
    new = {v: i for i, v in enumerate(perm)}
    return Graph(nv, [(new[a], new[b]) for a, b in edges]), new


def rand_lam(G, H, rng, sourced=None):
    lam = []
    for v in range(G.nv - 1):
        if sourced is None or v in sourced:
            lam.append([rng.randint(-4, 4) for _ in range(H)])
        else:
            lam.append([0] * H)
    return lam


def run_form(name, G, S, H, T, modes, rng, npts=3, sourced=None):
    sp = L.Split(G, S)
    res = []
    for mode in modes:
        ratios = set()
        for _ in range(npts):
            t = L.rand_times(G.ne, rng)
            lam = rand_lam(G, H, rng, sourced)
            if T == 0:
                Db, Df, info = L.factor_T0(sp, t, lam, mode)
                if Db == 0 and Df == 0:
                    ratios.add('0=0')
                    continue
                ratios.add(Db / Df if Df else 'Df=0,Db=%s' % Db)
            else:
                Db, c, U = L.factor_T1(sp, t, lam, mode)
                if Db.coeff == 0 and c == 0:
                    ratios.add('0=0')
                    continue
                r = L.same_density(Db, c, U)
                ratios.add(r if r is not None else 'MISMATCH %s vs %s/%s' % (Db, c, U))
        res.append((mode, ratios))
    nz = sum(1 for _ in range(1))
    print('%-34s H=%d T=%d k=%d m=%d nS=%d N_S=%d | %s' % (
        name, H, T, sp.k, sp.m, sp.nS, len(sp.S),
        '  '.join('%s: %s' % (mo, sorted(map(str, r))) for mo, r in res)), flush=True)
    return res


def run_chart(name, G, S, H, T, rng, npts=3, sourced=None):
    tt = L.TwoTerminal(G, S)
    out = set()
    for _ in range(npts):
        tS = L.rand_times(len(tt.sp.S), rng)
        tp = L.rand_times(len(tt.sp.O) + 1, rng)
        lam = rand_lam(G, H, rng, sourced)
        DG, JPhi, DS, Dp, sigma = tt.check_point(tS, tp, lam, T=T)
        lhs = DG.coeff * JPhi
        rhs = DS.coeff * Dp.coeff
        if lhs == 0 and rhs == 0:
            out.add('0=0')
        else:
            out.add(lhs / rhs if rhs else 'rhs=0')
    print('%-34s chart k=2 H=%d T=%d: D_G J / (D_S D\') = %s' % (name, H, T, sorted(map(str, out))), flush=True)
    return out


if __name__ == '__main__':
    rng = random.Random(20261123)
    print('== T = 0, form level (R: response differentials dZ, dh; P: without dh)')
    # D=2 K4 with one subdivided edge: bridge S = K4 - {01}, O = path 0-4-1
    nv, ed, ch = k4sub([2, 1, 1, 1, 1, 1], 2)
    G, new = relabel_root(nv, ed, 4)
    S = [e for c in ch[1:] for e in c]
    run_form('D2 K4[2,1,1,1,1,1] bridge', G, S, 2, 0, ['R', 'P'], rng)
    run_chart('D2 K4[2,1,1,1,1,1] bridge', G, S, 2, 0, rng)
    # D=3 n=3: K4 chain lengths, chain 01 of length 3 outside, S tight (7 edges)
    for chains in ([3, 2, 2, 1, 1, 1], [3, 1, 1, 2, 2, 1], [3, 2, 1, 1, 2, 1], [3, 1, 2, 1, 1, 2]):
        nv, ed, ch = k4sub(chains, 3)
        G, new = relabel_root(nv, ed, 4)
        S = [e for c in ch[1:] for e in c]
        if G.root in [G.edges[e][0] for e in S] + [G.edges[e][1] for e in S]:
            pass
        try:
            run_form('D3 K4%s bridge' % chains, G, S, 3, 0, ['R', 'P'], rng, npts=2)
            run_chart('D3 K4%s bridge' % chains, G, S, 3, 0, rng, npts=2)
        except ValueError as ex:
            print('  skip', chains, ex)

    print('== T = 1, form level and chart level, k = 2')
    nv, ed, ch = k4sub([2, 1, 1, 1, 1, 1], 2)
    G, new = relabel_root(nv, ed, 4)
    S = [e for c in ch[1:] for e in c]
    run_form('(1,1) K4[2,1,1,1,1,1] bridge', G, S, 1, 1, ['R', 'P'], rng, npts=3)
    for chains in ([3, 2, 1, 1, 2, 1], [3, 1, 2, 1, 1, 2]):
        nv, ed, ch = k4sub(chains, 3)
        G, new = relabel_root(nv, ed, 4)
        S = [e for c in ch[1:] for e in c]
        run_form('(2,1) K4%s bridge' % chains, G, S, 2, 1, ['R', 'P'], rng, npts=2)

    print('== random subnetworks of census graphs (S = edges at a vertex set W, root not in W)')
    from htk import catalog
    stats = {}
    for (H, T, n, nsamp) in ((2, 0, 3, 3), (2, 0, 4, 2), (3, 0, 3, 2), (1, 1, 3, 2), (1, 1, 4, 1), (2, 1, 3, 1)):
        Ds = H + T
        gl = catalog.tight_graphs(Ds, n)
        for gi, G in enumerate(gl):
            for _ in range(nsamp):
                cand = [v for v in range(G.nv - 1)]
                W = set(rng.sample(cand, rng.randint(1, min(3, len(cand)))))
                S = [e for e, (a, b) in enumerate(G.edges) if a in W or b in W]
                if len(S) == G.ne:
                    continue
                try:
                    sp = L.Split(G, S)
                except ValueError:
                    continue
                if sp.k < 2:
                    continue
                tight = (len(S) - 1 == sp.nS * Ds) if sp.k == 2 else None
                modes = ['R'] + (['P'] if sp.k == 2 else [])
                for mode in modes:
                    t = L.rand_times(G.ne, rng)
                    lam = rand_lam(G, H, rng)
                    if T == 0:
                        Db, Df, _ = L.factor_T0(sp, t, lam, mode)
                        ok = 'zero' if Db == 0 and Df == 0 else ('+1' if Db == Df else ('-1' if Db == -Df else 'FAIL'))
                    else:
                        Db, c, U = L.factor_T1(sp, t, lam, mode)
                        r = L.same_density(Db, c, U)
                        ok = 'zero' if Db.coeff == 0 and c == 0 else {1: '+1', -1: '-1', None: 'FAIL'}[r]
                    key = ((H, T), 'k=%d' % sp.k, mode, 'tight' if tight else ('nontight' if tight is False else '-'))
                    stats.setdefault(key, {}).setdefault(ok, 0)
                    stats[key][ok] += 1
    for key in sorted(stats):
        print('  (H,T)=%s %s mode %s %-8s : %s' % (key[0], key[1], key[2], key[3], stats[key]))

    print('== currents and displacement map, k = 2 (chart points): y_O = y\'_O, y_S = S standalone currents (lambda_S0, sigma)')
    from htk import brute as _b
    okc = []
    for chains in ([3, 2, 1, 1, 2, 1], [3, 1, 2, 1, 1, 2]):
        nv, ed, ch = k4sub(chains, 3)
        G, new = relabel_root(nv, ed, 4)
        S = [e for c in ch[1:] for e in c]
        tt = L.TwoTerminal(G, S)
        for _ in range(2):
            tS = L.rand_times(len(S), rng)
            tp = L.rand_times(len(tt.sp.O) + 1, rng)
            lam = rand_lam(G, 3, rng)
            t = tt.point(tS, tp, lam)
            y = _b.currents(G, t, lam)
            lamp, s, H = tt.sources_prime(tS, lam)
            yp = _b.currents(tt.Gp, tp, lamp)
            sigma = [yp[-1][i] - s[i] for i in range(H)]
            lf, _ = L.lam_full(G, lam)
            lamS = [lf[v] for v in tt.sp.S0] + [sigma]
            c = tt.point(tS, tp, lam)[tt.sp.S[0]] / Fr(tS[0])
            yS = _b.currents(tt.GS, [Fr(x) * c for x in tS], lamS)
            okc.append(all(y[e] == yp[j] for j, e in enumerate(tt.sp.O)) and
                       all(y[e] == yS[j] for j, e in enumerate(tt.sp.S)))
    print('  ', okc)

    print('== k = 3: det M_out * P_Delta = U(Gamma\' = O + triangle, rho_ab = 1/c_ab)')
    from htk import catalog as _cat
    okd = []
    for G in _cat.tight_graphs(2, 4)[:10]:
        for v in range(G.nv - 1):
            if G.degree(v) != 3:
                continue
            S = [e for e, (a, b) in enumerate(G.edges) if v in (a, b)]
            sp_ = L.Split(G, S)
            t = L.rand_times(G.ne, rng)
            lam = rand_lam(G, 2, rng)
            Db, Df, info = L.factor_T0(sp_, t, lam, 'R')
            K = sp_.kron(t)
            rho = {}
            for a in range(3):
                for b in range(a + 1, 3):
                    rho[(a, b)] = -1 / K['Lr'][a][b]
            vp = sorted(set(sp_.VO) | set(sp_.dS))
            idx = {x: i for i, x in enumerate(vp)}
            edges = [(idx[G.edges[e][0]], idx[G.edges[e][1]]) for e in sp_.O]
            tt_ = [Fr(t[e]) for e in sp_.O]
            for (a, b), r in rho.items():
                edges.append((idx[sp_.dS[a]], idx[sp_.dS[b]]))
                tt_.append(r)
            Gp = Graph(len(vp), edges)
            Up = _b.symanzik_U(Gp, tt_)
            okd.append(info['Uout'] * sum(rho.values()) == Up)
    print('  ', okd)
