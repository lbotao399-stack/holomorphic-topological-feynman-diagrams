"""Claim b1: the flow integrand sum_pi s_pi N_pi, re-evaluated by a different algorithm.

W_Gamma = D_{e_1} ... D_{e_13} (W_0 ... W_11),  D_e = sum_i dL/dmu^i_{tail(e)} d/dgamma_{i,head(e)}  (S3 eq-Delta;
the gamma derivative is even).  Each odd D_e acts on the single factor W_{tail(e)}; the W_v are even, so
    W_Gamma = (-1)^{inv(tail(e_1), ..., tail(e_13))} * prod_v [ (prod_{e: tail e = v, increasing} dL/dmu_e)
                                                               (prod_{e: head e = v} d/dgamma_e) W_v ],
inv = number of pairs k<l with tail(e_k) > tail(e_l); the remaining external odd letters then stand in increasing
vertex order.  Each bracket is computed here by a literal left-derivative on the ordered letter word of W_v
(S3 R-action): b: (mu^a, mu^b, d_a gamma_b);  c: (1/6) eps^{ijklm} (gamma_i, d_j gamma_k, d_l gamma_m).
A derivative on the contracted letter at e(0) inserts -y_e, at e(1) +y_e (S3 Section 4.3).
External: mu -> u_v, gamma -> vg_v, d on it -> lam_v.  The 5^13 index sum is contracted along the three paths.

The value is compared with the explorer's implementation flows.flow_integrand (per-letter-map signs s_pi and
12288 separate contractions) on flows of graphs that the explorer did not check pointwise.
"""
import sys, os, random, time
from fractions import Fraction as Fr
from itertools import product, permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, H, kirchhoff, perm_parity_sign, rq, rsimplex

rng = random.Random(4242)


def eps5(vals):
    if len(set(vals)) < 5:
        return 0
    return perm_parity_sign(list(vals))


def words(typ):
    """List of (coefficient, word); word = list of (field, index, derivative indices)."""
    out = []
    if typ == 'b':
        for a in range(H):
            for b in range(H):
                out.append((Fr(1), [('mu', a, ()), ('mu', b, ()), ('ga', b, (a,))]))
    else:
        for (i, j, k, l, m) in permutations(range(H)):
            out.append((Fr(eps5((i, j, k, l, m)), 6), [('ga', i, ()), ('ga', k, (j,)), ('ga', m, (l,))]))
    return out


def apply_derivs(coef, word, ops):
    """ops: list of ('mu', idx, edge) / ('ga', idx, edge), applied right to left.  Returns list of
    (coef, remaining word, {edge: derivative indices of the removed gamma letter})."""
    states = [(coef, word, {})]
    for (f, idx, e) in reversed(ops):
        new = []
        for c, w, dd in states:
            for k, (ff, ii, ders) in enumerate(w):
                if ff != f or ii != idx:
                    continue
                s = 1
                if f == 'mu':
                    nodd = sum(1 for (g, _, _) in w[:k] if g == 'mu')
                    s = -1 if nodd % 2 else 1
                d2 = dict(dd)
                if f == 'ga':
                    d2[e] = ders
                new.append((c * s, w[:k] + w[k + 1:], d2))
        states = new
    return states


def vertex_tensor(v, typ, tails, heads, edges, y, ext):
    """T[(i_e for e in sorted(tails+heads))] for vertex v."""
    inc = sorted(tails + heads)
    T = {}
    W = words(typ)
    for idx in product(range(H), repeat=len(inc)):
        I = dict(zip(inc, idx))
        ops = [('mu', I[e], e) for e in sorted(tails)] + [('ga', I[e], e) for e in heads]
        tot = Fr(0)
        for coef, w in W:
            if coef == 0:
                continue
            for c, rem, dd in apply_derivs(coef, w, ops):
                val = c
                for e, ders in dd.items():
                    sgn = 1 if edges[e][1] == v else -1
                    for d in ders:
                        val *= sgn * y[e][d]
                for (f, ii, ders) in rem:
                    val *= ext[v]['u'][ii] if f == 'mu' else ext[v]['vg'][ii]
                    for d in ders:
                        val *= ext[v]['lam'][d]
                tot += val
        T[idx] = tot
    return inc, T


def grouped_integrand(nv, edges, paths, fl, typ, y, ext):
    tails = {v: [] for v in range(nv)}; heads = {v: [] for v in range(nv)}
    for e, (t, h) in enumerate(fl):
        tails[t].append(e); heads[h].append(e)
    seq = [fl[e][0] for e in range(len(edges))]
    inv = sum(1 for k in range(len(seq)) for l in range(k + 1, len(seq)) if seq[k] > seq[l])
    gsign = -1 if inv % 2 else 1
    Tv = {v: vertex_tensor(v, typ[v], tails[v], heads[v], edges, y, ext) for v in range(nv)}
    # path matrices
    mats = []
    for pv, pe in paths:
        Mx = {(i, j): Fr(int(i == j)) for i in range(H) for j in range(H)}
        for k, v in enumerate(pv):
            inc, T = Tv[v]
            e_prev, e_next = pe[k], pe[k + 1]
            assert inc == sorted([e_prev, e_next])
            Tm = {(a, b): T[(a, b)] for a in range(H) for b in range(H)}   # inc sorted = (e_prev, e_next)
            Mx = {(i, j): sum(Mx[(i, a)] * Tm[(a, j)] for a in range(H)) for i in range(H) for j in range(H)}
        mats.append(Mx)
    incN, TN = Tv[0]
    incR, TR = Tv[nv - 1]
    firsts = [pe[0] for pv, pe in paths]; lasts = [pe[-1] for pv, pe in paths]
    assert incN == sorted(firsts) and incR == sorted(lasts)
    tot = Fr(0)
    for iF in product(range(H), repeat=3):
        dF = dict(zip(firsts, iF))
        tn = TN[tuple(dF[e] for e in incN)]
        if tn == 0:
            continue
        for iL in product(range(H), repeat=3):
            if firsts[0] == lasts[0] and iF[0] != iL[0]:
                continue
            dL = dict(zip(lasts, iL))
            tr = TR[tuple(dL[e] for e in incR)]
            if tr == 0:
                continue
            m = Fr(1)
            for P in range(3):
                if paths[P][1][0] == paths[P][1][-1]:
                    continue                      # single-edge path: identity
                m *= mats[P][(iF[P], iL[P])]
            tot += tn * tr * m
    return gsign * tot


def flows_m1(nv, edges):
    out = []
    deg = [0] * nv
    for a, b in edges:
        deg[a] += 1; deg[b] += 1
    for o in product((0, 1), repeat=len(edges)):
        fl = [(a, b) if o[e] == 0 else (b, a) for e, (a, b) in enumerate(edges)]
        ind = [0] * nv; outd = [0] * nv
        for (t, h) in fl:
            outd[t] += 1; ind[h] += 1
        typ = {}; ok = True
        for v in range(nv):
            io = (ind[v], outd[v])
            if deg[v] == 2:
                typ[v] = 'c' if io == (2, 0) else 'b'
            elif io == (1, 2):
                typ[v] = 'b'
            elif io == (3, 0):
                typ[v] = 'c'
            else:
                ok = False
        if ok and sum(1 for v in typ if typ[v] == 'c') == 1:
            out.append((fl, typ))
    return out


def check_vertex_formulas():
    """b1 bullet formulas from the literal Grassmann derivative (single vertex, edges e' < e'' < e)."""
    u = [rq(rng) for _ in range(H)]; J = [rq(rng) for _ in range(H)]; lam = [rq(rng) for _ in range(H)]
    vg = [rq(rng) for _ in range(H)]
    res = {}
    # pass-through: tail e'=0 (index i), head e=1 (index j), edge 1 = (x, v) so derivative inserts +y = J
    edges = [(0, 9), (5, 0)]
    y = [[Fr(0)] * H, J]
    inc, T = vertex_tensor(0, 'b', [0], [1], edges, y, {0: dict(u=u, vg=vg, lam=lam)})
    res['pass-through = u^j J_i - delta^j_i (J.u)'] = all(
        T[(i, j)] == u[j] * J[i] - (sum(J[a] * u[a] for a in range(H)) if i == j else 0)
        for i in range(H) for j in range(H))
    # trivalent b_0: tails e'=0 (i'), e''=1 (i''), head e=2 (j)
    edges = [(0, 7), (0, 8), (6, 0)]
    y = [[Fr(0)] * H, [Fr(0)] * H, J]
    inc, T = vertex_tensor(0, 'b', [0, 1], [2], edges, y, {0: dict(u=u, vg=vg, lam=lam)})
    f = lambda i1, i2, j: (J[i1] if j == i2 else 0) - (J[i2] if j == i1 else 0)
    res['trivalent b0: dL/dmu_{i\'} dL/dmu_{i\'\'} = -(delta^j_{i\'\'} J_{i\'} - delta^j_{i\'} J_{i\'\'})'] = all(
        T[(i1, i2, j)] == -f(i1, i2, j) for i1 in range(H) for i2 in range(H) for j in range(H))
    # source: tails e'=0 (i'), e''=1 (i''), external d_a gamma_b
    edges = [(0, 7), (0, 8)]
    y = [[Fr(0)] * H, [Fr(0)] * H]
    inc, T = vertex_tensor(0, 'b', [0, 1], [], edges, y, {0: dict(u=u, vg=vg, lam=lam)})
    res['source: = -(lam_{i\'} ga_{i\'\'} - lam_{i\'\'} ga_{i\'})'] = all(
        T[(i1, i2)] == -(lam[i1] * vg[i2] - lam[i2] * vg[i1]) for i1 in range(H) for i2 in range(H))
    # c sink: heads e=0 (j1), e'=1 (j2) both with derivative +J0/+J1; external letter
    J1 = [rq(rng) for _ in range(H)]
    edges = [(7, 0), (8, 0)]
    y = [J, J1]
    inc, T = vertex_tensor(0, 'c', [], [0, 1], edges, y, {0: dict(u=u, vg=vg, lam=lam)})
    def csink(j1, j2):
        s = Fr(0)
        for (i, j, k, l, m) in permutations(range(H)):
            ep = Fr(eps5((i, j, k, l, m)), 6)
            # maps of (gamma_i, d_j gamma_k, d_l gamma_m) to (e0, e1, ext)
            for pos in permutations(range(3)):
                lets = [(i, None), (k, j), (m, l)]
                (a0, d0), (a1, d1), (a2, d2) = lets[pos[0]], lets[pos[1]], lets[pos[2]]
                if a0 != j1 or a1 != j2:
                    continue
                val = ep * (J[d0] if d0 is not None else 1) * (J1[d1] if d1 is not None else 1)
                val *= vg[a2] * (lam[d2] if d2 is not None else 1)
                s += val
        return s
    res['c sink: (1/6) eps, J on contracted, lam on external (sum over 6 maps)'] = all(
        T[(j1, j2)] == csink(j1, j2) for j1 in range(H) for j2 in range(H))
    return res


if __name__ == '__main__':
    for k, v in check_vertex_formulas().items():
        print('vertex factor', k, ':', v)
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/sugra-11d-qme/code')
    import flows as EXPL                     # explorer's implementation, used only as the object under test
    from htk.graph import theta as htk_theta
    tests = [((2, 5, 6), 'c@interior'), ((2, 5, 6), "c@N"), ((4, 4, 5), "c@N'"), ((4, 4, 5), 'c@interior'),
             ((3, 4, 6), 'c@interior')]
    for pqr, want in tests:
        nv, edges, paths = theta_graph(*pqr)
        G, interior = htk_theta(*pqr)
        assert G.edges == edges and interior == [pv for pv, pe in paths]
        fls = flows_m1(nv, edges)
        assert len(fls) == 40
        cand = []
        for fl, typ in fls:
            cv = [v for v in typ if typ[v] == 'c'][0]
            role = 'c@N' if cv == 0 else ("c@N'" if cv == nv - 1 else 'c@interior')
            if role == want:
                cand.append((fl, typ))
        fl, typ = cand[rng.randrange(len(cand))]
        lam = {v: [rq(rng) for _ in range(H)] for v in range(1, nv - 1)}
        lam[0] = [Fr(0)] * H
        lam_root = [-sum(lam[v][i] for v in range(nv - 1)) for i in range(H)]
        lam_all = dict(lam); lam_all[nv - 1] = lam_root
        ext = {v: dict(u=[rq(rng) for _ in range(H)], vg=[rq(rng) for _ in range(H)], lam=lam_all[v]) for v in range(nv)}
        x = [rsimplex(rng, n_) for n_ in pqr]; w = rsimplex(rng, 3)
        t = [x[P][j] / w[P] for P in range(3) for j in range(pqr[P])]
        y = kirchhoff(nv, edges, t, lam)
        t0 = time.time()
        mine = grouped_integrand(nv, edges, paths, fl, typ, y, ext)
        t1 = time.time()
        pats = EXPL.patterns(G, fl, typ)
        theirs = EXPL.flow_integrand(G, interior, fl, typ, pats, y, ext)
        t2 = time.time()
        print(pqr, want, 'flow', fls.index((fl, typ)), '| grouped (mine) == explorer flow_integrand:', mine == theirs,
              '| nonzero:', mine != 0, '| letter maps', len(pats[0]),
              '| sec mine %.1f, explorer %.1f' % (t1 - t0, t2 - t1))
        if mine != theirs:
            print('   mine', mine, ' theirs', theirs, ' ratio', (mine / theirs) if theirs else None)
