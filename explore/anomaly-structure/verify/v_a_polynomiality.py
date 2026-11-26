"""VERIFY (a1): T = 0 polynomiality, independent of htk and of the explorer's code.

A. Pointwise (exact, own Kirchhoff solver vk.py, random spanning trees, forward-mode derivatives):
   A1  D_psi(t) does not depend on the positive cycle basis (two random trees);
   A2  D_psi(t) = sigma_T^H eps(n,H) D_dy(t), D_dy the density of wedge_r d^H y_{c_r} for a third
       random tree T, derivatives by exact forward mode (not by an R / dR formula);
   A3  D_psi(t; lam with component i scaled by s) = s^n D_psi(t; lam);
   A4  D_psi(t; g lam) = det(g)^n D_psi(t; lam);
   A5  (small cases) D_dy = sum_{S_1..S_n} prod_r [lam_{S_r}] * density(wedge_r dR_{c_r,S_r}),
       dR_{c,v} from forward mode with a unit source at v.
B. Integrated, exact: Grundmann-Moeller rules (exact rational nodes) in own series-parallel
   coordinates (sp_coords.py) for k parallel paths (bananas, thetas), and on the t-simplex for
   polygons; stabilisation in the GM index s; K(g lam) = det(g)^n K(lam), K(scaled) = s^n K;
   comparison with htk.kernel_exact (engine) as a cross-check only.
"""
import sys, time, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/verify')
import vk, sp_coords as spc

rng = random.Random(777)


def rand_t(N):
    return [Fr(rng.randint(1, 30), rng.randint(1, 30)) for _ in range(N)]


def rand_lam(nv, H, lo=-6, hi=6):
    return {v: [rng.randint(lo, hi) for _ in range(H)] for v in range(nv - 1)}


def rand_g(H):
    while True:
        g = [[rng.randint(-2, 2) for _ in range(H)] for _ in range(H)]
        if vk.det(g) != 0:
            return g


def act(g, lam):
    H = len(g)
    return {v: [sum(g[i][j] * l[j] for j in range(H)) for i in range(H)] for v, l in lam.items()}


def minor_expansion(nv, edges, t, lam, H, tree):
    N = len(edges)
    g = nv - 1
    cot = [e for e in range(N) if e not in tree]
    # dR[v][j][e] = d y_e / d t_j for unit source at v (H = 1)
    dR = {v: vk.current_derivs(nv, edges, t, {v: [1]}, 1) for v in range(g)}
    subsets = list(itertools.combinations(range(g), H))
    br = {S: vk.det([[lam[s][i] for s in S] for i in range(H)]) for S in subsets}
    tot = Fr(0)
    for Ss in itertools.product(subsets, repeat=len(cot)):
        coef = Fr(1)
        for S in Ss:
            coef *= br[S]
        if coef == 0:
            continue
        rows = []
        for c, S in zip(cot, Ss):
            for v in S:
                rows.append([dR[v][j][c][0] for j in range(1, N)])
        tot += coef * vk.det(rows)
    return tot / Fr(t[0])


def pointwise(name, nv, edges, H, do_minor):
    N = len(edges)
    n = N - nv + 1
    res = dict(A1=True, A2=True, A3=True, A4=True, A5=None, nonzero=False)
    for trial in range(3):
        t = rand_t(N)
        lam = rand_lam(nv, H)
        C1 = vk.positive_cycles(nv, edges, vk.random_tree(nv, edges, rng))
        C2 = vk.positive_cycles(nv, edges, vk.random_tree(nv, edges, rng))
        assert vk.check_cycles(nv, edges, C1) and vk.check_cycles(nv, edges, C2)
        q1, U = vk.density(nv, edges, t, lam, H, 0, C1, e0=rng.randrange(N))
        q2, _ = vk.density(nv, edges, t, lam, H, 0, C2, e0=rng.randrange(N))
        res['A1'] &= (q1 == q2)
        res['nonzero'] |= (q1 != 0)  # at least one of the trials nonzero (random sources can be degenerate)
        T3 = vk.random_tree(nv, edges, rng)
        sig = vk.sigma_tree(nv, edges, T3)
        Ddy = vk.dy_density(nv, edges, t, lam, H, T3)
        res['A2'] &= (q1 == sig ** H * vk.eps_nH(n, H) * Ddy)
        i0 = rng.randrange(H)
        s = rng.randint(2, 5)
        lam_s = {v: [x * (s if i == i0 else 1) for i, x in enumerate(l)] for v, l in lam.items()}
        qs, _ = vk.density(nv, edges, t, lam_s, H, 0, C1)
        res['A3'] &= (qs == s ** n * q1)
        g = rand_g(H)
        qg, _ = vk.density(nv, edges, t, act(g, lam), H, 0, C1)
        res['A4'] &= (qg == vk.det(g) ** n * q1)
        if do_minor:
            mex = minor_expansion(nv, edges, t, lam, H, T3)
            res['A5'] = (res['A5'] is not False) and (mex == Ddy)
    return res


def polygon_integral(m, lam, H, s):
    nv, edges = m, [(k, (k + 1) % m) for k in range(m)]
    C = vk.positive_cycles(nv, edges, vk.dfs_tree(nv, edges))
    tot = Fr(0)
    for w, pt in vk.gm_rule(m - 1, s):
        q, U = vk.density(nv, edges, list(pt), lam, H, 0, C)
        tot += w * q
    return tot


if __name__ == "__main__":
    T0 = time.time()
    sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
    from htk import catalog  # test-graph lists only
    import htk
    from htk import graph as hg
    print("=== A. pointwise identities (exact) ===")
    cases = []
    for k in (2, 3, 4, 5):
        cases.append(("banana%d" % k, 2, [(0, 1)] * k, 1))
    for (H, n) in [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2)]:
        gs = catalog.tight_graphs(H, n)
        for idx, G in enumerate(gs):
            cases.append(("D%d_n%d_%d" % (H, n, idx), G.nv, list(G.edges), H))
    allok = True
    for name, nv, edges, H in cases:
        t0 = time.time()
        N = len(edges)
        n = N - nv + 1
        g_ = nv - 1
        nprod = len(list(itertools.combinations(range(g_), H))) ** n
        r = pointwise(name, nv, edges, H, do_minor=(nprod <= 3000))
        ok = r['A1'] and r['A2'] and r['A3'] and r['A4'] and r['nonzero'] and (r['A5'] in (True, None))
        allok &= ok
        print("%-12s H=%d n=%d N=%2d  A1 %s A2 %s A3 %s A4 %s A5 %s nonzero %s  (%.1fs)" % (
            name, H, n, N, r['A1'], r['A2'], r['A3'], r['A4'], r['A5'], r['nonzero'], time.time() - t0), flush=True)
    print("all pointwise checks:", allok)

    print("\n=== B. exact integrals (GM rules in own SP coordinates) ===")
    # polygons
    for m in (3, 4, 5, 6):
        H = m - 1
        lam = rand_lam(m, H)
        vals = [polygon_integral(m, lam, H, s) for s in (1, 2)]
        g = rand_g(H)
        Kg = polygon_integral(m, act(g, lam), H, 2)
        Ge = hg.Graph(m, [(k, (k + 1) % m) for k in range(m)])
        Keng = htk.kernel_exact(Ge, [lam[v] for v in range(m - 1)], T=0)
        print("polygon%d H=%d: K(s=1) = %s, K(s=2) = %s, K(g lam) = det(g)^1 K: %s, engine %s, equal: %s" % (
            m, H, vals[0], vals[1], Kg == vk.det(g) * vals[1], Keng, Keng == vals[1]), flush=True)
    sp_cases = [((1, 1), 1), ((1, 1, 1), 1), ((1, 1, 1, 1), 1), ((1, 1, 1, 1, 1), 1),
                ((1, 2, 2), 2), ((1, 3, 3), 3), ((2, 2, 3), 3), ((1, 4, 4), 4), ((2, 3, 4), 4), ((3, 3, 3), 4)]
    for ns, H in sp_cases:
        t0 = time.time()
        V, edges, interior = spc.paths_graph(ns)
        n = len(edges) - V + 1
        lam = rand_lam(V, H)
        smax = 2 if H <= 3 else 1
        vals = [spc.integrate_T0(ns, lam, H, s) for s in range(1, smax + 1)]
        if H == 4:
            vals.append(spc.integrate_T0(ns, lam, H, 2))
        K = vals[-1]
        g = rand_g(H)
        Kg = spc.integrate_T0(ns, act(g, lam), H, 1)
        i0 = rng.randrange(H)
        sc = rng.randint(2, 4)
        lam_s = {v: [x * (sc if i == i0 else 1) for i, x in enumerate(l)] for v, l in lam.items()}
        Ks = spc.integrate_T0(ns, lam_s, H, 1)
        Keng = htk.kernel_exact(hg.Graph(V, edges), [lam[v] for v in range(V - 1)], T=0)
        print("paths%s H=%d n=%d: K(s=1..) = %s, stable %s; K(g lam) = det(g)^n K: %s; K(scaled comp %d by %d) = s^n K: %s; engine %s, equal: %s (%.1fs)" % (
            ns, H, n, [str(v) for v in vals], len(set(vals)) == 1, Kg == vk.det(g) ** n * K, i0, sc,
            Ks == sc ** n * K, Keng, Keng == K, time.time() - t0), flush=True)
    print("total %.1fs" % (time.time() - T0))
