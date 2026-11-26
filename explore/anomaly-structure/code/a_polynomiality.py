"""Part (a): T = 0, K_Gamma(lambda; 0) is a polynomial in the H x H minors of the source matrix,
homogeneous of degree n in each component, det^n under GL(H).

Checks (exact arithmetic):
 (i)   engine K with symbolic sources at all non-root vertices: degree n in each component;
 (ii)  K(g lambda) = det(g)^n K(lambda) for random integer g in GL(H);
 (iii) K lies in the span of the products of n minors [lambda_S] (exact linear solve), small cases;
 (iv)  pointwise identity at random rational t (independent brute force, mybrute.py):
         D_psi(t) = sigma_T^H eps(n,H) D_minor(t),
       D_psi: density of U^{-H} prod_{i,r} psi^i_r (positive cycle matrix),
       D_minor: density of sum_{S_1..S_n} prod_r [lambda_{S_r}] ^_r dR_{c_r,S_r}  (dR exact),
       including non-series-parallel graphs;
 (v)   engine oracle (direct lattice integration, no reduction) against the exact value.
"""
import sys, time, random, itertools
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import sympy as sp
import htk
from htk import catalog, graph as hg
import mybrute as mb

rng = random.Random(20261120)


def eps_nH(n, H):
    return (-1) ** (n * H + (n * (n - 1) // 2) * (H * (H - 1) // 2))


def tree_data(nv, edges):
    """BFS tree of mybrute.cycle_matrix (same construction), sigma_T = eps_T det B_T,
    unflipped fundamental cycles C_T and cotree in increasing order."""
    N = len(edges)
    adj = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        adj[a].append((b, e))
        adj[b].append((a, e))
    root = nv - 1
    seen = {root}
    order = [root]
    tree = []
    k = 0
    while k < len(order):
        u = order[k]
        k += 1
        for (w, e) in sorted(adj[u], key=lambda x: x[1]):
            if w not in seen:
                seen.add(w)
                order.append(w)
                tree.append(e)
    tree = sorted(tree)
    cot = [e for e in range(N) if e not in tree]
    perm = tree + cot
    inv = sum(1 for i in range(N) for j in range(i + 1, N) if perm[i] > perm[j])
    epsT = -1 if inv % 2 else 1
    B = mb.incidence(nv, edges)
    BT = [B[e] for e in tree]
    detB = mb.mdet(BT)
    return tree, cot, epsT * int(detB)


def minor_density(nv, edges, t, lam, H):
    """density (w.r.t. iota_E dt_0..dt_{N-1}) of sum_S prod_r [lam_{S_r}] ^_r dR_{c_r,S_r}."""
    N = len(edges)
    g = nv - 1
    tree, cot, sigma = tree_data(nv, edges)
    R, dR = mb.response_and_derivs(nv, edges, t)
    subsets = list(itertools.combinations(range(g), H))
    lamv = [[Fr(x) for x in lam.get(v, [0] * H)] for v in range(g)]
    br = {S: mb.mdet([[lamv[s][i] for s in S] for i in range(H)]) for S in subsets}
    # one-forms dR_{c,v} = sum_j dR[j][c][v] dt_j ; density of a wedge of N-1 one-forms:
    # [form]_{dt_1..dt_{N-1}} / t_0
    total = Fr(0)
    if len(subsets) ** len(cot) > 20000:
        # same identity without the minor expansion: rows sum_v lam^i_v dR_{c,v}
        rows = []
        for c in cot:
            for i in range(H):
                rows.append([sum(lamv[v][i] * dR[j][c][v] for v in range(g)) for j in range(1, N)])
        return mb.mdet(rows) / Fr(t[0]), sigma, "dy"
    for Ss in itertools.product(subsets, repeat=len(cot)):
        coef = Fr(1)
        for S in Ss:
            coef *= br[S]
        if coef == 0:
            continue
        rows = []
        for c, S in zip(cot, Ss):
            for v in S:
                rows.append([dR[j][c][v] for j in range(1, N)])
        total += coef * mb.mdet(rows) / Fr(t[0])
    return total, sigma, "minors"


def run_graph_numeric(G, H, name, trials=2):
    """engine with exact rational (numeric) sources: K(g lam) = det(g)^n K(lam), and scaling of one
    component by s multiplies K by s^n."""
    n = G.ne - G.nv + 1
    out = {"name": name, "H": H, "n": n, "N": G.ne, "nv": G.nv, "mode": "numeric sources"}
    ok_pt = True
    for k in range(2):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(G.ne)]
        lam = {v: [rng.randint(-7, 7) for _ in range(H)] for v in range(G.nv - 1)}
        q, U = mb.density(G.nv, G.edges, t, lam, 0)
        dm, sigma, how = minor_density(G.nv, G.edges, t, lam, H)
        ok_pt &= (q == sigma ** H * eps_nH(n, H) * dm) and q != 0
    out["pointwise"] = ok_pt
    out["pointwise_via"] = how
    red = htk.sp_reduce(G)
    dens = red.density([[1] * H] * (G.nv - 1), 0) if False else None
    ok_gl, ok_h, vals = True, True, []
    for k in range(trials):
        lam = [[Fr(rng.randint(-6, 6)) for _ in range(H)] for v in range(G.nv - 1)]
        K0 = htk.kernel_exact(G, lam, T=0)
        g = [[rng.randint(-2, 2) for _ in range(H)] for _ in range(H)]
        while mb.mdet(g) == 0:
            g = [[rng.randint(-2, 2) for _ in range(H)] for _ in range(H)]
        lg = [[sum(g[i][j] * l[j] for j in range(H)) for i in range(H)] for l in lam]
        Kg = htk.kernel_exact(G, lg, T=0)
        ok_gl &= (Kg == mb.mdet(g) ** n * K0)
        i0 = rng.randrange(H)
        s_ = Fr(rng.randint(2, 5))
        ls = [[x * (s_ if i == i0 else 1) for i, x in enumerate(l)] for l in lam]
        Ks = htk.kernel_exact(G, ls, T=0)
        ok_h &= (Ks == s_ ** n * K0)
        vals.append(str(K0))
    out["GL_cov"] = ok_gl
    out["homog"] = ok_h
    out["K_samples"] = vals
    return out


def run_graph(G, H, name, do_span=True, do_engine=True):
    n = G.ne - G.nv + 1
    out = {"name": name, "H": H, "n": n, "N": G.ne, "nv": G.nv}
    # (iv) pointwise
    ok_pt = True
    for k in range(2):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(G.ne)]
        lam = {v: [rng.randint(-7, 7) for _ in range(H)] for v in range(G.nv - 1)}
        q, U = mb.density(G.nv, G.edges, t, lam, 0)
        dm, sigma, how = minor_density(G.nv, G.edges, t, lam, H)
        ok_pt &= (q == sigma ** H * eps_nH(n, H) * dm) and q != 0
    out["pointwise"] = ok_pt
    out["pointwise_via"] = how
    if not do_engine:
        return out
    L = {v: list(sp.symbols('x%d_1:%d' % (v, H + 1))) for v in range(G.nv - 1)}
    t0 = time.time()
    K = sp.expand(htk.kernel_exact(G, L, T=0))
    out["t_engine"] = round(time.time() - t0, 2)
    out["nonzero"] = K != 0
    # (i) degree n in each component
    degs = []
    for i in range(H):
        vars_i = [L[v][i] for v in L]
        P = sp.Poly(K, *vars_i)
        degs.append(sorted(set(sum(m) for m in P.monoms())))
    out["component_degrees"] = degs
    out["homog"] = all(d == [n] for d in degs)
    # (ii) GL(H) covariance
    ok_gl = True
    for k in range(2):
        while True:
            g = sp.Matrix(H, H, lambda i, j: rng.randint(-3, 3))
            if g.det() != 0:
                break
        sub = {}
        for v in L:
            vec = g * sp.Matrix(L[v])
            for i in range(H):
                sub[L[v][i]] = vec[i]
        Kg = sp.expand(K.subs(sub, simultaneous=True))
        ok_gl &= sp.expand(Kg - g.det() ** n * K) == 0
    out["GL_cov"] = ok_gl
    # (iii) span of products of minors
    if do_span:
        verts = sorted(L)
        subsets = list(itertools.combinations(verts, H))
        minors = {S: sp.Matrix([[L[s][i] for s in S] for i in range(H)]).det() for S in subsets}
        prods = list(itertools.combinations_with_replacement(subsets, n))
        cs = sp.symbols('c0:%d' % len(prods))
        expr = sp.expand(sum(c * sp.prod([minors[S] for S in pr]) for c, pr in zip(cs, prods)) - K)
        allv = [x for v in verts for x in L[v]]
        eqs = sp.Poly(expr, *allv).coeffs()
        sol = sp.linsolve(eqs, cs)
        out["span_minors"] = sol != sp.S.EmptySet
        out["n_products"] = len(prods)
    out["K"] = K
    return out


def main():
    graphs = []
    # H = 1: bananas (|V| = 2, |I| = n + 1)
    for n in (1, 2, 3, 4):
        graphs.append((hg.banana(n + 1), 1, "banana%d" % (n + 1), True, True))
    # H = 2, T = 0: all tight graphs n = 1..3 (SP and non-SP), n = 4 SP only + prism pointwise
    for n in (1, 2, 3, 4):
        for G in catalog.tight_graphs(2, n):
            sp_ = catalog.is_series_parallel(G)
            graphs.append((G, 2, G.name + ("" if sp_ else "(nonSP)"), n <= 3, sp_))
    # H = 3: n = 1, 2 (symbolic); n = 3: numeric sources (see below)
    for n in (1, 2):
        for G in catalog.tight_graphs(3, n):
            sp_ = catalog.is_series_parallel(G)
            graphs.append((G, 3, G.name + ("" if sp_ else "(nonSP)"), True, sp_))
    # H = 4: n = 1 (symbolic); n = 2: numeric sources (see below)
    for G in catalog.tight_graphs(4, 1):
        graphs.append((G, 4, G.name, True, True))
    rows = []
    if "--resume" in sys.argv:
        graphs = [x for x in graphs if x[1] == 4]
    for G, H, name, span, eng in graphs:
        t0 = time.time()
        r = run_graph(G, H, name, do_span=span, do_engine=eng)
        r["t"] = round(time.time() - t0, 1)
        rows.append(r)
        print({k: v for k, v in r.items() if k != "K"}, flush=True)
    # H = 3, n = 3 and H = 4, n = 2: numeric sources (symbolic sources at 7 vertices exceed the time bound)
    for G in catalog.tight_graphs(4, 2):
        t0 = time.time()
        r = run_graph_numeric(G, 4, G.name)
        r["t"] = round(time.time() - t0, 1)
        print(r, flush=True)
    for G in catalog.tight_graphs(3, 3):
        sp_ = catalog.is_series_parallel(G)
        t0 = time.time()
        if sp_:
            r = run_graph_numeric(G, 3, G.name)
        else:
            r = run_graph(G, 3, G.name + "(nonSP)", do_span=False, do_engine=False)
        r["t"] = round(time.time() - t0, 1)
        print(r, flush=True)
    # (v) oracle checks at numeric sources
    print("oracle checks (direct lattice integration of the brute-force density):")
    for G, H, lam in [(catalog.tight_graphs(2, 2)[0], 2, None), (catalog.tight_graphs(3, 1)[0], 3, None),
                      (catalog.tight_graphs(2, 3)[0], 2, None), (catalog.tight_graphs(3, 2)[0], 3, None)]:
        lamn = [[rng.randint(-5, 5) for _ in range(H)] for v in range(G.nv - 1)]
        Kex = htk.kernel_exact(G, lamn, T=0)
        val, err = htk.kernel_numeric(G, lamn, T=0)
        print("  %s H=%d: exact %s = %.10f, oracle %.10f +- %.1e, rel.dev %.1e" % (
            G.name, H, Kex, float(Kex), val, err, abs(val - float(Kex)) / max(abs(float(Kex)), 1e-300)), flush=True)
    return rows


if __name__ == "__main__":
    main()
