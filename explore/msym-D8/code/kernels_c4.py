"""Kernels of the pentagon and of the theta graphs (1,4,4), (2,3,4), (3,3,3) on C^4 (H = 4, T = 0).

Engine (htk.spred) with explicit plans of tight chain moves between O and N; own brute force
(brute_c4.py, own sections t(xi)) at random rational points; own exact integration; QMC cross-check.
Output: kernels_c4.log, kernels_c4.pkl (densities, displacement maps, domains for ops_c41.py).
"""
import sys, time, random, pickle, itertools, math
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code')
import sympy as sp
import numpy as np
from htk import sp_reduce
from graphs_c4 import pentagon, theta_O, THETAS, AUT
import brute_c4 as bf

H = 4
OUT = None


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if OUT is not None:
        OUT.write(s + '\n')
        OUT.flush()


def lam_symbols(outputs):
    return {v: [sp.Symbol('l%d_%d' % (v, a)) for a in range(H)] for v in outputs}


# ------------------------------------------------------------------ own sections t(xi)

def section_pentagon_r(r):
    """S2-type variables: r_j = (1/T) sum_{k>j} t_k, T = 1."""
    return [1 - r[0], r[0] - r[1], r[1] - r[2], r[2] - r[3], r[3]]


def own_section(name, G, paths, xi):
    """Own t(xi) for the plans used below; xi: dict with the engine's coordinate names."""
    m = G.ne
    t = [None] * m
    if name == 'pent':
        # one move, chains [e0..e3] (O->4) and [e4] (reversed); gauge t_red = 1
        w2 = xi['w1_2']; w1 = 1 - w2
        xs = [1 - xi['x1_1_2'] - xi['x1_1_3'] - xi['x1_1_4'], xi['x1_1_2'], xi['x1_1_3'], xi['x1_1_4']]
        R1, R2 = 1 / w1, 1 / w2
        for j in range(4):
            t[j] = xs[j] * R1
        t[4] = R2
        return t
    (P1, _), (P2, _), (P3, _) = paths
    if len(P1) + len(P2) + len(P3) == 9 and len(P1) == 3 and len(P2) == 3:
        w = [1 - xi['w1_2'] - xi['w1_3'], xi['w1_2'], xi['w1_3']]
        for i, P in enumerate((P1, P2, P3)):
            xs = [1 - xi['x1_%d_2' % (i + 1)] - xi['x1_%d_3' % (i + 1)], xi['x1_%d_2' % (i + 1)], xi['x1_%d_3' % (i + 1)]]
            for j, e in enumerate(P):
                t[e] = xs[j] / w[i]
        return t
    # two moves: (P1, P2) -> a ; (a, P3) -> final
    u2 = xi['w2_2']; u1 = 1 - u2
    ta, R3 = 1 / u1, 1 / u2
    w2 = xi['w1_2']; w1 = 1 - w2
    for i, (P, w) in enumerate(((P1, w1), (P2, w2))):
        L = len(P)
        if L == 1:
            xs = [Fr(1)]
        else:
            rest = [xi['x1_%d_%d' % (i + 1, j)] for j in range(2, L + 1)]
            xs = [1 - sum(rest)] + rest
        for j, e in enumerate(P):
            t[e] = xs[j] * ta / w
    rest = [xi['x2_2_%d' % j] for j in range(2, len(P3) + 1)]
    xs = [1 - sum(rest)] + rest
    for j, e in enumerate(P3):
        t[e] = xs[j] * R3
    return t


def det_t_dxi(tfun, names, pt):
    """det[t; dt/dxi_1; ...] at pt (exact, by sympy differentiation of the section)."""
    syms = [sp.Symbol(n) for n in names]
    texpr = tfun({n: s for n, s in zip(names, syms)})
    rows = [[sp.nsimplify(x) for x in texpr]]
    for s in syms:
        rows.append([sp.diff(x, s) for x in texpr])
    Mx = sp.Matrix(rows)
    val = Mx.subs({s: sp.Rational(pt[n].numerator, pt[n].denominator) for s, n in zip(syms, names)})
    d = val.det()
    return Fr(int(sp.Rational(d).p), int(sp.Rational(d).q)), [Fr(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in val.row(0)]


# ------------------------------------------------------------------ own exact integration

def integrate_poly(expr, blocks, extra_syms=()):
    """int over the product of open simplices; blocks: list of lists of independent coordinate names
    (each block: b_0 = 1 - sum b_i eliminated).  Monomial prod b_i^{a_i} -> prod a_i! / (k + sum a_i)!."""
    coords = [c for blk in blocks for c in blk]
    syms = [sp.Symbol(c) for c in coords]
    P = sp.Poly(sp.expand(expr), *syms)
    pos = {c: i for i, c in enumerate(coords)}
    tot = sp.Integer(0)
    for mon, c in P.terms():
        w = sp.Integer(1)
        for blk in blocks:
            a = [mon[pos[x]] for x in blk]
            w *= sp.Mul(*[sp.factorial(x) for x in a]) / sp.factorial(len(blk) + sum(a))
        tot += c * w
    return sp.expand(tot)


# ------------------------------------------------------------------ brackets

def bracket(lam, S):
    return sp.Matrix([lam[v] for v in S]).det()


def standard_pairs(outputs):
    """Pairs (S, S') of 4-subsets with S_i <= S'_i (standard bitableaux of shape 4 x 2)."""
    subs = list(itertools.combinations(outputs, 4))
    return [(S, T) for S in subs for T in subs if S <= T and all(a <= b for a, b in zip(S, T))]


def exact_lstsq(rows, rhs):
    """Unique exact solution of the overdetermined consistent system rows x = rhs (Fractions); asserts full
    column rank and consistency."""
    n = len(rows[0])
    A = [[Fr(x) for x in r] + [Fr(b)] for r, b in zip(rows, rhs)]
    piv_cols, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv_cols.append(c)
        r += 1
    assert r == n, 'rank deficit %d < %d' % (r, n)
    assert all(A[i][n] == 0 for i in range(r, len(A))), 'inconsistent'
    return [A[i][n] for i in range(n)]


def to_brackets_deg8(poly_fn, outputs, rng):
    """Express an SL(4)-invariant polynomial (callable on numeric lam dicts) of degree 8 as
    sum c_{S,S'} <S><S'> over standard pairs (S_i <= S'_i); exact solve at random integer points."""
    pairs = standard_pairs(outputs)
    rows, rhs = [], []
    for _ in range(len(pairs) + 12):
        lamn = {v: [rng.randint(-6, 6) for _ in range(H)] for v in outputs}
        Mv = {S: int(sp.Matrix([lamn[v] for v in S]).det()) for S in itertools.combinations(outputs, 4)}
        rows.append([Mv[S] * Mv[T] for (S, T) in pairs])
        rhs.append(Fr(int(sp.Rational(poly_fn(lamn)).p), int(sp.Rational(poly_fn(lamn)).q)))
    sol = exact_lstsq(rows, rhs)
    return {pairs[i]: sol[i] for i in range(len(pairs)) if sol[i] != 0}


def to_brackets_deg4(poly_fn, outputs, rng):
    subs = list(itertools.combinations(outputs, 4))
    rows, rhs = [], []
    for _ in range(len(subs) + 8):
        lamn = {v: [rng.randint(-6, 6) for _ in range(H)] for v in outputs}
        rows.append([sp.Matrix([lamn[v] for v in S]).det() for S in subs])
        rhs.append(poly_fn(lamn))
    A = sp.Matrix(rows); b = sp.Matrix(rhs)
    sol = A.solve_least_squares(b)
    assert A * sol == b
    return {subs[i]: sol[i] for i in range(len(subs)) if sol[i] != 0}


# ------------------------------------------------------------------ main per graph

def analyse(name, G, spokes, outputs, paths, plan, nlam_points=4, qmc_n=2 ** 17):
    t0 = time.time()
    log('=' * 100)
    log('graph', name, ' nv =', G.nv, ' edges =', G.edges, ' loops =', G.loops, ' root(O) =', G.root)
    log('spokes (edge: (input j, sign of z_j)):', spokes, ' outputs:', outputs)
    red = sp_reduce(G, plan=plan)
    log(red.describe())
    blocks = [list(blk.indep) for blk in red.blocks]
    log('blocks (independent coordinates; first barycentric ratio of each block eliminated):', [list(b.bary) for b in red.blocks])
    log('t_e = prod ratio^exp:', {e: red.texp[e] for e in range(G.ne)})
    lam = lam_symbols(outputs)
    symbolic = G.loops == 1 or name != (3, 3, 3)
    if symbolic:
        dens = red.density(lam, T=0, check=3)
        log('engine density method:', dens.method, ' selfcheck(3 points) passed')
    else:
        dens = None
        log('engine density: evaluated with numeric sources at each point below (symbolic 8x8 local determinant not expanded)')
    big = (dens is None) or any(len(f.terms()) > 3000 for f in dens.factors)
    gexpr = None if big else sp.Mul(*[sp.factor(f.as_expr()) for f in dens.factors])
    log('density g(xi; lam) =', gexpr if not big else '(engine local density not expanded here)')
    lamflat = [s_ for v in outputs for s_ in lam[v]]
    for f in ([] if big else dens.factors):
        fe = sp.factor(f.as_expr())
        for ff in (sp.Mul.make_args(fe)):
            base = ff.as_base_exp()[0]
            Pf = sp.Poly(sp.expand(base), *lamflat)
            if Pf.total_degree() == 4:
                cS = {}
                for S in itertools.combinations(outputs, 4):
                    mon = [0] * len(lamflat)
                    for i, v in enumerate(S):
                        mon[lamflat.index(lam[v][i])] = 1
                    c = Pf.coeff_monomial(tuple(mon))
                    if c != 0:
                        cS[S] = sp.factor(c)
                rebuilt = sum(c * bracket(lam, S) for S, c in cS.items())
                assert sp.expand(rebuilt - base) == 0
                log('   factor (degree 4 in lam) =', ' + '.join('(%s)<%s>' % (c, ''.join(str(v) for v in S)) for S, c in cS.items()),
                    '' if ff.as_base_exp()[1] == 1 else '  [power %s]' % ff.as_base_exp()[1])
            elif Pf.total_degree() > 0:
                log('   factor of degree', Pf.total_degree(), 'in lam,', len(Pf.terms()), 'terms')
            else:
                log('   factor (no lam):', ff)
    # displacement maps: X_v = sum_e Xmat[v][e] z_e ; z_e = sign * z_j on spokes
    Xm = red.displacement_expr()
    X = {}
    for v in outputs:
        X[v] = {}
        for e, (j, s) in spokes.items():
            X[v][j] = sp.factor(s * Xm[v][e])
    log('displacement maps X_v = sum_j c_{v,j} z_j:')
    for v in outputs:
        log('   X_%d =' % v, ' + '.join('(%s) z_%d' % (X[v][j], j) for j in sorted(X[v])))
        assert sp.simplify(sum(X[v].values()) - 1) == 0
    log('   sum_j c_{v,j} = 1 for every v: True')

    # ---------------- own brute force, pointwise
    rng = random.Random(7)
    names = list(red.coords)
    cot = bf.cotree_of(G.nv, G.edges)
    log('own cotree (Kruskal, edge order):', cot)
    zsym = {j: [sp.Symbol('z%d_%d' % (j, a)) for a in range(H)] for j in sorted({j for j, s in spokes.values()})}
    for it in range(nlam_points):
        pt = red.random_point(rng)
        xi = {n: pt[n] for n in names}
        lamn = {v: [rng.randint(-9, 9) for _ in range(H)] for v in outputs}
        tf = lambda X_: own_section(name, G, paths, X_)
        J, tvals = det_t_dxi(tf, names, xi)
        t_eng = red.tmap(pt)
        assert [Fr(x) for x in t_eng] == tvals, ('sections differ', t_eng, tvals)
        Dt, y = bf.density_D(G.nv, G.edges, tvals, lamn, cot)
        bf.check_kirchhoff(G.nv, G.edges, tvals, lamn, y)
        g_own = Dt * abs(J)
        subs = {sp.Symbol(c): sp.Rational(xi[c].numerator, xi[c].denominator) for c in names}
        subs.update({s: lamn[v][a] for v in outputs for a, s in enumerate(lam[v])})
        if dens is not None:
            psyms = sorted([s_ for v in outputs for s_ in lam[v]], key=str)
            g_eng, _rad = dens.evaluate(red._vals(pt) + tuple(Fr(int(subs[s_])) for s_ in psyms))
        else:
            dn = red.density([lamn.get(v, [0] * H) for v in range(G.nv - 1)], T=0, check=1)
            g_eng, _rad = dn.evaluate(red._vals(pt))
        g_eng = sp.Rational(g_eng.numerator, g_eng.denominator)
        assert sp.Rational(g_eng) == sp.Rational(g_own.numerator, g_own.denominator), (g_eng, g_own)
        # displacement: - sum_e y_e . z_e == sum_v lam_v . X_v at random rational z
        zn = {j: [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(H)] for j in zsym}
        lhs = Fr(0)
        for e, (j, s) in spokes.items():
            lhs -= sum(y[e][a] * s * zn[j][a] for a in range(H))
        rhs = sp.Integer(0)
        for v in outputs:
            for j, c in X[v].items():
                cv = c.subs(subs)
                rhs += cv * sum(lamn[v][a] * sp.Rational(zn[j][a].numerator, zn[j][a].denominator) for a in range(H))
        assert sp.Rational(rhs) == sp.Rational(lhs.numerator, lhs.denominator)
        log('  point %d: own section == engine tmap; own Kirchhoff density D(t)|det[t;dt/dxi]| = %s == engine g; '
            'sign det[t;dt/dxi] = %+d; -sum y.z == sum lam.X: True' % (it, g_own, 1 if J > 0 else -1))

    # ---------------- zero-displacement value
    t1 = time.time()
    lamflat = [s_ for v in outputs for s_ in lam[v]]
    weights = {}
    if G.loops == 1:
        K0_eng = red.integrate(lam, T=0, check=0, dens=dens)
        PK = sp.Poly(K0_eng, *lamflat)
        fn = lambda lamn: sp.Rational(PK.eval(tuple(sp.Integer(lamn[v][a]) for v in outputs for a in range(H))))
        K0_store = K0_eng
        log('K(lam;0) (engine, symbolic) = %s (%.1fs)' % (sp.factor(K0_eng), time.time() - t1))
    else:
        from weights_c4 import ThetaWeights
        W = ThetaWeights(name, log=log)
        zero = ((0,) * H,) * 3
        weights[zero] = W.weights(zero)
        for j in range(3):
            for a in range(H):
                ks = tuple(tuple(1 if (jj == j and aa == a) else 0 for aa in range(H)) for jj in range(3))
                weights[ks] = W.weights(ks)
        w0 = weights[zero]
        idx = {v: i for i, v in enumerate(outputs)}

        def fn(lamn, w0=w0):
            tot = Fr(0)
            for m, c in w0.items():
                term = c
                for v in outputs:
                    for aa in range(H):
                        if m[idx[v]][aa]:
                            term *= Fr(lamn[v][aa]) ** m[idx[v]][aa]
                tot += term
            return sp.Rational(tot.numerator, tot.denominator)
        K0_store = None
        log('K(lam;0) by Cauchy-Binet minors (own): %d monomials in lam; first-order z-weights: %s monomials (%.1fs)'
            % (len(w0), [len(weights[k]) for k in weights if k != zero], time.time() - t1))
        rk = random.Random(21)
        for _ in range(3):
            lamn = {v: [rk.randint(-7, 7) for _ in range(H)] for v in outputs}
            subsl = {s_: lamn[v][aa] for v in outputs for aa, s_ in enumerate(lam[v])}
            own = integrate_poly(gexpr.subs(subsl), blocks) if gexpr is not None else None
            eng = red.integrate([lamn.get(v, [0] * H) for v in range(G.nv - 1)], T=0, check=0)
            cb = fn(lamn)
            assert cb == sp.Rational(eng.numerator, eng.denominator), (cb, eng)
            assert own is None or sp.Rational(own) == cb
            log('  lam = %s: engine K(lam;0) (numeric sources) = %sCauchy-Binet value = %s' % (lamn, '' if own is None else 'own integration of the engine density = ', cb))
        # first-order weights against the engine Taylor coefficient at random z
        lamn = {v: [rk.randint(-7, 7) for _ in range(H)] for v in outputs}
        zn = {j: [rk.randint(-5, 5) for _ in range(H)] for j in (1, 2, 3)}
        zl = [[0] * H for _ in range(G.ne)]
        for e, (j, sgn) in spokes.items():
            zl[e] = [sgn * x for x in zn[j]]
        K1_eng = red.taylor([lamn.get(v, [0] * H) for v in range(G.nv - 1)], zl, 1, T=0, check=0)[1]
        K1_cb = Fr(0)
        for j in range(3):
            for aa in range(H):
                ks = tuple(tuple(1 if (jj == j and ab == aa) else 0 for ab in range(H)) for jj in range(3))
                tot = Fr(0)
                for m, c in weights[ks].items():
                    term = c
                    for v in outputs:
                        for ab in range(H):
                            if m[idx[v]][ab]:
                                term *= Fr(lamn[v][ab]) ** m[idx[v]][ab]
                    tot += term
                K1_cb += Fr(zn[j + 1][aa]) * tot
        assert Fr(K1_eng) == K1_cb, (K1_eng, K1_cb)
        log('  first order in z at lam = %s, z = %s: engine taylor K_1 = Cauchy-Binet first-order weights = %s' % (lamn, zn, K1_cb))
    if G.loops == 1:
        br = to_brackets_deg4(fn, outputs, random.Random(3))
    else:
        br = to_brackets_deg8(fn, outputs, random.Random(3))
    log('K(lam;0) in brackets <ijkl> = det(lam_i,lam_j,lam_k,lam_l) (vertex labels):')
    log('   ', br)

    # ---------------- QMC of the brute-force density in the own section (floating point)
    from scipy.stats import qmc
    rq = random.Random(11)
    lamn = {v: [rq.randint(-5, 5) for _ in range(H)] for v in outputs}
    exact = fn(lamn)
    dim = len(names)
    sob = qmc.Sobol(d=dim, scramble=True, seed=5)
    U01 = sob.random(qmc_n) * (1 - 2e-9) + 1e-9
    # map the unit cube to the product of simplices (sorted spacings) with Jacobian 1/prod k!
    xi_arr = {}
    jac = 1.0
    col = 0
    for blk in blocks:
        k = len(blk)
        u = np.sort(U01[:, col:col + k], axis=1)
        sp_ = np.diff(np.concatenate([np.zeros((qmc_n, 1)), u, np.ones((qmc_n, 1))], axis=1))  # k+1 spacings
        for i, c in enumerate(blk):
            xi_arr[c] = sp_[:, i + 1]
        jac /= math.factorial(k)
        col += k
    tsym = own_section(name, G, paths, {c: sp.Symbol(c) for c in names})
    csyms = [sp.Symbol(c) for c in names]
    tf = sp.lambdify(csyms, tsym, 'numpy')
    Jexpr = sp.factor(sp.Matrix([tsym] + [[sp.diff(x, s) for x in tsym] for s in csyms]).det())
    Jf = sp.lambdify(csyms, Jexpr, 'numpy')
    args = [xi_arr[c] for c in names]
    T = np.array([np.broadcast_to(np.asarray(x, dtype=float), (qmc_n,)) for x in tf(*args)]).T
    Jv = np.abs(np.broadcast_to(np.asarray(Jf(*args), dtype=float), (qmc_n,)))
    D = bf.density_float_batch(G.nv, G.edges, T, lamn, cot)
    vals = D * Jv
    est = vals.mean() / math.factorial(1) * (1.0 / (1.0 / jac)) if False else vals.mean() * jac
    rel = abs(est - float(exact)) / abs(float(exact))
    log('QMC (Sobol, %d points, brute-force float Kirchhoff density in own section), lam = %s:' % (qmc_n, lamn))
    log('   QMC = %.12g, exact = %s = %.12g, relative difference %.2e (%d agreeing digits)'
        % (est, exact, float(exact), rel, int(-math.log10(rel)) if 0 < rel < 1 else (16 if rel == 0 else 0)))
    log('time %.1fs' % (time.time() - t0))
    return dict(name=name, G=(G.nv, G.edges), spokes=spokes, outputs=outputs, coords=names, blocks=blocks,
                density=gexpr, X=X, K0=K0_store, weights=weights, K0_brackets=br, lam=lam, orientation=red.orientation)


def main(which):
    tag = '_'.join(w.strip('()').replace(',', '') for w in which)
    pk = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/kernels_c4_%s.pkl' % tag
    res = {}
    if 'pent' in which:
        main_pent(res)
        with open(pk, 'wb') as f:
            pickle.dump(res, f)
        log('saved', pk, 'pent')
    for pqr in THETAS:
        if str(pqr).replace(' ', '') in which:
            G, spokes, outputs, paths = theta_O(*pqr)
            O, N = G.root, 0
            P = [p_[0] for p_ in paths]
            if pqr == (3, 3, 3):
                plan = [(O, N, [P[0], P[1], P[2]])]
            else:
                plan = [(O, N, [P[0], P[1]]), (O, N, [[G.ne], P[2]])]
            res[pqr] = analyse(pqr, G, spokes, outputs, paths, plan)
            res[pqr]['paths'] = paths
            with open(pk, 'wb') as f:
                pickle.dump(res, f)
            log('saved', pk, pqr)


def main_pent(res):
    G, spokes, outputs = pentagon()
    plan = [(4, 3, [[0, 1, 2, 3], [4]])]
    res['pent'] = analyse('pent', G, spokes, outputs, None, plan)
    # pentagon in the S2 variables r_j
    rng = random.Random(1)
    lamn = {v: [rng.randint(-9, 9) for _ in range(H)] for v in outputs}
    rs = sp.symbols('r1:5')
    tex = section_pentagon_r(rs)
    Jr = sp.Matrix([tex] + [[sp.diff(x, s) for x in tex] for s in rs]).det()
    detl = sp.Matrix([lamn[v] for v in outputs]).det()
    for _ in range(4):
        rv = sorted([Fr(rng.randint(1, 999), 1000) for _ in range(4)], reverse=True)
        tv = [Fr(int(sp.Rational(x.subs(dict(zip(rs, rv)))).p), int(sp.Rational(x.subs(dict(zip(rs, rv)))).q)) for x in tex]
        Dt, y = bf.density_D(G.nv, G.edges, tv, lamn, [0])
        Jv = sp.Rational(Jr.subs(dict(zip(rs, rv))))
        y1 = [-sum(r * lamn[j][a] for r, j in zip(rv, outputs)) for a in range(H)]
        assert y[0] == y1
        log('pentagon, r-variables: y_{e0} = -sum_j r_j lam_j: True; det[t;dt/dr] = %s; D(t) det[t;dt/dr] = %s; det(lam_1..lam_4) = %s'
            % (Jv, Dt * Fr(int(Jv.p), int(Jv.q)), detl))


if __name__ == '__main__':
    which = sys.argv[1:] or ['pent', '(1,4,4)', '(2,3,4)', '(3,3,3)']
    OUT = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/code/kernels_c4_%s.log' % '_'.join(w.strip('()').replace(',', '') for w in which), 'w')
    main(which)
