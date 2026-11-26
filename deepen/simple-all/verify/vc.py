"""Part (c): G-invariant H^3_(a,b), 0 <= a,b <= 2, for simple g with (g^{*(x)3})^g = C f, and d E_f.

Universal (sympy, exact):  cochains of arity k: P(p_0..p_{k-1}) T(X_0..), sum p = 0 (p_{k-1} eliminated).
  C^3_(a,b) = { Sym(m) f },  Sym(m)(q0,q1,q2) = sum_{s in S3} m(q_s0, q_s1)   (alt(m f) = Sym(m) f, f alternating)
  C^2_(a,b) = { (m(q0) - m(q1)) kappa }   (alt(m(q0) kappa), q1 = -q0)
  d c(x_0..x_k) = sum_{i<j} (-1)^{i+j} c([x_i,x_j], x_rest),   [x_i, x_j] carries p_i + p_j
  d(P f)(x0..x3) = sum_{i<j} (-1)^{i+j} P(p_i+p_j, p_k, p_l) ff(i,j,k,l),  ff(i,j,k,l) = kappa([X_i,X_j],[X_k,X_l])
     reduced with ff(0,1,2,3) = T1, ff(0,3,1,2) = T3, ff(0,2,1,3) = T1 + T3 (Jacobi), pair symmetries
     -> c1 T1 + c3 T3;  T1, T3 linearly independent for simple g (checked numerically below) => rank universal.
  d((m(q0)-m(q1)) kappa)(x0,x1,x2) = sum_{i<j} (-1)^{i+j} (m(p_i+p_j) - m(p_k)) kappa([X_i,X_j],X_k) -> c f.
Numeric (exact mod p, explicit structure constants): the same cochains evaluated at random points, ranks of
d2, d3 from evaluation matrices; d E_f at random points; T1, T3 independence; the formula
  d E_f = 2 [ B(p0-p3, p2-p1) T1 + B(p0-p1, p2-p3) T3 ],  B(p,q) = p^1 q^2 + p^2 q^1.
"""
import sys, itertools
import numpy as np
import sympy as sp
import vlie, vpoint

P = [sp.symbols('p%d_1 p%d_2' % (i, i)) for i in range(5)]


def sub_cons(expr, k):
    """eliminate p_{k-1} = -(p_0+..+p_{k-2})."""
    s = {P[k - 1][c]: -sum(P[i][c] for i in range(k - 1)) for c in range(2)}
    return sp.expand(expr.subs(s, simultaneous=True))


def monos(a, b, nv):
    """monomials of bidegree (a,b) in nv vectors."""
    out = []
    for ea in itertools.product(range(a + 1), repeat=nv):
        if sum(ea) != a:
            continue
        for eb in itertools.product(range(b + 1), repeat=nv):
            if sum(eb) != b:
                continue
            out.append((ea, eb))
    return out


def mono_expr(m, vecs):
    ea, eb = m
    e = sp.Integer(1)
    for v, x, y in zip(vecs, ea, eb):
        e *= v[0] ** x * v[1] ** y
    return e


def poly_vec(expr, gens):
    if expr == 0:
        return {}
    return sp.Poly(expr, *gens).as_dict()


def span_basis(exprs, gens):
    """basis (list of exprs) of the span, by exact rank."""
    keys, rows = {}, []
    for e in exprs:
        dct = poly_vec(e, gens)
        rows.append(dct)
        for k in dct:
            keys.setdefault(k, len(keys))
    M = sp.zeros(len(exprs), max(len(keys), 1))
    for i, dct in enumerate(rows):
        for k, v in dct.items():
            M[i, keys[k]] = v
    basis, Mb = [], sp.zeros(0, M.shape[1])
    for i in range(len(exprs)):
        cand = Mb.col_join(M[i, :])
        if cand.rank() > Mb.rank():
            Mb = cand; basis.append(exprs[i])
    return basis


def C3_basis(a, b):
    gens3 = [P[0][0], P[0][1], P[1][0], P[1][1]]
    exprs = []
    for m in monos(a, b, 2):
        e = sum(mono_expr(m, [P[s[0]], P[s[1]]]) for s in itertools.permutations(range(3)))
        e = sub_cons(e, 3)
        if e != 0:
            exprs.append(e)
    return span_basis(exprs, gens3) if exprs else []


def C2_basis(a, b):
    m = ((a,), (b,))
    e = sub_cons(mono_expr(m, [P[0]]) - mono_expr(m, [P[1]]), 2)
    return [e] if e != 0 else []


def subs_args(expr, qs, k_old):
    """expr in P[0..k_old-2] (P[k_old-1] eliminated) evaluated at momenta qs (list of sympy 2-vectors)."""
    s = {P[i][c]: qs[i][c] for i in range(k_old - 1) for c in range(2)}
    return expr.subs(s, simultaneous=True)


# ff(i,j,k,l) in the basis (T1, T3):  canonical pair forms
def ff_coeff(i, j, k, l):
    sg = 1
    if i > j:
        i, j = j, i; sg = -sg
    if k > l:
        k, l = l, k; sg = -sg
    if i > k:
        i, j, k, l = k, l, i, j
    key = (i, j, k, l)
    tab = {(0, 1, 2, 3): (1, 0), (0, 3, 1, 2): (0, 1), (0, 2, 1, 3): (1, 1)}
    c1, c3 = tab[key]
    return sg * c1, sg * c3


def d3_poly(Pf):
    """d(Pf) -> (c1, c3), polynomials in p0,p1,p2 (p3 eliminated)."""
    c1 = c3 = sp.Integer(0)
    vec = lambda i: (P[i][0], P[i][1])
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        qs = [(P[i][0] + P[j][0], P[i][1] + P[j][1]), vec(k), vec(l)]
        val = subs_args(Pf, qs, 3)
        a1, a3 = ff_coeff(i, j, k, l)
        c1 += (-1) ** (i + j) * a1 * val
        c3 += (-1) ** (i + j) * a3 * val
    return sub_cons(c1, 4), sub_cons(c3, 4)


def d2_poly(Pk):
    """d((m(q0)-m(q1)) kappa) -> c (coefficient of f(X0,X1,X2))."""
    c = sp.Integer(0)
    for i, j in itertools.combinations(range(3), 2):
        k = [r for r in range(3) if r not in (i, j)][0]
        qs = [(P[i][0] + P[j][0], P[i][1] + P[j][1])]
        val = subs_args(Pk, qs, 2)
        # kappa([X_i,X_j], X_k) = f(X_i, X_j, X_k) = sgn(i,j,k) f(X0,X1,X2)
        c += (-1) ** (i + j) * vlie.perm_sign((i, j, k)) * val
    return sub_cons(c, 3)


def universal():
    gens4 = [P[i][c] for i in range(3) for c in range(2)]
    gens3 = [P[i][c] for i in range(2) for c in range(2)]
    rows = []
    res = {}
    for a in range(3):
        for b in range(3):
            B3 = C3_basis(a, b); B2 = C2_basis(a, b)
            # rank d3
            vecs = []
            for e in B3:
                c1, c3 = d3_poly(e)
                vecs.append((c1, c3))
            keys = {}
            for c1, c3 in vecs:
                for t, cc in enumerate((c1, c3)):
                    for k in poly_vec(cc, gens4):
                        keys.setdefault((t, k), len(keys))
            M = sp.zeros(len(vecs), max(1, len(keys)))
            for r_, (c1, c3) in enumerate(vecs):
                for t, cc in enumerate((c1, c3)):
                    for k, v in poly_vec(cc, gens4).items():
                        M[r_, keys[(t, k)]] = v
            rk3 = M.rank() if len(vecs) else 0
            # rank d2
            d2 = [d2_poly(e) for e in B2]
            rk2 = sum(1 for c in d2 if c != 0)
            Z3 = len(B3) - rk3
            H3 = Z3 - rk2
            res[(a, b)] = (len(B2), len(B3), rk2, rk3, Z3, H3)
            rows.append('   (%d,%d)  dim C^2 %d  dim C^3 %d  rank d2 %d  rank d3 %d  dim Z^3 %d  dim H^3 %d   C^3 basis: %s'
                        % (a, b, len(B2), len(B3), rk2, rk3, Z3, H3, [sp.factor(e) for e in B3]))
    return res, rows


def dEf_universal():
    Ef = sp.expand(2 * sum(P[i][0] * P[i][1] for i in range(3)))
    Ef = sub_cons(Ef, 3)
    c1, c3 = d3_poly(Ef)
    Bf = lambda u, v: u[0] * v[1] + u[1] * v[0]
    v = lambda i: (P[i][0], P[i][1])
    sub = lambda x, y: (x[0] - y[0], x[1] - y[1])
    t1 = sub_cons(2 * Bf(sub(v(0), v(3)), sub(v(2), v(1))), 4)
    t3 = sub_cons(2 * Bf(sub(v(0), v(1)), sub(v(2), v(3))), 4)
    # alternative sum formula: 2 sum_{i<j} (-1)^{i+j} (p_i^1 p_j^2 + p_i^2 p_j^1) ff(i,j,k,l)
    s1 = s3 = 0
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        a1, a3 = ff_coeff(i, j, k, l)
        val = 2 * Bf(v(i), v(j))
        s1 += (-1) ** (i + j) * a1 * val; s3 += (-1) ** (i + j) * a3 * val
    s1, s3 = sub_cons(s1, 4), sub_cons(s3, 4)
    alt_ok = sub_cons(sum(P[s[0]][0] * P[s[0]][1] for s in itertools.permutations(range(3))) - Ef, 3) == 0
    return (sp.expand(c1 - t1) == 0 and sp.expand(c3 - t3) == 0, sp.expand(c1 - s1) == 0 and sp.expand(c3 - s3) == 0,
            sp.factor(c1), sp.factor(c3), alt_ok)


# ------------------------------------------------------------------ numeric checks on explicit algebras
def numeric(key, p=1048573, npts=6, seed=21):
    L = vlie.build(key)
    E = vpoint.Engine(L, p)
    rng = np.random.default_rng(seed)
    d = L.d

    def rv():
        return rng.integers(0, p, d).astype(np.float64)

    def fval(X, Y, Z):
        return E.kap(E.bracket(X, Y), Z)

    # T1, T3 independence: values at random points
    pts = [[rv() for _ in range(4)] for _ in range(npts)]
    T = [(E.kap(E.bracket(x[0], x[1]), E.bracket(x[2], x[3])), E.kap(E.bracket(x[0], x[3]), E.bracket(x[1], x[2]))) for x in pts]
    rkT = vlie.rank_mod(np.array(T, dtype=np.int64), p)
    gens4 = [P[i][c] for i in range(3) for c in range(2)]
    out = ['%-16s p=%d: rank of (T1, T3) values at %d random points = %d' % (L.name, p, npts, rkT)]
    # explicit evaluation of d on C^3 basis at random points (momenta numeric)
    for a in range(3):
        for b in range(3):
            B3 = C3_basis(a, b); B2 = C2_basis(a, b)
            rowsd3 = []
            for _ in range(npts):
                X = [rv() for _ in range(4)]
                mom = [tuple(int(t) for t in rng.integers(0, p, 2)) for _ in range(3)]
                mom.append(tuple((-sum(m[c] for m in mom)) % p for c in range(2)))
                row = []
                for e in B3:
                    tot = 0
                    for i, j in itertools.combinations(range(4), 2):
                        k, l = [r for r in range(4) if r not in (i, j)]
                        q = [((mom[i][0] + mom[j][0]) % p, (mom[i][1] + mom[j][1]) % p), mom[k], mom[l]]
                        Pv = int(sp.Poly(e, *[P[0][0], P[0][1], P[1][0], P[1][1]]).eval({P[0][0]: q[0][0], P[0][1]: q[0][1], P[1][0]: q[1][0], P[1][1]: q[1][1]}) % p) if e.free_symbols else int(e) % p
                        val = fval(E.bracket(X[i], X[j]), X[k], X[l])
                        tot = (tot + (-1) ** (i + j) * Pv * val) % p
                    row.append(tot)
                rowsd3.append(row)
            rk3 = vlie.rank_mod(np.array(rowsd3, dtype=np.int64), p) if B3 else 0
            rk2 = 0
            for e in B2:
                X = [rv() for _ in range(3)]
                mom = [tuple(int(t) for t in rng.integers(0, p, 2)) for _ in range(2)]
                mom.append(tuple((-sum(m[c] for m in mom)) % p for c in range(2)))
                tot = 0
                for i, j in itertools.combinations(range(3), 2):
                    k = [r for r in range(3) if r not in (i, j)][0]
                    q0 = ((mom[i][0] + mom[j][0]) % p, (mom[i][1] + mom[j][1]) % p)
                    Pv = int(sp.Poly(e, P[0][0], P[0][1]).eval({P[0][0]: q0[0], P[0][1]: q0[1]}) % p)
                    tot = (tot + (-1) ** (i + j) * Pv * E.kap(E.bracket(X[i], X[j]), X[k])) % p
                rk2 += 1 if tot else 0
            Z3 = len(B3) - rk3
            out.append('   (%d,%d): dim C^3 %d, rank d3 (evaluation at %d points) %d, rank d2 %d, Z^3 %d, H^3 %d'
                       % (a, b, len(B3), npts, rk3, rk2, Z3, Z3 - rk2))
    # d E_f at a random point, and the T1/T3 formula
    X = [rv() for _ in range(4)]
    mom = [tuple(int(t) for t in rng.integers(0, p, 2)) for _ in range(3)]
    mom.append(tuple((-sum(m[c] for m in mom)) % p for c in range(2)))
    Ef = lambda q: 2 * sum(qq[0] * qq[1] for qq in q)
    tot = 0
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        q = [((mom[i][0] + mom[j][0]) % p, (mom[i][1] + mom[j][1]) % p), mom[k], mom[l]]
        tot = (tot + (-1) ** (i + j) * Ef(q) * fval(E.bracket(X[i], X[j]), X[k], X[l])) % p
    Bf = lambda u, v: (u[0] * v[1] + u[1] * v[0]) % p
    sb = lambda u, v: ((u[0] - v[0]) % p, (u[1] - v[1]) % p)
    T1 = E.kap(E.bracket(X[0], X[1]), E.bracket(X[2], X[3])); T3 = E.kap(E.bracket(X[0], X[3]), E.bracket(X[1], X[2]))
    form = 2 * (Bf(sb(mom[0], mom[3]), sb(mom[2], mom[1])) * T1 + Bf(sb(mom[0], mom[1]), sb(mom[2], mom[3])) * T3) % p
    out.append('   d E_f at a random point = %d (nonzero: %s) ; equals 2[B(p0-p3,p2-p1)T1 + B(p0-p1,p2-p3)T3]: %s' % (tot, tot != 0, tot == form))
    return '\n'.join(out)


if __name__ == '__main__':
    if sys.argv[1] == 'universal':
        res, rows = universal()
        print('Universal table (T = f only; ranks exact over Q):')
        print('\n'.join(rows))
        ok1, ok2, c1, c3, altok = dEf_universal()
        print('d E_f: alt(q0^1 q0^2 f) = 2 sum q_i^1 q_i^2 f: %s ; d E_f = 2[B(p0-p3,p2-p1) T1 + B(p0-p1,p2-p3) T3]: %s ; '
              'sum formula: %s ; c1 = %s ; c3 = %s' % (altok, ok1, ok2, c1, c3))
    else:
        for key in sys.argv[1:]:
            print(numeric(key), flush=True)
