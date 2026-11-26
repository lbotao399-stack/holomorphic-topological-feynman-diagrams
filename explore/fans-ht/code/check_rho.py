"""rho_n(M) of the fan.

Triangle cycles c_r = a_r + b_{r+1} - b_r give  psi^tau_r = (1/2)[tau_r alpha_r - tau_{r+1} beta_r],
alpha_r = d log(a_r/b_r), beta_r = d log(a_r/b_{r+1});  <tau_u tau_v> = 2 G_uv, G = Lambda^{-1} on the
path vertices, G_uv = g_min(u,v) h_max(u,v).  Claim:

  Pf(-1/2 dM M^{-1} dM) = 2^{-n/2} sum_{matchings pi of 1..n} wedge_{r=1}^n C^pi_r,
  C_r = A_r = g_r alpha_r - g_{r+1} beta_r  (r the smaller element of its pair),
  C_r = B_r = h_r alpha_r - h_{r+1} beta_r  (r the larger element),
  rho_n(M) = (2^n U)^{-1/2} sum_pi wedge_r C^pi_r ;   g_i h_{i+1} - g_{i+1} h_i = -a_i.

 R1  identity of forms in the dt_e at random rational t (exact), n = 2, 4, 6; semiseparability of G and
     the Wronskian identity.
 R2  the same identity pulled back to the S1 coordinates with the explicit r-expressions
     (U monomial, g_u = 1/P_u, g_{n+1} = (1-v_n)/(P_n u_n), h_v = P_v^{-1} sum_{k>=v} P_k^2 c_k,
     h_{n+1} = P_n u_n v_n, alpha_r, beta_r as d log of the ratios below), n = 2, 4, exact:
       a_r/b_r = (1-u_r) v_{r-1}/u_r  (v_0 := 1),  a_r/b_{r+1} = (1-u_r) v_r/(1-v_r)    (r < n),
       a_n/b_n = (1-u_n-v_n) v_{n-1}/u_n,          a_n/b_{n+1} = (1-u_n-v_n)/v_n.
"""
import sys
import random
import itertools
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import sympy as sp
import fanlib as F

rng = random.Random(31)


def brute_pf(n, t):
    N = 2 * n + 1
    C = F.cycles_triangles(n)
    # use the plain triangle basis (no sign flip): compare with the same basis below
    C = []
    for i in range(1, n + 1):
        row = [0] * N
        row[2 * (i - 1)] = 1
        row[2 * i + 1 if i < n else 2 * n] = 1
        row[2 * (i - 1) + 1] = -1
        C.append(row)
    M = [[sum(Fr(C[r][e] * C[s][e]) * t[e] for e in range(N)) for s in range(n)] for r in range(n)]
    Minv = F.finv(M)
    dM = [[{(e,): Fr(C[r][e] * C[s][e]) for e in range(N) if C[r][e] * C[s][e] != 0} for s in range(n)] for r in range(n)]
    A = [[{} for _ in range(n)] for _ in range(n)]
    for r in range(n):
        for s in range(r + 1, n):
            acc = {}
            for a in range(n):
                for b in range(n):
                    if Minv[a][b] != 0:
                        acc = F.fadd(acc, F.wedge(dM[r][a], dM[s][b]), -Minv[a][b] / 2)
            A[r][s] = acc
            A[s][r] = {k: -x for k, x in acc.items()}
    return F.pfaff(A, list(range(n))), F.fdet(M)


def formula_pf(n, alpha, beta, g, h, wedge, add, scale, one):
    tot = None
    for m in F.matchings(list(range(1, n + 1))):
        op = {i for (i, j) in m}
        w = one
        for r in range(1, n + 1):
            if r in op:
                c = add(scale(alpha[r], g[r]), scale(beta[r], -g[r + 1]))
            else:
                c = add(scale(alpha[r], h[r]), scale(beta[r], -h[r + 1]))
            w = wedge(w, c)
        tot = w if tot is None else add(tot, w)
    return scale(tot, Fr(1, 2 ** (n // 2)))


# ---------------------------------------------------------------- R1
for n in (2, 4, 6):
    for _ in range(2):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 13)) for _ in range(2 * n + 1)]
        pf, U = brute_pf(n, t)
        a = {i: t[2 * (i - 1)] for i in range(1, n + 1)}
        b = {i: t[2 * (i - 1) + 1] for i in range(1, n + 1)}
        b[n + 1] = t[2 * n]
        ia = {i: 2 * (i - 1) for i in range(1, n + 1)}
        ib = {i: 2 * (i - 1) + 1 for i in range(1, n + 1)}
        ib[n + 1] = 2 * n
        _, _, L = F.kirchhoff(n, t, [0] * (n + 1))
        G = F.finv(L)                           # vertices 0..n = S1 1..n+1
        Gs = lambda u, v: G[u - 1][v - 1]
        h = {v: Gs(1, v) for v in range(1, n + 2)}
        g = {u: Gs(u, n + 1) / Gs(1, n + 1) for u in range(1, n + 2)}
        for u in range(1, n + 2):
            for v in range(u, n + 2):
                assert g[u] * h[v] == Gs(u, v), "semiseparable"
        for i in range(1, n + 1):
            assert g[i] * h[i + 1] - g[i + 1] * h[i] == -a[i], "Wronskian"
        alpha = {r: F.fadd({(ia[r],): 1 / a[r]}, {(ib[r],): -1 / b[r]}) for r in range(1, n + 1)}
        beta = {r: F.fadd({(ia[r],): 1 / a[r]}, {(ib[r + 1],): -1 / b[r + 1]}) for r in range(1, n + 1)}
        fp = formula_pf(n, alpha, beta, g, h, F.wedge, F.fadd, lambda f, c: {k: c * x for k, x in f.items() if c * x != 0}, {(): Fr(1)})
        assert fp == pf, ("R1", n)
print("R1: Pf(-1/2 dM M^-1 dM) = 2^{-n/2} sum_pi wedge_r C_r as forms in dt (exact), G semiseparable, g_i h_{i+1} - g_{i+1} h_i = -a_i; n = 2, 4, 6")

# ---------------------------------------------------------------- R2
for n in (2, 4):
    us = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    vs = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    coords = []
    for k in range(1, n + 1):
        coords += [us[k], vs[k]]
    texpr = F.t_of_r(n, us, vs)
    # explicit r-expressions
    P = {1: sp.Integer(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * us[k - 1]
    R = {n: sp.Integer(1)}
    for k in range(n - 1, 0, -1):
        R[k] = R[k + 1] * us[k + 1] / (1 - vs[k])
    c = {k: us[k] * (1 - us[k]) * R[k] for k in range(1, n + 1)}
    gr = {u: 1 / P[u] for u in range(1, n + 1)}
    gr[n + 1] = (1 - vs[n]) / (P[n] * us[n])
    hr = {v: sum(P[k] ** 2 * c[k] for k in range(v, n + 1)) / P[v] for v in range(1, n + 1)}
    hr[n + 1] = P[n] * us[n] * vs[n]
    Ur = sp.Mul(*[us[k + 1] ** k * (1 - vs[k]) ** (-k) / vs[k] for k in range(1, n)])
    ratA = {r: (1 - us[r]) * (vs[r - 1] if r > 1 else 1) / us[r] for r in range(1, n)}
    ratB = {r: (1 - us[r]) * vs[r] / (1 - vs[r]) for r in range(1, n)}
    ratA[n] = (1 - us[n] - vs[n]) * (vs[n - 1] if n > 1 else 1) / us[n]
    ratB[n] = (1 - us[n] - vs[n]) / vs[n]
    # check the ratios against the chart
    for r in range(1, n + 1):
        assert sp.simplify(texpr[2 * (r - 1)] / texpr[2 * (r - 1) + 1] - ratA[r]) == 0
        nb = 2 * r + 1 if r < n else 2 * n
        assert sp.simplify(texpr[2 * (r - 1)] / texpr[nb] - ratB[r]) == 0
    for _ in range(2):
        uu, vv = F.random_r(n, rng)
        sub = {us[k]: sp.Rational(uu[k].numerator, uu[k].denominator) for k in us}
        sub.update({vs[k]: sp.Rational(vv[k].numerator, vv[k].denominator) for k in vs})
        tq = [Fr(int(sp.Rational(x.subs(sub)).p), int(sp.Rational(x.subs(sub)).q)) for x in texpr]
        pf, U = brute_pf(n, tq)
        assert sp.Rational(Ur.subs(sub)) == sp.Rational(U.numerator, U.denominator), "U monomial"
        # chart Jacobian dt_e/dr_j
        Jt = [[sp.Rational(sp.diff(te, x).subs(sub)) for x in coords] for te in texpr]
        # pull back pf (forms in dt) to forms in dr
        pull = {}
        for key, val in pf.items():
            for S in itertools.combinations(range(2 * n), n):
                d = sp.Matrix([[Jt[e][j] for j in S] for e in key]).det()
                if d != 0:
                    pull[S] = pull.get(S, 0) + sp.Rational(val.numerator, val.denominator) * d
        pull = {k: v for k, v in pull.items() if v != 0}
        # formula in dr
        def dlog(expr):
            return {(j,): sp.Rational(sp.diff(sp.log(expr), x).subs(sub)) for j, x in enumerate(coords) if sp.diff(sp.log(expr), x).subs(sub) != 0}
        alpha = {r: dlog(ratA[r]) for r in range(1, n + 1)}
        beta = {r: dlog(ratB[r]) for r in range(1, n + 1)}
        gq = {k: sp.Rational(x.subs(sub)) for k, x in gr.items()}
        hq = {k: sp.Rational(x.subs(sub)) for k, x in hr.items()}
        fp = formula_pf(n, alpha, beta, gq, hq, F.wedge, F.fadd, lambda f, cc: {k: cc * x for k, x in f.items() if cc * x != 0}, {(): sp.Integer(1)})
        diff = F.fadd(fp, pull, -1)
        assert all(sp.simplify(v) == 0 for v in diff.values()), ("R2", n)
print("R2: the r-coordinate expression of Pf(-1/2 dM M^-1 dM) equals the pullback of the brute-force form, n = 2, 4; U = prod_k u_{k+1}^k (1-v_k)^{-k} v_k^{-1}")
