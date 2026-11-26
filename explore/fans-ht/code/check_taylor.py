"""(b) Taylor coefficients at endpoint sources.

 T1  polynomial currents mu_ik, nu_i (S1 eq-1-38 at lam_2..lam_n = 0) = Kirchhoff currents, and
     -sum_e y_e z_e = -lam_1 z_{1,O} - lam_{n+1} z_{n+1,O} - sum_i y_{i,i+1} Z_i   (random rational r, z).
 T2  order-m part of K(lam; z): coefficient formula (taylor.coeff_11, assembled over all multi-indices)
     = direct exact integral of density * (-sum_e y_e z_e)^m / m!  (symbolic z), n = 2, 4.
 T3  engine (htk spred.taylor, its own coordinates) = formula, n = 2 (orders 0..3), n = 4 (orders 0..2).
 T4  (2,0): taylor.coeff_20 (S1 eq-sum-N4-coeff) = direct integral of S1 eq-1-41 density, n = 2, 3.
"""
import sys
import math
import random
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import sympy as sp
import fanlib as F
import taylor as TY
from integrate_general import beta2

rng = random.Random(7)


def symbols(n):
    u = {k: sp.Symbol("u%d" % k) for k in range(1, n + 1)}
    v = {k: sp.Symbol("v%d" % k) for k in range(1, n + 1)}
    V = {k: sp.Symbol("V%d" % k) for k in range(1, n + 1)}
    return u, v, V


def mu_nu(n, u, v, V):
    mu, nu = {}, {}
    for i in range(1, n + 1):
        for k in range(i, n + 1):
            e = sp.Mul(*[V[m] for m in range(i, k)]) * sp.Mul(*[u[l] for l in range(1, k + 1)])
            mu[(i, k)] = e * (v[k] if k < n else 1)
        nu[i] = v[n] * sp.Mul(*[V[m] for m in range(i, n)])
    return mu, nu


def edge_z(n):
    """z symbols on the engine edge order a_1, b_1, ..., a_n, b_n, b_{n+1}."""
    za = {i: sp.Symbol("za%d" % i) for i in range(1, n + 1)}
    zb = {i: sp.Symbol("zb%d" % i) for i in range(1, n + 2)}
    zlist = []
    for i in range(1, n + 1):
        zlist += [za[i], zb[i]]
    zlist.append(zb[n + 1])
    return za, zb, zlist


# ---------------------------------------------------------------- T1
for n in (2, 4, 6):
    u, v, V = symbols(n)
    mu, nu = mu_nu(n, u, v, V)
    for _ in range(3):
        uu, vv = F.random_r(n, rng)
        t = F.t_of_r(n, uu, vv)
        l1, ln1 = Fr(rng.randint(-9, 9)), Fr(rng.randint(-9, 9))
        lam = [l1] + [Fr(0)] * (n - 1) + [ln1]
        y, _, _ = F.kirchhoff(n, t, lam)
        sub = {u[k]: sp.Rational(uu[k].numerator, uu[k].denominator) for k in u}
        sub.update({v[k]: sp.Rational(vv[k].numerator, vv[k].denominator) for k in v})
        sub.update({V[k]: 1 - sub[v[k]] for k in V})
        for i in range(1, n + 1):
            ya = sum(mu[(i, k)] for k in range(i, n + 1)).subs(sub) * l1 - nu[i].subs(sub) * ln1
            assert sp.Rational(ya) == sp.Rational(y[2 * (i - 1)].numerator, y[2 * (i - 1)].denominator)
        z = [Fr(rng.randint(-5, 5)) for _ in range(2 * n + 1)]
        lhs = -sum(ye * ze for ye, ze in zip(y, z))
        za = {i: z[2 * (i - 1)] for i in range(1, n + 1)}
        zb = {i: z[2 * (i - 1) + 1] for i in range(1, n + 1)}
        zb[n + 1] = z[2 * n]
        rhs = -l1 * zb[1] - ln1 * zb[n + 1] - sum(y[2 * (i - 1)] * (za[i] + zb[i + 1] - zb[i]) for i in range(1, n + 1))
        assert lhs == rhs
print("T1: currents (eq-1-38 at endpoint sources) and exponent decomposition hold, n = 2, 4, 6")


# ---------------------------------------------------------------- direct integration helper
def integrate_poly_11(n, expr, u, v, V):
    """int over the S1 domain of (n-1)!!/2^{n/2} * mono(r) * expr, expr polynomial in u, v, V
    (V = 1 - v) and other symbols; mono = v_n^{n/2-1} prod_{k<n} v_k^{1/2} V_k^{k/2-1} prod u_m^{(n-1-m)/2}."""
    gens = [u[k] for k in range(1, n + 1)] + [v[k] for k in range(1, n + 1)] + [V[k] for k in range(1, n + 1)]
    P = sp.Poly(sp.expand(expr), *gens)
    out = 0
    for mon, c in P.terms():
        eu = mon[:n]
        ev = mon[n:2 * n]
        eV = mon[2 * n:]
        q, k = Fr(TY.dfact(n - 1), 2 ** (n // 2)), 0
        for j in range(1, n):
            a2 = (n - 1 - j) + 2 * eu[j - 1]          # doubled exponent of u_j
            q *= Fr(2, a2 + 2)
            dq, dk = beta2(1 + 2 * ev[j - 1] + 2, (j - 2) + 2 * eV[j - 1] + 2)
            q *= dq
            k += dk
        assert eV[n - 1] == 0
        a2 = -1 + 2 * eu[n - 1]
        c2 = (n - 2) + 2 * ev[n - 1]
        dq, dk = beta2(a2 + 2, c2 + 4)                # int_{Delta^2} u^a v^c = B(a+1, c+2)/(c+1)
        q *= dq / Fr(c2 + 2, 2)
        k += dk
        out += c * sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k
    return sp.expand(out)


def formula_order_m(n, m, l1, ln1, za, zb):
    Z = {i: za[i] + zb[i + 1] - zb[i] for i in range(1, n + 1)}
    tot = 0
    for r0, s0, r, s in TY.multi_indices(n, m):
        q, k = TY.coeff_11(n, r, s)
        term = (-l1 * zb[1]) ** r0 / sp.factorial(r0) * (-ln1 * zb[n + 1]) ** s0 / sp.factorial(s0)
        for (i, kk), p in r.items():
            term *= (-l1 * Z[i]) ** p / sp.factorial(p)
        for i, p in s.items():
            term *= (ln1 * Z[i]) ** p / sp.factorial(p)
        tot += term * sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k
    return sp.expand((l1 * ln1) ** (n // 2) * tot)


# ---------------------------------------------------------------- T2, T3
from htk import fan, sp_reduce
for n, mmax_direct, mmax_engine in ((2, 3, 3), (4, 3, 2)):
    u, v, V = symbols(n)
    mu, nu = mu_nu(n, u, v, V)
    za, zb, zlist = edge_z(n)
    l1, ln1 = sp.Integer(3), sp.Integer(-2)
    Z = {i: za[i] + zb[i + 1] - zb[i] for i in range(1, n + 1)}
    E = -l1 * zb[1] - ln1 * zb[n + 1] - sum((sum(mu[(i, k)] for k in range(i, n + 1)) * l1 - nu[i] * ln1) * Z[i] for i in range(1, n + 1))
    red = sp_reduce(fan(n), D=2, prefer="large")
    lam_eng = [[l1]] + [[0]] * (n - 1) + [[ln1]]
    eng = red.taylor(lam_eng, [[x] for x in zlist], mmax_engine, T=1, check=1)
    for m in range(0, mmax_direct + 1):
        direct = integrate_poly_11(n, (l1 * ln1) ** (n // 2) * E ** m / sp.factorial(m), u, v, V)
        form = formula_order_m(n, m, l1, ln1, za, zb)
        assert sp.expand(direct - form) == 0, ("T2", n, m)
        if m <= mmax_engine:
            e = sp.sympify(eng[m])
            assert sp.expand(e - form) == 0, ("T3", n, m, e, form)
        print("n=%d order %d: formula = direct integral%s; %d monomials in z" % (
            n, m, " = engine" if m <= mmax_engine else "", len(sp.Add.make_args(form))))
    print("  n=%d, order 1 term:" % n, sp.factor(formula_order_m(n, 1, sp.Symbol("l1"), sp.Symbol("lN"), za, zb)))

# ---------------------------------------------------------------- T4 (2,0) control of the format
for n in (2, 3):
    u, v, V = symbols(n)
    for m in range(0, 3):
        for r0, s0, r, s in TY.multi_indices(n, m):
            mu, nu = mu_nu(n, u, v, V)
            expr = sp.Integer(1)
            for (i, k), p in r.items():
                expr *= mu[(i, k)] ** p
            for i, p in s.items():
                expr *= nu[i] ** p
            mono = v[n] ** (n - 1) * sp.Mul(*[u[j] ** (n - j) * v[j] * V[j] ** (j - 1) for j in range(1, n)])
            P = sp.Poly(sp.expand(expr * mono), *([u[k] for k in u] + [v[k] for k in v] + [V[k] for k in V]))
            tot = Fr(0)
            for mon, c in P.terms():
                q = Fr(int(c))
                for j in range(1, n):
                    q *= Fr(1, mon[j - 1] + 1)
                    a, b = mon[n + j - 1], mon[2 * n + j - 1]
                    q *= Fr(math.factorial(a) * math.factorial(b), math.factorial(a + b + 1))
                a, c_ = mon[n - 1], mon[2 * n - 1]
                q *= Fr(math.factorial(a) * math.factorial(c_), math.factorial(a + c_ + 2))
                tot += q
            assert tot == TY.coeff_20(n, r, s), (n, r, s)
print("T4: S1 eq-sum-N4-coeff reproduced by direct integration of eq-1-41, n = 2, 3, orders 0..2")
print("all Taylor checks passed")
