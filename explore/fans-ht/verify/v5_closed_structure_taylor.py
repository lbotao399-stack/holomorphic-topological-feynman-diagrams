"""V8/V9 closed forms (a7, d1 integrals), ensemble form (a8), structure (a11), reflection (a12),
Taylor coefficients (b1).  Exact arithmetic; own Beta/Gamma code (v4_integrate.gamma_q)."""
import itertools
import math
import random
import time
from fractions import Fraction as Fr

import sympy as sp

from vlib import (fan_graph, fan_cycles, t_of_r, r_of_t, rand_r, jac_det, brute_density, claim_density,
                  claim_LM, claim_gh, half_power, same_halfint, currents_mesh, perm_sign, ia, ib, haf_LM)
from v4_integrate import gamma_q, beta_q, kernels_file

rng = random.Random(5)


def dfact(k):
    out = 1
    while k > 1:
        out *= k
        k -= 2
    return out


# ---------------------------------------------------------------- a7, d1: closed forms
print("a7 / d1 closed forms")
for n in range(1, 31):
    # (2,0): int v_n^{n-1} prod_{j<n} u_j^{n-j} v_j (1-v_j)^{j-1}
    q = Fr(1)
    for j in range(1, n):
        q *= Fr(1, n - j + 1) * beta_q(2, j)[0]
    g1, _ = gamma_q(1)
    g2, _ = gamma_q(n)
    g3, _ = gamma_q(n + 2)
    q *= g1 * g2 / g3
    assert q == Fr(1, math.factorial(n) ** 3 * (n + 1)), n
    if n % 2:
        continue
    # (1,1) endpoint: (n-1)!!/2^{n/2} int v_n^{n/2-1} prod_{k<n} v_k^{1/2}(1-v_k)^{k/2-1} prod_m u_m^{(n-1-m)/2}
    qq, kk = Fr(dfact(n - 1), 2 ** (n // 2)), 0
    for m in range(1, n):
        qq *= 1 / (Fr(n - 1 - m, 2) + 1)
    for k in range(1, n):
        b, kb = beta_q(Fr(3, 2), Fr(k, 2))
        qq *= b
        kk += kb
    a1, k1 = gamma_q(Fr(1, 2))
    a2, k2 = gamma_q(Fr(n, 2))
    a3, k3 = gamma_q(Fr(n, 2) + Fr(3, 2))
    qq *= a1 * a2 / a3
    kk += k1 + k2 - k3
    val = sp.Rational(qq.numerator, qq.denominator) * sp.sqrt(sp.pi) ** kk
    cf1 = (4 * sp.pi) ** (n // 2) / ((n + 1) * sp.factorial(n) ** 2)
    cf2 = sp.pi ** (n // 2 + 1) / (2 ** (n + 1) * sp.gamma(sp.Rational(n, 2) + 1) ** 2 * sp.gamma(sp.Rational(n + 1, 2)) * sp.gamma(sp.Rational(n + 3, 2)))
    assert sp.simplify(val - cf1) == 0 and sp.simplify(val - cf2) == 0, n
    if n <= 10:
        print("  n=%2d: int g_end = %s" % (n, sp.simplify(val)))
print("  (2,0): int eq-1-41 density = 1/((n!)^3 (n+1)) exact, n = 1..30;  (1,1) endpoint: own Beta product = (4 pi)^{n/2}/((n+1)(n!)^2)"
      " = pi^{n/2+1}/(2^{n+1} Gamma(n/2+1)^2 Gamma((n+1)/2) Gamma((n+3)/2)), n = 2..30 even")

# ---------------------------------------------------------------- a8 ensemble form (pointwise)
print("a8 ensemble form: g_n = w(r) Haf_n(l, m)/(n-1)!!, w = endpoint density at lam_1 = lam_{n+1} = 1")
for n in (2, 4, 6, 8):
    for p in range(3):
        u, v = rand_r(n, rng)
        lam = {i: Fr(rng.randint(-7, 7)) for i in range(1, n + 2)}
        g = claim_density(n, u, v, lam)
        L, M = claim_LM(n, u, v, lam)
        P = {1: Fr(1)}
        for k in range(2, n + 2):
            P[k] = P[k - 1] * u[k - 1]
        mm = {}
        for j in range(1, n + 1):
            s = Fr(0)
            for k in range(j + 1, n + 1):
                lamjk = sum((P[k] / P[w] * lam[w] for w in range(j + 1, k + 1)), Fr(0))
                den = Fr(1)
                for m_ in range(k, n):
                    den *= 1 - v[m_]
                s += (1 - u[k]) * lamjk / den
            mm[j] = lam[n + 1] + s / v[n]
            assert mm[j] * P[n] * u[n] * v[n] == M[j]
        pairs = [(v[n], n - 2)]
        for k in range(1, n):
            pairs += [(v[k], 1), (1 - v[k], k - 2)]
        for m_ in range(1, n + 1):
            pairs.append((u[m_], n - 1 - m_))
        wr, wrad = half_power(pairs)
        wr = wr * dfact(n - 1) / 2 ** (n // 2)
        rhs = wr * haf_LM(n, L, mm) / dfact(n - 1)
        assert same_halfint(g[0], g[1], rhs, wrad) or (g[0] == 0 and rhs == 0)
    # normalization of the ensemble: prod of Beta/Dirichlet normalizers = int w
    print("  n=%d: m_j P_n u_n v_n = M_j and g = w Haf(l,m)/(n-1)!! exact at 3 points" % n)

# ---------------------------------------------------------------- a11 structure (symbolic, n = 2, 4)
print("a11 structure")
from v4_integrate import build_p, integrate_p
KF = kernels_file()
for n in (2, 4):
    ls = {i: sp.Symbol("l%d" % i) for i in range(1, n + 2)}
    s_ = sp.Symbol("s")
    p, us, vs = build_p(n, ls)
    q = sp.cancel(p / ls[1])
    q2 = sp.expand(sp.cancel(q.subs(ls[2], s_ - us[1] * ls[1])))
    dep = q2.has(us[1]) or q2.has(ls[1])
    print("  n=%d: (g/m_0)/lam_1 with lam_2 = s - u_1 lam_1 is free of u_1 and lam_1: %s (m_0 is free of u_1)" % (n, not dep))
for n in (2, 4, 6, 8):
    ls = {i: sp.Symbol("l%d" % i) for i in range(1, n + 2)}
    K = KF[n]
    Kc = (4 * sp.pi) ** (n // 2) / ((n + 1) * sp.factorial(n) ** 2)
    sub = {ls[i]: 0 for i in range(3, n + 1)}
    K3 = sp.expand(K.subs(sub))
    f3 = sp.expand(Kc * ls[n + 1] ** (n // 2) * ((ls[1] + ls[2]) ** (n // 2) - ls[2] ** (n // 2)))
    div = sp.cancel(K / (ls[1] * ls[n + 1]))
    poly_div = sp.fraction(sp.together(div))[1].free_symbols == set()
    print("  n=%d (general_kernels.txt): lam_1 lam_{n+1} divides K: %s; K(lam_1, lam_2, 0.., lam_{n+1}) = C_n lam_{n+1}^{n/2}[(lam_1+lam_2)^{n/2} - lam_2^{n/2}]: %s"
          % (n, poly_div, sp.simplify(K3 - f3) == 0))

# ---------------------------------------------------------------- a12 reflection (pointwise, brute force)
print("a12 reflection")


def refl_perm(n):
    """pi(p) = old position whose time sits at new position p:  t'_p = t_{pi(p)}."""
    pi = [None] * (2 * n + 1)
    for j in range(1, n + 1):
        pi[ia(j)] = ia(n + 1 - j)
    for j in range(1, n + 2):
        pi[ib(j, n)] = ib(n + 2 - j, n)
    return pi


for n in (2, 4, 6, 8):
    G, C = fan_cycles(n, "pathtree")
    pi = refl_perm(n)
    sig = []
    for p in range(2 if n < 8 else 1):
        t = [Fr(rng.randint(1, 40), rng.randint(1, 40)) for _ in range(G.N)]
        lam = [Fr(rng.choice([-1, 1]) * rng.randint(1, 6)) for _ in range(n + 1)]
        tp = [t[pi[q]] for q in range(G.N)]
        lamp = list(reversed(lam))
        D1, U1 = brute_density(G, C, t, [lam], 1)
        D2, U2 = brute_density(G, C, tp, [lamp], 1)
        assert U1 == U2
        assert D1 != 0
        sig.append(str(D2 / D1))
    print("  n=%d: D(phi(t); lam reversed) / D(t; lam) = %s at random rational t (brute force);  sgn(pi_E) = %+d" % (n, sig, perm_sign(pi)))
print("  sgn(pi_E) for even n = 2..40:", sorted(set(perm_sign(refl_perm(n)) for n in range(2, 41, 2))))
print("  sgn(pi_E) for odd  n = 1..39:", [perm_sign(refl_perm(n)) for n in range(1, 40, 2)][:6], "...")
for n in (2, 4, 6, 8):
    ls = {i: sp.Symbol("l%d" % i) for i in range(1, n + 2)}
    K = KF[n]
    Kr = K.subs({ls[i]: sp.Symbol("m%d" % i) for i in ls}, simultaneous=True).subs({sp.Symbol("m%d" % i): ls[n + 2 - i] for i in ls}, simultaneous=True)
    print("  n=%d: K(lam_{n+1},..,lam_1) - K(lam_1,..,lam_{n+1}) = %s (general_kernels.txt)" % (n, sp.simplify(Kr - K)))

# ---------------------------------------------------------------- b1 Taylor
print("b1 Taylor (endpoint sources)")


def mu_nu(n, u, v):
    mu, nu = {}, {}
    for i in range(1, n + 1):
        for k in range(i, n + 1):
            x = Fr(1)
            for m in range(i, k):
                x *= 1 - v[m]
            for l in range(1, k + 1):
                x *= u[l]
            if k < n:
                x *= v[k]
            mu[(i, k)] = x
        y = v[n]
        for m in range(i, n):
            y *= 1 - v[m]
        nu[i] = y
    return mu, nu


for n in (2, 4, 6):
    G = fan_graph(n)
    for p in range(3):
        u, v = rand_r(n, rng)
        t = t_of_r(n, u, v)
        l1, ln1 = Fr(rng.randint(-5, 5)), Fr(rng.randint(-5, 5))
        lam = [l1] + [Fr(0)] * (n - 1) + [ln1]
        y = currents_mesh(G, t, lam)
        mu, nu = mu_nu(n, u, v)
        for i in range(1, n + 1):
            assert y[ia(i)] == sum(mu[(i, k)] for k in range(i, n + 1)) * l1 - nu[i] * ln1
        z = [Fr(rng.randint(-9, 9), rng.randint(1, 9)) for _ in range(G.N)]
        lhs = sum(y[e] * z[e] for e in range(G.N))
        Z = {i: z[ia(i)] + z[ib(i + 1, n)] - z[ib(i, n)] for i in range(1, n + 1)}
        rhs = l1 * z[ib(1, n)] + ln1 * z[ib(n + 1, n)] + sum(y[ia(i)] * Z[i] for i in range(1, n + 1))
        assert lhs == rhs
    print("  n=%d: y_{i,i+1} = sum_k mu_ik lam_1 - nu_i lam_{n+1} (mesh solver) and sum_e y_e z_e decomposition: exact at 3 points" % n)


def c_formula(n, rr, ss):
    """c^{(n)} of claim b1 as (q, k): q sqrt(pi)^k."""
    q, k = Fr(dfact(n - 1), 2 ** (n // 2)), 0
    Rn = sum(rr.get((i, n), 0) for i in range(1, n + 1))
    S = sum(ss.get(i, 0) for i in range(1, n + 1))
    a, ka = gamma_q(Fr(1, 2) + Rn)
    b, kb = gamma_q(Fr(n, 2) + S)
    c, kc = gamma_q(Fr(n + 3, 2) + Rn + S)
    q *= a * b / c
    k += ka + kb - kc
    for j in range(1, n):
        X = sum(rr.get((i, kk), 0) for kk in range(j, n + 1) for i in range(1, kk + 1))
        q *= Fr(2, n + 1 - j + 2 * X)
    for j in range(1, n):
        rj = sum(rr.get((i, j), 0) for i in range(1, j + 1))
        Y = sum(ss.get(i, 0) + sum(rr.get((i, kk), 0) for kk in range(j + 1, n + 1)) for i in range(1, j + 1))
        a, ka = gamma_q(Fr(3, 2) + rj)
        b, kb = gamma_q(Fr(j, 2) + Y)
        c, kc = gamma_q(Fr(j + 3, 2) + rj + Y)
        q *= a * b / c
        k += ka + kb - kc
    return q, k


def c_direct(n, rr, ss):
    """int over the domain of g_end/(lam_1 lam_{n+1})^{n/2} prod mu_ik^{r_ik} prod nu_i^{s_i}: exponents
    collected from the monomials mu, nu, g_end and integrated with own Beta code."""
    eu = {m: Fr(n - 1 - m, 2) for m in range(1, n + 1)}
    ev = {k: Fr(1, 2) for k in range(1, n)}
    ew = {k: Fr(k, 2) - 1 for k in range(1, n)}      # (1 - v_k)
    evn = Fr(n, 2) - 1
    for (i, k), r_ in rr.items():
        for m in range(i, k):
            ew[m] += r_
        for l in range(1, k + 1):
            eu[l] += r_
        if k < n:
            ev[k] += r_
    for i, s_ in ss.items():
        evn += s_
        for m in range(i, n):
            ew[m] += s_
    q, kk = Fr(dfact(n - 1), 2 ** (n // 2)), 0
    for m in range(1, n):
        q *= 1 / (eu[m] + 1)
    for k in range(1, n):
        b, kb = beta_q(ev[k] + 1, ew[k] + 1)
        q *= b
        kk += kb
    a, ka = gamma_q(eu[n] + 1)
    b, kb = gamma_q(evn + 1)
    c, kc = gamma_q(eu[n] + evn + 3)
    q *= a * b / c
    kk += ka + kb - kc
    return q, kk


for n in (2, 4, 6):
    cnt = 0
    keys = [(i, k) for k in range(1, n + 1) for i in range(1, k + 1)]
    for trial in range(60):
        rr = {key: rng.choice([0, 0, 0, 1, 2]) for key in keys}
        ss = {i: rng.choice([0, 0, 1, 2]) for i in range(1, n + 1)}
        if trial == 0:
            rr = {key: 0 for key in keys}
            ss = {i: 0 for i in range(1, n + 1)}
        cf = c_formula(n, rr, ss)
        cd = c_direct(n, rr, ss)
        assert cf == cd, (n, rr, ss, cf, cd)
        cnt += 1
    q0, k0 = c_formula(n, {}, {})
    print("  n=%d: c^(n) formula = direct Beta integration of g_end prod mu^r nu^s for %d multi-indices (incl. 0); c_0 = %s" %
          (n, cnt, sp.simplify(sp.Rational(q0.numerator, q0.denominator) * sp.sqrt(sp.pi) ** k0)))

# order-1 coefficients n = 2, 4 from int g_end y_e (y_e from mu, nu; spokes by conservation)
for n in (2, 4):
    l1, l3 = sp.symbols("l1 lN")
    za = {i: sp.Symbol("z_a%d" % i) for i in range(1, n + 1)}
    zb = {i: sp.Symbol("z_b%d" % i) for i in range(1, n + 2)}

    def I(rr, ss):
        q, k = c_direct(n, rr, ss)
        return sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k
    base = I({}, {})
    ya = {}
    for i in range(1, n + 1):
        ya[i] = sum(I({(i, k): 1}, {}) for k in range(i, n + 1)) * l1 - I({}, {i: 1}) * l3
    yb = {1: l1 * base - ya[1]}
    for i in range(2, n + 1):
        yb[i] = ya[i - 1] - ya[i]
    yb[n + 1] = l3 * base + ya[n]
    order1 = -(sum(ya[i] * za[i] for i in range(1, n + 1)) + sum(yb[i] * zb[i] for i in range(1, n + 2)))
    order1 = sp.expand(order1 * (l1 * l3) ** (n // 2))
    if n == 2:
        claim = -(sp.pi / 30) * l1 * l3 * (l1 * (4 * za[1] + za[2] + 6 * zb[1] + 3 * zb[2] + zb[3]) + l3 * (-za[1] - 4 * za[2] + zb[1] + 3 * zb[2] + 6 * zb[3]))
    else:
        claim = -(sp.pi ** 2 / 6300) * l1 ** 2 * l3 ** 2 * (l1 * (20 * za[1] + 10 * za[2] + 4 * za[3] + za[4] + 15 * zb[1] + 10 * zb[2] + 6 * zb[3] + 3 * zb[4] + zb[5])
                                                           + l3 * (-za[1] - 4 * za[2] - 10 * za[3] - 20 * za[4] + zb[1] + 3 * zb[2] + 6 * zb[3] + 10 * zb[4] + 15 * zb[5]))
    print("  n=%d order 1: -sum_e z_e int g_end y_e = %s;  equals claim: %s" % (n, sp.factor(order1), sp.simplify(order1 - sp.expand(claim)) == 0))
