"""C5 (pentagon kernel), pentagon part of C6.  Fresh mesh-current solver (common.py), no htk, no explorer code.

(a) exact, random rational t (not normalized): D(t) (sum t)^5 == <0123>; y_{e0} = -sum_j r_j lam_j,
    r_j = (sum of t on the edges after vertex j)/T; -sum_e y_e.z_e == sum_j lam_j.(r_j z_1 + (1-r_j) z_2),
    z_{e0} = z_1, z_{e4} = -z_2.
(b) exact: Taylor coefficients of K(lam;z) through order 3 in z from the t-simplex (Dirichlet moments in t,
    y linear in t on sum t = 1) == <0123> * ordered-r moments (moments by iterated sympy integration).
(c) float: direct quadrature of D(t) exp(-sum y.z) over the 4-simplex in t (fresh float solver) vs
    <0123> e^{b} exp[0,s_1,..,s_4] (divided difference, mpmath 40 digits), s_k = sum_{j<=k} lam_j.(z_1 - z_2),
    b = sum_j lam_j.z_2.
(d) claimed moment formula int r1^a r2^b r3^c r4^d = 1/((d+1)(c+d+2)(b+c+d+3)(a+b+c+d+4)) vs sympy, all a+b+c+d<=6.
"""
import sys, random, itertools, math
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
import numpy as np
import sympy as sp
import mpmath as mp
from common import Network, pentagon, bracket, gauss_simplex, H

nv, edges = pentagon()
net = Network(nv, edges)
print('pentagon edges', edges, ' BFS tree', net.tree, ' cotree', net.cotree)
rng = random.Random(2024)
outs = [0, 1, 2, 3]

# (a)
ok_a = True
for it in range(6):
    t = [Fr(rng.randint(1, 50), rng.randint(1, 7)) for _ in range(5)]
    lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in outs}
    D, y = net.density(t, lam)
    assert net.check_kcl(y, lam)
    T = sum(t)
    b0123 = bracket(lam, (0, 1, 2, 3))
    c1 = (D * T ** 5 == b0123)
    r = [sum(t[v + 1:]) / T for v in outs]          # vertex v sits between e_v and e_{v+1}
    c2 = all(y[0][a] == -sum(r[v] * lam[v][a] for v in outs) for a in range(H))
    z1 = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(H)]
    z2 = [Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(H)]
    lhs = -sum(y[0][a] * z1[a] for a in range(H)) - sum(y[4][a] * (-z2[a]) for a in range(H))
    rhs = sum(lam[v][a] * (r[v] * z1[a] + (1 - r[v]) * z2[a]) for v in outs for a in range(H))
    c3 = lhs == rhs
    print('(a) point %d: D(t) T^5 = %s, <0123> = %s: %s; y_e0 = -sum r_j lam_j: %s; -sum y.z = sum lam.x: %s'
          % (it, D * T ** 5, b0123, c1, c2, c3))
    ok_a &= c1 and c2 and c3
print('(a) all:', ok_a)

# (d) moment formula
r1, r2, r3, r4 = sp.symbols('r1:5', positive=True)


def mom_sympy(a, b, c, d):
    e = r1 ** a * r2 ** b * r3 ** c * r4 ** d
    e = sp.integrate(e, (r4, 0, r3))
    e = sp.integrate(e, (r3, 0, r2))
    e = sp.integrate(e, (r2, 0, r1))
    e = sp.integrate(e, (r1, 0, 1))
    return sp.Rational(e)


ok_d = True
cnt = 0
for tot in range(7):
    for (a, b, c, d) in itertools.product(range(tot + 1), repeat=4):
        if a + b + c + d != tot:
            continue
        claim = sp.Rational(1, (d + 1) * (c + d + 2) * (b + c + d + 3) * (a + b + c + d + 4))
        ok_d &= (mom_sympy(a, b, c, d) == claim)
        cnt += 1
print('(d) moment formula == iterated sympy integral for %d exponent tuples (a+b+c+d <= 6): %s' % (cnt, ok_d))

# (b) Taylor coefficients through order 3, exact.
# On sum t = 1 (t4 = 1 - t0 - t1 - t2 - t3) y_e0 = -sum_e t_e Q_e, Q_e = lam_0 + .. + lam_{e-1}; y_e4 = y_e0 + sum lam.
# K_k(lam; z) := (1/k!) int_simplex D (-y_e0.z1 + y_e4.z2)^k dt ;  D = <0123> on the simplex (from (a)).
ts = sp.symbols('t0:5')
lam_s = {v: [sp.Symbol('L%d_%d' % (v, a)) for a in range(H)] for v in outs}
lamn = {v: [rng.randint(-6, 6) for _ in range(H)] for v in outs}
z1n = [sp.Rational(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(H)]
z2n = [sp.Rational(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(H)]
Q = [[sum(lamn[v][a] for v in range(e)) for a in range(H)] for e in range(5)]
ye0 = [-sum(ts[e] * Q[e][a] for e in range(5)) for a in range(H)]
ye4 = [ye0[a] + sum(lamn[v][a] for v in outs) for a in range(H)]
expo = -sum(ye0[a] * z1n[a] for a in range(H)) + sum(ye4[a] * z2n[a] for a in range(H))


def dirichlet_t(poly):
    """int over {t >= 0, sum t = 1} (Lebesgue dt_0..dt_3) of a polynomial in t0..t4: prod a_i! / (4 + sum a)!"""
    P = sp.Poly(sp.expand(poly), *ts)
    tot = sp.Integer(0)
    for mon, c in P.terms():
        tot += c * sp.Mul(*[sp.factorial(x) for x in mon]) / sp.factorial(4 + sum(mon))
    return tot


b0123 = bracket(lamn, (0, 1, 2, 3))
xs = [[r * z1n[a] + (1 - r) * z2n[a] for a in range(H)] for r in (r1, r2, r3, r4)]
ex_r = sum(lamn[v][a] * xs[v][a] for v in outs for a in range(H))
ok_b = True
for k in range(4):
    lhs = b0123 * dirichlet_t(expo ** k) / sp.factorial(k)
    integrand = sp.expand(ex_r ** k / sp.factorial(k))
    e = sp.integrate(integrand, (r4, 0, r3))
    e = sp.integrate(e, (r3, 0, r2))
    e = sp.integrate(e, (r2, 0, r1))
    e = sp.integrate(e, (r1, 0, 1))
    rhs = b0123 * e
    print('(b) order %d in z: t-simplex = %s, <0123> int_r (sum lam.x)^k/k! = %s, equal: %s' % (k, lhs, sp.nsimplify(rhs), sp.simplify(lhs - rhs) == 0))
    ok_b &= sp.simplify(lhs - rhs) == 0
print('(b) all orders 0..3:', ok_b)

# (c) float quadrature in t vs divided difference
mp.mp.dps = 40
B, W = gauss_simplex(4, 14)
Tm = np.concatenate([B, 1 - B.sum(axis=1, keepdims=True)], axis=1)
for trial in range(3):
    lamf = {v: [rng.randint(-4, 4) for _ in range(H)] for v in outs}
    z1 = np.array([rng.uniform(-0.4, 0.4) for _ in range(H)])
    z2 = np.array([rng.uniform(-0.4, 0.4) for _ in range(H)])
    D, y = net.density_float(Tm, lamf)
    ex = -(y[:, 0, :] @ z1) - (y[:, 4, :] @ (-z2))
    val = np.sum(W * D * np.exp(ex))
    lam_arr = [mp.mpf(0)] * 0
    a = [sum(mp.mpf(lamf[v][i]) * (mp.mpf(z1[i]) - mp.mpf(z2[i])) for i in range(H)) for v in outs]
    bb = sum(mp.mpf(lamf[v][i]) * mp.mpf(z2[i]) for v in outs for i in range(H))
    s = [mp.mpf(0)] + [sum(a[:k + 1]) for k in range(4)]
    dd = mp.mpf(0)
    for i in range(5):
        den = mp.mpf(1)
        for j in range(5):
            if j != i:
                den *= (s[i] - s[j])
        dd += mp.e ** s[i] / den
    exact = bracket(lamf, (0, 1, 2, 3)) * mp.e ** bb * dd
    rel = abs((mp.mpf(val) - exact) / exact)
    print('(c) trial %d: quadrature (14^4 Gauss nodes, fresh float solver) = %.15g, <0123> e^b exp[0,s1..s4] = %s, rel diff %.1e (%d digits)'
          % (trial, val, mp.nstr(exact, 16), float(rel), int(-mp.log10(rel)) if rel > 0 else 16))
