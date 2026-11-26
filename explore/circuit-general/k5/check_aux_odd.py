"""Exact evaluation of the isolated odd-level class (q odd, k = (q-3)/2):
   b = (-k,-k,lam,lam,lam,-k,lam,lam,lam,lam,lam,lam,-k,-k,-k),  lam = q-1   (numerator triangles {0,1,2}, {3,4,5},
   harmonic K33).  q=3: class (0,0,2,2,2,0,...) = 1; q=5: 23/69984 (DS from level 6, fix_u5.log).
Frame 5 -> oo, 3 -> 0, 4 -> e; d = 2q; dmu(x) = d^dx / pi^(d/2) |x|^(-2 lam) |x-e|^(-2 lam):
   J = int dmu(x0) dmu(x1) dmu(x2) (|x01|^2 |x02|^2 |x12|^2)^k
x_a = t_a e + y_a, y_a in R^m (m = d-1), w_ij = y_i.y_j, |x_ij|^2 = c_ij - 2 w_ij, c_ij = (t_i-t_j)^2 + r_i^2 + r_j^2.
Angular average (independent isotropic unit vectors u_a = y_a/r_a):
   E[(u0.u1)^a (u0.u2)^b (u1.u2)^c] = sum_j C(a,j) C(b,j) j! (a-j-1)!! (b-j-1)!! E[s^(j+c)] / prod_{i<(a+b)/2} (m+2i)
   E[s^(2l)] = (2l-1)!! / prod_{i<l} (m+2i),  E[s^odd] = 0.
One-point moments M(j,k') = int dmu t^j (r^2)^k', r^2 = x^2 - t^2, t = (x^2 + 1 - |x-e|^2)/2  ->  sums of
   G_d(lam-a, lam-b) = Gamma(d/2-A)Gamma(d/2-B)Gamma(A+B-d/2)/(Gamma(A)Gamma(B)Gamma(d-A-B)), d = 2q - 2 eps, eps -> 0.
usage: python3 check_aux_odd.py Q"""
import sys
from math import comb, factorial

import sympy as sp

q = int(sys.argv[1])
assert q % 2 == 1 and q >= 3
k = (q - 3) // 2
lam = q - 1
dd = 2 * q
m = dd - 1
eps = sp.symbols('eps')
d = dd - 2 * eps
X2, Y2 = sp.symbols('X2 Y2')


def dfact(n):
    return 1 if n <= 0 else n * dfact(n - 2)


def Es(n):
    if n % 2:
        return sp.Integer(0)
    l = n // 2
    den = 1
    for i in range(l):
        den *= (m + 2 * i)
    return sp.Rational(dfact(2 * l - 1), den)


def ang(a, b, c):
    if (a + b) % 2:
        return sp.Integer(0)
    den = 1
    for i in range((a + b) // 2):
        den *= (m + 2 * i)
    tot = sp.Integer(0)
    for j in range(0, min(a, b) + 1):
        if (a - j) % 2 or (b - j) % 2:
            continue
        tot += comb(a, j) * comb(b, j) * factorial(j) * dfact(a - j - 1) * dfact(b - j - 1) * Es(j + c)
    return tot / den


def G(A, B):
    return sp.gamma(d / 2 - A) * sp.gamma(d / 2 - B) * sp.gamma(A + B - d / 2) / (
        sp.gamma(A) * sp.gamma(B) * sp.gamma(d - A - B))


_M = {}


def M(j, kk):
    if (j, kk) in _M:
        return _M[(j, kk)]
    t = (X2 + 1 - Y2) / 2
    poly = sp.Poly(sp.expand(t ** j * (X2 - t ** 2) ** kk), X2, Y2)
    tot = 0
    for (a, b), c in poly.terms():
        tot += c * sp.series(G(lam - a, lam - b), eps, 0, 2).removeO()
    v = sp.nsimplify(sp.limit(sp.simplify(tot), eps, 0))
    _M[(j, kk)] = v
    return v


T = sp.symbols('t0:3')
U = sp.symbols('u0:3')            # r_a^2
W = sp.symbols('w01 w02 w12')
c01 = (T[0] - T[1]) ** 2 + U[0] + U[1]
c02 = (T[0] - T[2]) ** 2 + U[0] + U[2]
c12 = (T[1] - T[2]) ** 2 + U[1] + U[2]
P = sp.expand(((c01 - 2 * W[0]) * (c02 - 2 * W[1]) * (c12 - 2 * W[2])) ** k)
poly = sp.Poly(P, *W)
# angular average: w01^a w02^b w12^c -> r0^(a+b) r1^(a+c) r2^(b+c) ang(a,b,c); exponents of r_a even overall
Pt = sp.Integer(0)
for (a, b, c), coef in poly.terms():
    e0, e1, e2 = a + b, a + c, b + c
    if e0 % 2 or e1 % 2 or e2 % 2:
        continue
    A = ang(a, b, c)
    if A == 0:
        continue
    Pt += coef * A * U[0] ** (e0 // 2) * U[1] ** (e1 // 2) * U[2] ** (e2 // 2)
Pt = sp.Poly(sp.expand(Pt), *T, *U)
J = sp.Integer(0)
for mon, coef in Pt.terms():
    pr = sp.Integer(1)
    for a in range(3):
        pr *= M(mon[a], mon[3 + a])
    J += coef * pr
J = sp.nsimplify(sp.simplify(J))
print('q=%d k=%d: M(0,0) = %s (G_%d(%d,%d)); J = %s' % (q, k, _M.get((0, 0)), dd, lam, lam, J))
