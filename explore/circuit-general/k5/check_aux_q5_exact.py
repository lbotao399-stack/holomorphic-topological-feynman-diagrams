"""Exact independent evaluation of the isolated level-5 class
   b = (-1,-1,4,4,4,-1,4,4,4,4,4,4,-1,-1,-1)  (numerator triangles {0,1,2}, {3,4,5}; harmonic K33, beta = 4)
Frame 5 -> oo, 3 -> 0, 4 -> e, integrated x_0, x_1, x_2 in R^10, dmu(x) = d^10 x / pi^5 |x|^-8 |x-e|^-8:
   J = int dmu^3 |x_01|^2 |x_02|^2 |x_12|^2
x_a = t_a e + y_a, y_a in R^9; direction average (independent, m = 9):  E[w_ij] = 0,  E[w01 w02] = 0,
   E[w01 w02 w12] = r0^2 r1^2 r2^2 / m^2   (w_ij = y_i . y_j)
   P~ = c01 c02 c12 - 8 r0^2 r1^2 r2^2 / m^2,   c_ij = (t_i - t_j)^2 + r_i^2 + r_j^2
One-point moments M(j,k) = int dmu t^j (r^2)^k,  r^2 = x^2 - t^2,  t = (x^2 + 1 - |x-e|^2)/2:
   polynomial in x^2, |x-e|^2  ->  sum of G_d(4-a, 4-b) = Gamma(d/2-A)Gamma(d/2-B)Gamma(A+B-d/2)/(Gamma(A)Gamma(B)Gamma(d-A-B))
   at d = 10 - 2 eps, eps -> 0 (each M(j,k) with j + 2k <= 4 is convergent at d = 10)."""
import sympy as sp

eps = sp.symbols('eps')
X2, Y2 = sp.symbols('X2 Y2')            # x^2, |x - e|^2
d = 10 - 2 * eps


def G(A, B):
    return sp.gamma(d / 2 - A) * sp.gamma(d / 2 - B) * sp.gamma(A + B - d / 2) / (
        sp.gamma(A) * sp.gamma(B) * sp.gamma(d - A - B))


_Gc = {}


def Gval(A, B):
    """finite limit eps -> 0 is not needed termwise; return Laurent series to O(eps)"""
    key = (A, B)
    if key not in _Gc:
        if A <= 0 or B <= 0:
            # 1/Gamma(A) = 0 for A non-positive integer at fixed A: G vanishes unless Gamma(A+B-d/2) or
            # Gamma(d/2 - A) ... are singular; keep full eps dependence by shifting A -> A + 0 (A is fixed):
            pass
        _Gc[key] = sp.series(G(A, B), eps, 0, 2).removeO()
    return _Gc[key]


_M = {}


def M(j, k):
    if (j, k) in _M:
        return _M[(j, k)]
    t = (X2 + 1 - Y2) / 2
    poly = sp.Poly(sp.expand(t ** j * (X2 - t ** 2) ** k), X2, Y2)
    tot = 0
    for (a, b), c in poly.terms():
        A, B = 4 - a, 4 - b
        if A <= 0 or B <= 0:
            # integer non-positive power: 1/Gamma(A) = 0 exactly; G(A,B) at fixed integer A <= 0 vanishes
            # identically in eps (scaleless-type polynomial insertion has no pole in Gamma(d/2 - A)).
            g = sp.simplify(G(A, B))
            tot += c * sp.series(g, eps, 0, 2).removeO()
            continue
        tot += c * Gval(A, B)
    v = sp.simplify(sp.limit(tot, eps, 0))
    _M[(j, k)] = v
    return v


T = sp.symbols('t0:3')
U = sp.symbols('u0:3')
m = 9
c01 = (T[0] - T[1]) ** 2 + U[0] + U[1]
c02 = (T[0] - T[2]) ** 2 + U[0] + U[2]
c12 = (T[1] - T[2]) ** 2 + U[1] + U[2]
Pt = sp.Poly(sp.expand(c01 * c02 * c12 - sp.Rational(8, m * m) * U[0] * U[1] * U[2]), *T, *U)
J = 0
for mon, coef in Pt.terms():
    pr = 1
    for a in range(3):
        pr *= M(mon[a], mon[3 + a])
    J += coef * pr
J = sp.nsimplify(sp.simplify(J))
print('moments M(j,k):', {k: v for k, v in sorted(_M.items())})
print('M(0,0) = %s (G_10(4,4) = 1/18)' % _M[(0, 0)])
print('J =', J, '   (DS from level 6: 23/69984)   equal:', J == sp.Rational(23, 69984))
