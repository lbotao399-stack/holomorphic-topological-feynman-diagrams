"""(a), (d): monomial densities and closed forms.

 C1  (1,1) endpoint density = 2^{-n/2} Pi(r) (n-1)!! lam_1^{n/2} (P_n u_n v_n lam_{n+1})^{n/2}
     = (n-1)!!/2^{n/2} (lam_1 lam_{n+1})^{n/2} v_n^{n/2-1} prod_{k<n} v_k^{1/2}(1-v_k)^{k/2-1} prod_m u_m^{(n-1-m)/2}
     (symbolic identity in r, n = 2..12 even);  exact integral = (4 pi lam_1 lam_{n+1})^{n/2}/((n+1) n!^2), n = 2..30.
 C2  (2,0) control, same derivation: density (-1)^{n(n+1)/2} det(lam_1,lam_{n+1})^n J(r) prod_i a_i h_{n+1}^n
     equals (-1)^{n(n+1)/2} x S1 eq-1-41 (symbolic in r, n = 1..12); exact integral 1/((n!)^3 (n+1)), n = 1..30.
"""
import sys
import math
from fractions import Fraction as Fr

sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/fans-ht/code")
import sympy as sp
import fanlib as F
from integrate_general import beta2


def rsyms(n):
    u = {k: sp.Symbol("u%d" % k, positive=True) for k in range(1, n + 1)}
    v = {k: sp.Symbol("v%d" % k, positive=True) for k in range(1, n + 1)}
    return u, v


def pieces(n, u, v):
    P = {1: sp.Integer(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * u[k - 1]
    t = F.t_of_r(n, u, v)
    a = {i: t[2 * (i - 1)] for i in range(1, n + 1)}
    B = {n: u[n]}
    R = {}
    for i in range(n - 1, 0, -1):
        R[i] = B[i + 1] / (1 - v[i])
        B[i] = u[i] * R[i]
    J = sp.Integer(1)
    for i in range(1, n):
        J *= v[i] * R[i] / (a[i] * B[i] * B[i + 1])
    J /= a[n] * B[n] * v[n]
    U = sp.Mul(*[u[k + 1] ** k * (1 - v[k]) ** (-k) / v[k] for k in range(1, n)])
    return P, a, J, U


for n in range(1, 13):
    u, v = rsyms(n)
    P, a, J, U = pieces(n, u, v)
    h = P[n] * u[n] * v[n]
    lhs = J * sp.Mul(*[a[i] for i in range(1, n + 1)]) * h ** n
    s1 = v[n] ** (n - 1) * sp.Mul(*[u[j] ** (n - j) * v[j] * (1 - v[j]) ** (j - 1) for j in range(1, n)])
    assert sp.simplify(lhs / s1) == 1, ("C2", n)
    if n % 2 == 0:
        Pi = sp.sqrt(U) * J * sp.Mul(*[a[i] for i in range(1, n + 1)])
        Pi_claim = v[n] ** -1 * sp.Mul(*[sp.sqrt(v[k]) * (1 - v[k]) ** sp.Rational(k - 2, 2) for k in range(1, n)]) * \
            sp.Mul(*[u[k] ** sp.Rational(-(k + 1), 2) for k in range(1, n + 1)])
        for trial in range(3):
            sub = {x: sp.Rational(3 + 5 * trial + k, 41) for k, x in enumerate(list(u.values()) + list(v.values()))}
            val = (Pi / Pi_claim).subs(sub)
            assert sp.nsimplify(sp.simplify(val ** 2)) == 1 and float(val) > 0, ("Pi", n)
        mono = Pi_claim * h ** (n // 2)
        claim = v[n] ** sp.Rational(n - 2, 2) * sp.Mul(*[sp.sqrt(v[k]) * (1 - v[k]) ** sp.Rational(k - 2, 2) for k in range(1, n)]) * \
            sp.Mul(*[u[m] ** sp.Rational(n - 1 - m, 2) for m in range(1, n + 1)])
        for trial in range(3):
            sub = {x: sp.Rational(2 + 3 * trial + k, 53) for k, x in enumerate(list(u.values()) + list(v.values()))}
            val = (mono / claim).subs(sub)
            assert sp.simplify(val ** 2) == 1 and float(val) > 0, ("C1", n)
print("C2: symbolic identity, n = 1..12; C1 and Pi(r) = U^{1/2} J prod a_i: exact at 3 rational points each (square and sign), even n = 2..12")


def dfact(k):
    o = 1
    while k > 1:
        o *= k
        k -= 2
    return o


for n in range(1, 31):
    # (2,0): integral of S1 eq-1-41 monomial
    q = Fr(1)
    for j in range(1, n):
        q *= Fr(1, n - j + 1) * Fr(math.factorial(1) * math.factorial(j - 1), math.factorial(j + 1))
    q *= Fr(math.factorial(0) * math.factorial(n - 1), math.factorial(n + 1))
    assert q == Fr(1, math.factorial(n) ** 3 * (n + 1)), n
    if n % 2 == 0:
        qq, kk = Fr(dfact(n - 1), 2 ** (n // 2)), 0
        for m in range(1, n):
            qq *= Fr(2, n + 1 - m)
        for k in range(1, n):
            dq, dk = beta2(3, k)
            qq *= dq
            kk += dk
        dq, dk = beta2(1, n + 2)               # int_{Delta^2} u^{-1/2} v^{n/2-1} = B(1/2, n/2+1)/(n/2)
        qq *= dq / Fr(n, 2)
        kk += dk
        assert kk == n and qq == Fr(4 ** (n // 2), (n + 1) * math.factorial(n) ** 2), (n, qq, kk)
print("closed forms: (1,1) K_end = (4 pi)^{n/2}/((n+1) n!^2), n = 2..30 even; (2,0) 1/((n!)^3 (n+1)), n = 1..30")
for n in (2, 4, 6, 8, 10):
    print("  n=%2d: (1,1) %s * pi^%d   (2,0) %s" % (n, Fr(4 ** (n // 2), (n + 1) * math.factorial(n) ** 2), n // 2,
                                             Fr(1, math.factorial(n) ** 3 * (n + 1))))
