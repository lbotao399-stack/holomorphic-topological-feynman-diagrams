"""S1: reflection symmetry K(lam_1..lam_{n+1}) = K(lam_{n+1}..lam_1), n = 2, 4, 6, 8 (exact polynomials).
   S2: first elimination: g_n(r; lam) = lam_1 F(r without u_1; lam_2 + u_1 lam_1, lam_3, ...): the density
       with lam_2 -> lam_2 - u_1 lam_1, divided by lam_1, is independent of u_1 and lam_1 (n = 2, 4, 6;
       symbolic in lam, exact at rational r with u_1 varied).  Consequence:
       K(lam) = Psi(lam_1 + lam_2, lam_3, ..) - Psi(lam_2, lam_3, ..),  d Psi/d lam_2 = int g/lam_1 (u_1-independent).
   S3: sources (lam_1, lam_2, lam_{n+1}):  K = K_end(n) lam_{n+1}^{n/2} [(lam_1 + lam_2)^{n/2} - lam_2^{n/2}], n = 2..8."""
import sys, pickle, random
from fractions import Fraction as Fr
sys.path.insert(0, ".")
import sympy as sp
import fanlib as F
data = pickle.load(open("../general_kernels.pkl", "rb"))
for n in (2, 4, 6, 8):
    lam = sp.symbols("l1:%d" % (n + 2))
    K = sp.sympify(data[n], locals={str(x): x for x in lam})
    Kr = K.subs({lam[i]: lam[n - i] for i in range(n + 1)}, simultaneous=True)
    assert sp.expand(K - Kr) == 0
    Kend = sp.Integer(4) ** (n // 2) * sp.pi ** (n // 2) / ((n + 1) * sp.factorial(n) ** 2)
    sub = {x: 0 for i, x in enumerate(lam) if i not in (0, 1, n)}
    assert sp.expand(K.subs(sub) - Kend * lam[n] ** (n // 2) * ((lam[0] + lam[1]) ** (n // 2) - lam[1] ** (n // 2))) == 0
print("S1 reflection symmetry and S3 (lam_1, lam_2, lam_{n+1}) formula hold, n = 2, 4, 6, 8")
rng = random.Random(3)


def LMsym(n, u, v, lam):
    q = lambda x: sp.Rational(x.numerator, x.denominator)
    P = {1: sp.Integer(1)}
    for k in range(2, n + 2):
        P[k] = P[k - 1] * q(u[k - 1])
    R = {n: sp.Integer(1)}
    for k in range(n - 1, 0, -1):
        R[k] = R[k + 1] * q(u[k + 1]) / (1 - q(v[k]))
    c = {k: q(u[k]) * (1 - q(u[k])) * R[k] for k in range(1, n + 1)}
    L, acc = {}, 0
    for i in range(1, n + 1):
        acc += lam[i - 1] / P[i]
        L[i] = acc
    M = {j: sum((P[k] ** 2 * c[k] * (L[k] - L[j]) for k in range(j + 1, n + 1)), sp.Integer(0)) + P[n] * q(u[n]) * q(v[n]) * lam[n]
         for j in range(1, n + 1)}
    return L, M
for n in (2, 4, 6):
    lam = list(sp.symbols("l1:%d" % (n + 2)))
    u, v = F.random_r(n, rng)
    vals = []
    for u1 in (Fr(1, 7), Fr(3, 5), Fr(9, 10)):
        uu = dict(u)
        uu[1] = u1
        L, M = LMsym(n, uu, v, [lam[0], lam[1] - sp.Rational(u1.numerator, u1.denominator) * lam[0]] + lam[2:])
        rat, rad = F.Pi_prefactor(n, uu, v)
        g = sp.expand(sp.Rational(rat.numerator, rat.denominator) * F.hafnian_LM(n, L, M) / 2 ** (n // 2))
        vals.append((g, rad))
    for (g, rad) in vals:
        q = sp.expand(sp.cancel(g / lam[0]))
        assert not q.has(lam[0])
    q0 = sp.expand(sp.cancel(vals[0][0] / lam[0])) * sp.sqrt(vals[0][1] * vals[0][1].denominator ** 2) / vals[0][1].denominator
    for (g, rad) in vals[1:]:
        q = sp.expand(sp.cancel(g / lam[0])) * sp.sqrt(rad * rad.denominator ** 2) / rad.denominator
        assert sp.simplify(q - q0) == 0
print("S2 first-elimination structure holds, n = 2, 4, 6")
