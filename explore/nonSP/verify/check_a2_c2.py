"""Checks of the general formulas quoted in a2 and of the weight statement in c2 (mpmath, 40 digits)."""
import mpmath as mp
mp.mp.dps = 40

# (1) Feynman/Dirichlet formula on Delta^2 (measure dgA dgB, gC = 1 - gA - gB), |beta| = 7
def lhs(beta, rho):
    f = lambda a, b: (a ** beta[0] * b ** beta[1] * (1 - a - b) ** beta[2]) / (rho[0] * a + rho[1] * b + rho[2] * (1 - a - b)) ** 10
    return mp.quad(lambda a: mp.quad(lambda b: f(a, b), [0, 1 - a]), [0, 1])
for beta, rho in [((3, 2, 2), (mp.mpf(1) / 3, mp.mpf(3) / 4, mp.mpf(1) / 2)), ((7, 0, 0), (mp.mpf(2) / 5, mp.mpf(1) / 7, mp.mpf(9) / 10)),
                  ((1, 4, 2), (mp.mpf(5) / 6, mp.mpf(1) / 2, mp.mpf(1) / 9))]:
    rhs = mp.factorial(beta[0]) * mp.factorial(beta[1]) * mp.factorial(beta[2]) / mp.factorial(9) * \
        rho[0] ** (-beta[0] - 1) * rho[1] ** (-beta[1] - 1) * rho[2] ** (-beta[2] - 1)
    print("Feynman formula beta=%s: lhs/rhs - 1 = %s" % (beta, mp.nstr(lhs(beta, rho) / rhs - 1, 3)))

# (2) 1/(Q^10 G^b) = Gamma(10+b)/(Gamma(10)Gamma(b)) int_0^1 xi^9 (1-xi)^(b-1) (sum_v (1 - xi r_v) g_v)^(-10-b)
g = [mp.mpf(2) / 7, mp.mpf(3) / 7, mp.mpf(2) / 7]; r = [mp.mpf(1) / 3, mp.mpf(5) / 6, mp.mpf(1) / 8]
Q = sum((1 - rv) * gv for rv, gv in zip(r, g)); Gs = sum(g)
for b in (1, 2, 3, 4):
    val = mp.gamma(10 + b) / (mp.gamma(10) * mp.gamma(b)) * mp.quad(
        lambda x: x ** 9 * (1 - x) ** (b - 1) * sum((1 - x * rv) * gv for rv, gv in zip(r, g)) ** (-10 - b), [0, 1])
    print("Feynman parametrization b=%d: rel dev %s" % (b, mp.nstr(val * Q ** 10 * Gs ** b - 1, 3)))

# (3) c2: int_0^1 xi^(-m) log^k(1-xi) dxi, 1 <= m <= k <= 3, in Q + Q z2 + Q z3 + Q z4
basis = [mp.mpf(1), mp.zeta(2), mp.zeta(3), mp.zeta(4)]
for k in (1, 2, 3):
    for m in range(1, k + 1):
        v = mp.quad(lambda x: x ** (-m) * mp.log(1 - x) ** k, [0, mp.mpf(1) / 2, 1])
        rel = mp.pslq([v] + basis, maxcoeff=10 ** 6, maxsteps=10 ** 6)
        print("int xi^-%d log^%d(1-xi) = %s   relation [v,1,z2,z3,z4]: %s" % (m, k, mp.nstr(v, 20), rel))
