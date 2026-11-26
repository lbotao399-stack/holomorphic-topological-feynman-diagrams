"""Closed forms and specializations of T_1, T_3 (single-x evaluation, p_k = tr x^k or str x^k):
gl(t) from assemble_n3.log (formal ribbon), so(N) from the explicit contractions N = 5..9 (degree <= 3 in N),
sp(2n) by N -> -2n with tr x^4 -> -tr x^4 (checked against sp4, sp6)."""
import sympy as sp
t, N = sp.symbols('t N')
p1, p2, p3, p4 = sp.symbols('p1 p2 p3 p4')
T1 = sp.Rational(1, 8)*p2**2 - sp.Rational(1, 6)*p3*p1 + t/24*p4
T3 = (-sp.Rational(1, 17280)*p1**4 - sp.Rational(1, 5760)*p2**2 + sp.Rational(1, 4320)*p3*p1 + t/8640*p2*p1**2
      - t/17280*p4 - t**2/34560*p2**2 - t**2/25920*p3*p1 + t**3/103680*p4)
print('gl(t): T_1 =', sp.factor(T1 * 48), '/48   (= Str_ad x^4 / 48)')
print('gl(t): T_3 =', sp.collect(sp.expand(T3), [p1, p2, p3, p4]))
print('gl(t) = gl(N|N) (t = 0): T_3 + T_1/720 =', sp.factor(sp.expand(T3 + T1 / 720).subs(t, 0)))
sl = lambda e: sp.expand(e.subs(p1, 0))
print('sl(t): T_1 =', sl(T1), ';  T_3 =', sp.collect(sl(T3), [p2, p4]))
print('sl(t): single-trace parts: T_1 ->', sp.factor(sl(T1).coeff(p4)), ', T_3 ->', sp.factor(sl(T3).coeff(p4)))
for n_ in (2, 3):
    # traceless n x n: p4 = p2^2/2 for n = 2, 3 (Cayley-Hamilton)
    a1 = sl(T1).subs(t, n_).subs(p4, p2**2 / 2); a3 = sl(T3).subs(t, n_).subs(p4, p2**2 / 2)
    print('sl%d: T_1 = %s, T_3 = %s, T_3/T_1 = %s' % (n_, sp.factor(a1), sp.factor(a3), sp.factor(a3 / a1)))
q4, q22 = sp.symbols('trx4 trx2sq')
alpha1, beta1 = (N - 8) / sp.Integer(48), sp.Rational(1, 16)
alpha3 = (N - 8) * (N - 5) * (N + 4) / sp.Integer(829440)
beta3 = -(N**2 - 11 * N + 36) / sp.Integer(276480)
data = {5: (0, sp.Rational(-1, 46080)), 6: (sp.Rational(-1, 41472), sp.Rational(-1, 46080)), 7: (sp.Rational(-11, 414720), sp.Rational(-1, 34560)),
        8: (0, sp.Rational(-1, 23040)), 9: (sp.Rational(13, 207360), sp.Rational(-1, 15360))}
print('so(N) closed forms reproduce N = 5..9:', all(alpha3.subs(N, k) == a and beta3.subs(N, k) == b for k, (a, b) in data.items()))
print('sp(2n) via N -> -2n, trx4 -> -trx4: sp4', (-alpha3.subs(N, -4), beta3.subs(N, -4)), ' sp6', (-alpha3.subs(N, -6), beta3.subs(N, -6)),
      ' [computed: sp4 (0, -1/2880), sp6 (77/207360, -23/46080)]')
print('so(N): T_3^irr / T_1^irr =', sp.factor(alpha3 / alpha1), ';  so8: T_3/T_1 =', beta3.subs(N, 8) / beta1)
a = sp.Symbol('a')
print('gl(1) (abelian) check: T_1(t=1, p_k = a^k) =', sp.expand(T1.subs({t: 1, p1: a, p2: a**2, p3: a**3, p4: a**4})),
      '; T_3 =', sp.expand(T3.subs({t: 1, p1: a, p2: a**2, p3: a**3, p4: a**4})))
print('psl(N|N) (t = 0, p1 = 0): T_1 =', sl(T1).subs(t, 0), ', T_3 =', sl(T3).subs(t, 0))
print('normalized cochains: omega_1 = 24 N_1 Delta T_1 = (3/pi^3) Delta T_1 ; omega_3 = 24 N_3 Delta^3 T_3 = -(3/(64 pi^9)) Delta^3 T_3')
print('gl(N|N) quadruple-trace term of omega_3: ', sp.Rational(-3, 64) * sp.Rational(-1, 17280), '/pi^9 * Delta^3 prod_i str(x_i)')
# so6 = sl4 (Lambda^2 C^4): tr_6 = 2 tr_4, W_Gamma ~ kappa^{1-n}  =>  T_n^{sl4} = 2^{n-1} T_n^{so6}|_{tr_6 x^4 = -4 p4 + 3 p2^2, (tr_6 x^2)^2 = 4 p2^2}
for nn, (al, be) in ((1, (alpha1, beta1)), (3, (alpha3, beta3))):
    so6 = al.subs(N, 6) * (-4 * p4 + 3 * p2**2) + be.subs(N, 6) * 4 * p2**2
    sl4 = sl(T1 if nn == 1 else T3).subs(t, 4)
    print('so6 = sl4 check, n=%d:' % nn, sp.expand(2**(nn - 1) * so6 - sl4) == 0)
