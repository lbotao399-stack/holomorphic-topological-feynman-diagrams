"""Arithmetic of claim d10 from the verified class data:
 [o_2] = c_O [F_223] with [F_133] = -[F_223] (vformal.log);  q(J) = -eps [o_2];
 sl_N: q(aD) = a^2/(12N^2) [F_223];  gl_N: lam(alpha,beta,gamma) = (4alpha^2 - 2N^2 beta^2 + N^3 beta gamma)/(48N^2)."""
import sympy as sp

N, eps, nu, a, b, s, t = sp.symbols('N epsilon nu a b s t')
pi = sp.pi
cO = {'c3 sign': (pi * nu / 240) * (3 + 2), 'K_133 flipped': (pi * nu / 240) * (-3 + 2),
      'K_223 flipped': (pi * nu / 240) * (3 - 2)}
claims = {'c3 sign': -eps * N ** 2 * pi * nu / 4, 'K_133 flipped': eps * N ** 2 * pi * nu / 20,
          'K_223 flipped': -eps * N ** 2 * pi * nu / 20}
for k, c in cO.items():
    sol = sp.solve(sp.Eq(a ** 2 / (12 * N ** 2), -eps * c), a ** 2)[0]
    print(f"{k}: c_O = {sp.simplify(c)};  a^2 = {sp.simplify(sol)};  matches d10: {sp.simplify(sol - claims[k]) == 0}")
sol = sp.simplify(-eps * cO['c3 sign'] * 12 * N ** 2)
for sg in (1, -1):
    print(f"  Nnorm = {'+' if sg > 0 else '-'}hbar^2/(64 pi^5): a^2 = {sp.simplify(sol.subs(nu, sg / (64 * pi ** 5)))}")
al, be, ga = sp.symbols('alpha beta gamma')
lam = (4 * al ** 2 - 2 * N ** 2 * be ** 2 + N ** 3 * be * ga) / (48 * N ** 2)
mu = -N * be * (N * ga - be) / 48
# b(S1 - (2/N)S0): (alpha,beta,gamma) = (0, b, b/N)
print("J = b(S1-(2/N)S0): mu =", sp.simplify(mu.subs({al: 0, be: b, ga: b / N})), " lam =", sp.simplify(lam.subs({al: 0, be: b, ga: b / N})))
# s j_M + t pi^*D_sl: (alpha,beta,gamma) = (s/2 + t, s/N, s/N^2)
sub = {al: s / 2 + t, be: s / N, ga: s / N ** 2}
print("J = s j_M + t pi^*D_sl: mu =", sp.simplify(mu.subs(sub)), " 12N^2 lam =", sp.factor(sp.simplify(12 * N ** 2 * lam.subs(sub))))
print("branch beta = N gamma: lam - (alpha^2/(12N^2) - beta^2/48) =", sp.simplify(lam.subs(ga, be / N) - (al ** 2 / (12 * N ** 2) - be ** 2 / 48)))
print("branch beta = 0:       lam - alpha^2/(12N^2) =", sp.simplify(lam.subs(be, 0) - al ** 2 / (12 * N ** 2)))
