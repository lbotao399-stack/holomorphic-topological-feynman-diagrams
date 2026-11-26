"""V5: assembly arithmetic of O_2 and the counting of (e) (verifier's own code; exact sympy)."""
import sys
import sympy as sp
import networkx as nx
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
import htk
from htk import catalog

print("V5  assembly and counting")
pi, hbar, U = sp.pi, sp.Symbol('hbar'), sp.Symbol('U')
Delta = sp.Symbol('Delta')                  # [d1 d2][d3 d4]
W133, W223, F133, F223 = sp.symbols('W133 W223 F133 F223')
K133, K223 = -pi / 20 * Delta, pi / 30 * Delta
Aut = 4
G_over_K = sp.Rational(-1, 64) / pi**5      # V2 pointwise constant, (2,1) thetas
Nnorm = U * G_over_K * hbar**2              # U = universal sign (open)
O2 = Nnorm * (K133 * W133 + K223 * W223) / Aut
c = sp.simplify(O2.coeff(W223) / Delta / 2)
print("  O_2 = Nnorm/4 (K133 W133 + K223 W223) =", sp.factor(O2))
print("  O_2 / (Delta (-3 W133 + 2 W223)) =", sp.simplify(O2 / (Delta * (-3 * W133 + 2 * W223))),
      "; pi Nnorm/240 =", sp.simplify(pi * Nnorm / 240))
cls = sp.expand((O2 / Delta).subs({W133: F133, W223: F223}).subs(F133, -F223))
print("  [O_2] with [F133] = -[F223]:", sp.factor(cls), "; pi Nnorm/48 [F223] =", sp.factor(pi * Nnorm / 48 * F223),
      "; equal:", sp.simplify(cls - pi * Nnorm / 48 * F223) == 0)
alt = sp.expand((Nnorm * (K133 * F133 - K223 * F223) / Aut / Delta).subs(F133, -F223))
alt2 = sp.expand((Nnorm * (-K133 * F133 + K223 * F223) / Aut / Delta).subs(F133, -F223))
print("  opposite relative sign: K223 flipped ->", sp.factor(alt), "; K133 flipped ->", sp.factor(alt2),
      "(= -+ pi Nnorm/240 [F223])")
print("  |c| = hbar^2/%s pi^4 ; |c'| = hbar^2/%s pi^4" % (sp.simplify(1 / (pi / (240 * 64 * pi**5))) / pi**4,
                                                        sp.simplify(1 / (pi / (48 * 64 * pi**5))) / pi**4))

print("(e) counting: tight |I| = Dn+1 (Schwinger degree), position degree (|V|-1) D = (D-1)|I| + 1, cubic: E = 3|V| - 2|I|")
n, D = sp.symbols('n D')
I_ = D * n + 1
V_ = sp.solve(sp.Eq((sp.Symbol('V') - 1) * D, (D - 1) * I_ + 1), sp.Symbol('V'))[0]
E_ = sp.expand(3 * V_ - 2 * I_)
print("  |V| =", sp.expand(V_), "; E =", sp.factor(E_), "; gh = E - D =", sp.factor(E_ - D))
for Dv in (2, 3):
    for nv in range(1, 7):
        Iv, Vv = Dv * nv + 1, (Dv - 1) * nv + 2
        Ev = 3 * Vv - 2 * Iv
        print("  D=%d n=%d: |I|=%2d |V|=%2d E=%2d gh=%2d%s" % (Dv, nv, Iv, Vv, Ev, Ev - Dv,
              "  (n odd: rho_n = 0 at T = 1)" if nv % 2 else ""))
for Dv, nv in [(2, 1), (2, 2), (3, 1), (3, 2)]:
    gs = catalog.enumerate_tight(Dv, nv)
    rows = []
    for G in gs:
        degs = sorted(d for _, d in G.degree())
        if max(degs) <= 3:
            rows.append((degs, sum(3 - d for d in degs)))
    print("  catalog D=%d n=%d: %d tight graphs; max degree <= 3: %s" % (Dv, nv, len(gs), rows))
