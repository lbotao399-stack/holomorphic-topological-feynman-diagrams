"""dim C^4_inv (G x SL2-invariant alternating 4-cochains of weight (2,2)) from characters:
P = SL2-invariant quartics of the 4 momenta = S^2(B), B = span of brackets = V(2,1,1) (character of (01) on B: -1),
so P = triv + V(3,1) + V(2,2) and P (x) sgn = sgn + V(2,1,1) + V(2,2);
dim C^4_inv = <chi_P chi_T, sgn> = m_T(sgn) + m_T(V(2,1,1)) + m_T(V(2,2)),  T = (g^{(x)4})^g with S4 permuting slots,
chi_T(sigma) = sum_w sgn(w) [coefficient at rho - w rho of prod_{cycles c of sigma} psi^{|c|} chi_adj] (Weyl).
The decomposition of P is checked: chi_P computed directly from the 6 quadratic monomials in u, v, w."""
import itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from t_inv_weyl import CARTAN, NAMES, positive_roots, weyl_rho_points, mul, scale

classes = [('e', (0, 1, 2, 3), 1), ('(01)', (1, 0, 2, 3), 6), ('(01)(23)', (1, 0, 3, 2), 3), ('(012)', (1, 2, 0, 3), 8),
           ('(0123)', (1, 2, 3, 0), 6)]
table = {'triv': [1, 1, 1, 1, 1], 'sgn': [1, -1, 1, 1, -1], 'V31': [3, 1, -1, 0, -1], 'V211': [3, -1, -1, 0, 1],
         'V22': [2, 0, 2, -1, 0]}

# chi_P directly: action of sigma on quadratics in brackets [p_a p_b] (p3 = -(p0+p1+p2))
p = [sp.symbols('x%d y%d' % (i, i)) for i in range(3)]
p.append((-p[0][0] - p[1][0] - p[2][0], -p[0][1] - p[1][1] - p[2][1]))
br = lambda a, b: sp.expand(p[a][0] * p[b][1] - p[a][1] * p[b][0])
u, v, w = br(0, 1), br(0, 2), br(1, 2)
basis = [u * u, v * v, w * w, u * v, u * w, v * w]
gens = [x for pp in p[:3] for x in pp]
polys = [sp.Poly(sp.expand(b), *gens) for b in basis]
keys = sorted({m for P in polys for m in P.as_dict()})
Bm = sp.Matrix([[P.as_dict().get(k, 0) for P in polys] for k in keys])
chiP = []
for _, s, _ in classes:
    sub = {}
    for i in range(3):
        for c in range(2):
            sub[p[i][c]] = p[s[i]][c]
    imgs = [sp.Poly(sp.expand(b.subs(sub, simultaneous=True)), *gens) for b in basis]
    Im = sp.Matrix([[P.as_dict().get(k, 0) for P in imgs] for k in keys])
    X = (Bm.T * Bm).inv() * Bm.T * Im          # matrix of sigma in the basis
    assert Bm * X == Im
    chiP.append(int(X.trace()))
mP = {r: Fr(sum(c[2] * chiP[i] * table[r][i] for i, c in enumerate(classes)), 24) for r in table}
print("chi_P on classes", [c[0] for c in classes], "=", chiP, "; multiplicities:", {k: str(v) for k, v in mP.items()})

for t, A in CARTAN.items():
    r = len(A)
    pos = positive_roots(A)
    chi = {}
    for a in pos:
        chi[a] = chi.get(a, 0) + 1
        chi[tuple(-x for x in a)] = chi.get(tuple(-x for x in a), 0) + 1
    chi[(0,) * r] = chi.get((0,) * r, 0) + r
    pts = weyl_rho_points(A, pos)
    x2 = mul(chi, chi)
    prods = {'e': mul(x2, x2), '(01)': mul(scale(chi, 2), x2), '(01)(23)': mul(scale(chi, 2), scale(chi, 2)),
             '(012)': mul(scale(chi, 3), chi), '(0123)': scale(chi, 4)}
    chiT = [sum(sg * prods[c[0]].get(pt, 0) for pt, sg in pts) for c in classes]
    mT = {k: Fr(sum(c[2] * chiT[i] * table[k][i] for i, c in enumerate(classes)), 24) for k in table}
    dimC = sum(mT[k] * mP_k for k, mP_k in [('sgn', mP['triv']), ('V211', mP['V31']), ('V22', mP['V22'])])
    # <chi_P chi_T, sgn> directly
    direct = Fr(sum(c[2] * chiP[i] * chiT[i] * table['sgn'][i] for i, c in enumerate(classes)), 24)
    print("%-4s %-10s chi_T = %s ; m_T = %s ; dim C^4_inv = %s (direct inner product %s)"
          % (t, NAMES[t], chiT, {k: int(v) for k, v in mT.items()}, dimC, direct))
