"""gl_N, all N >= 4: exactness of omega(N) = omega_0 + N omega_1 + N^2 omega_2 (formal trace basis).
A = d(C^3) has rank r = 6 (N-independent).  With r independent columns J and rows I of A with
det A[I, J] != 0, omega(N) is in col(A) iff every bordered minor det [A[I+i, J] | omega(N)[I+i]]
vanishes (i over the remaining rows).  These minors are polynomials in N; their gcd is printed.
Also for F_223 = alt([p1p2][p3p4] W_223) alone (the class of omega is (pi/48)[F_223]/(pi/60)... see report)."""
import sys
import numpy as np
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import cochain as cc, formal_gl as fg
from c_obstruction import kernel_poly

mom = cc.Momenta(6, 2)
P = mom.p
Nsym = sp.Symbol('N')
monos3 = cc.all_monos(mom, 3, (2, 2))
Wf, Kp = {}, {}
for pqr in [(1, 3, 3), (2, 2, 3)]:
    G, ext, K = kernel_poly(pqr)
    Kp[pqr] = K
    Wf[pqr] = fg.colour_weight_formal(G.nv, G.edges, ext)
om = {pqr: fg.alternate_formal(mom, Kp[pqr], Wf[pqr], 4) for pqr in Wf}
omega = {}
for pqr in om:
    for k_, d in om[pqr].items():
        fg.add(omega, k_, d)
omega = fg.clean(omega)
cols3 = []
for tau in fg.all_taus(3):
    for mu in monos3:
        eta = fg.alternate_formal(mom, mom.mono_expr(mu, 3, [P[0], P[1]]), {tau: {0: 1}}, 3)
        cols3.append(fg.at_N(fg.d_formal(mom, eta, 3), 0) if eta else {})
keys = sorted({k for c in cols3 for k in c} | set(omega) | set(om[(2, 2, 3)]), key=str)
A = sp.Matrix([[c.get(k, 0) for c in cols3] for k in keys])
r = A.rank()
print("rank d(C^3) =", r)
# pivot columns and rows
_, piv_cols = A.rref()
J = list(piv_cols)
_, piv_rows = A[:, J].T.rref()
I = list(piv_rows)
print("pivot columns", J, "pivot rows", I, "det A[I,J] =", A[I, J].det())


def poly_vec(c):
    return sp.Matrix([sum(v * Nsym ** p for p, v in c.get(k, {}).items()) for k in keys])


for name, vec in (("omega", poly_vec(omega)), ("F_223 = omega_223/2", poly_vec(om[(2, 2, 3)]) / 2),
                  ("omega_133 - (3/2) omega_223", poly_vec(om[(1, 3, 3)]) - sp.Rational(3, 2) * poly_vec(om[(2, 2, 3)]))):
    g = sp.Integer(0)
    nonzero = 0
    for i in range(len(keys)):
        if i in I:
            continue
        rows = I + [i]
        M = A[rows, J].row_join(vec[rows, 0])
        m = sp.expand(M.det())
        if m != 0:
            nonzero += 1
            g = sp.gcd(g, m)
    print("%s: %d nonzero bordered minors; gcd over all bordered minors = %s" % (name, nonzero, sp.factor(g)))
