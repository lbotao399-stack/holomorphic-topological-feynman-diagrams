"""Symbolic comparison of the reduction of q3lib with S2 eqs. (1.6)-(1.9), (6.4)-(6.5), (7.10)-(7.14)."""
import sympy as sp
from q3lib import Reduction, LAM, Z, R, integrate_domain

r1, r2, r3, r4, r5, r6 = R[:6]
L = lambda v: sp.Matrix(LAM[v])
br = lambda i, j, k: sp.Matrix.hstack(L(i), L(j), L(k)).det()
z = [sp.Matrix(Z[j]) for j in range(3)]


def kernel_density(name, zero):
    red = Reduction(name)
    sgn = {'box': 1, 'A': -1, 'B': 1}[name]          # from check_kernels.py (brute force)
    f = sp.prod(red.factors())
    sub = {LAM[v][a]: 0 for v in zero for a in range(3)}
    return red, sp.expand(sgn * f.subs(sub)), sub


# box
red, J, sub = kernel_density('box', [])
print('box: I density - det(l1,l2,l3) =', sp.simplify(J - br(0, 1, 2)))
X, _ = red.displacement()
ok = all(sp.expand(X[v][a] - (red.steps[0]['r'][v] * Z[0][a] + (1 - red.steps[0]['r'][v]) * Z[1][a])) == 0
         for v in range(3) for a in range(3))
print('box: X_j = r_j z_1 + (1 - r_j) z_2:', ok)
print('box: int =', sp.factor(integrate_domain(J, red.domain())))

# (A): vertices 1..5 -> 0..4, lambda_3 = 0
red, J, sub = kernel_density('A', [2])
JA = r1 ** 2 * (r5 * br(0, 1, 3) + r6 * br(0, 1, 4)) * ((1 - r2) * br(0, 3, 4) + (1 - r3) * br(1, 3, 4))
print('(A): I density - J_A =', sp.expand(J - JA))
X, _ = red.displacement()
W = (1 - r1) * z[0] + r1 * z[1] - z[2]
XA = {0: r2 * z[0] + (1 - r2) * (z[2] + r4 * W), 1: r3 * z[0] + (1 - r3) * (z[2] + r4 * W),
      3: z[2] + r5 * W, 4: z[2] + r6 * W}
print('(A): X_v of (1.8):', all(sp.expand(X[v][a] - XA[v][a]) == 0 for v in XA for a in range(3)))
val = integrate_domain(J, red.domain())
print('(A): int J_A - (2<124>+<125>)(<145>+2<245>)/432 =',
      sp.expand(val - (2 * br(0, 1, 3) + br(0, 1, 4)) * (br(0, 3, 4) + 2 * br(1, 3, 4)) / 432))

# (B): u, w, q, a, b -> 0..4, lambda_q = 0
red, J, sub = kernel_density('B', [2])
JB = r1 * (1 - r1) * (r5 * br(0, 1, 3) + r6 * br(0, 1, 4)) * ((1 - r2) * br(0, 3, 4) + (1 - r3) * br(1, 3, 4))
print('(B): I density + J_B =', sp.expand(J + JB))
X, _ = red.displacement()
XB = {0: r2 * z[0] + (1 - r2) * (z[2] + r4 * W), 1: r3 * z[1] + (1 - r3) * (z[2] + r4 * W),
      3: z[2] + r5 * W, 4: z[2] + r6 * W}
print('(B): X_v of (1.8):', all(sp.expand(X[v][a] - XB[v][a]) == 0 for v in XB for a in range(3)))
val = integrate_domain(J, red.domain())
print('(B): int + (2<uwa>+<uwb>)(<uab>+<wab>)/288 =',
      sp.expand(val + (2 * br(0, 1, 3) + br(0, 1, 4)) * (br(0, 3, 4) + br(1, 3, 4)) / 288))
