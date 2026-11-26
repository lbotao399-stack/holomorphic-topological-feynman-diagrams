"""Explicit-matrix checks (no formal d / bracket / contraction rules used on the left-hand sides):
 (1) Moyal: (1/2){j_M,j_M} = d i_M on gl2, gl3 at random points, bracket by dual-basis shuffle sum,
     d by the CE sum with explicit commutators.
 (2) sl3: (1/2){D,D}(x) - F_223(x)/(12N^2) - (d eta_N)(x) = 0 at random traceless points, bracket explicit.
 (3) basis change identities on gl_N, N = 3,4,5: D_sl(pi X) = (D - (2/N)S1 + (4/N^2)S0)(X),
     M1 = [p1p2] sum_cyc trX1 tr(pi X2 pi X3) = (S1 - (3/N)S0)(X),  j_M = (1/2)pi^*D_sl + M1/N + U/N^2.
 (4) gl2: tr(X{Y,Z}) = trX tr(YZ) + trY tr(XZ) + trZ tr(XY) - trX trY trZ  (symbolic, sympy)."""
import random
from fractions import Fraction as Fr
import sympy
from vce import *
from vce import _polymul
from vnum import ev, d_explicit, bracket_explicit, basis_gl, basis_sl, inv_gram, rand_point, matmul, tr
from vw import W223

rng = random.Random(99)
pq3 = bracket_pq(3, 0, 1)
D = from_poly_tau(pq3, [(1, (1, 2, 0)), (1, (2, 0, 1))])
S1 = from_poly_tau(pq3, [(1, (0, 2, 1)), (1, (2, 1, 0)), (1, (1, 0, 2))])
S0 = from_poly_tau(pq3, [(1, (0, 1, 2))])
jM = scale(D, Fr(1, 2))
iM = scale(from_poly_tau(_polymul(pq3, pq3), [(1, (1, 2, 0)), (-1, (2, 0, 1))]), Fr(1, 8))

for N in (2, 3):
    B = basis_gl(N); Gi = inv_gram(B)
    ok = True
    for _ in range(3):
        X, p = rand_point(N, 4, False, rng)
        lhs = Fr(1, 2) * bracket_explicit(jM, 3, jM, 3, X, p, N, B, Gi)
        rhs = d_explicit(iM, X, p, N)
        ok &= (lhs == rhs)
        nz = lhs != 0
    print(f"gl{N}: (1/2){{j_M,j_M}} == d i_M at 3 random points (explicit): {ok}; value nonzero at last point: {nz}")

# (2) sl3 coboundary
N = 3
B = basis_sl(N); Gi = inv_gram(B)
Dsl = restrict_sl(D)
P4 = polyprod(bracket_pq(4, 0, 1), bracket_pq(4, 2, 3))
F = {}
for mm, c in P4.items():
    for (sig, n), v in W223(True).items():
        add(F, (sig, mm, n), Fr(c) * v)
F223 = restrict_sl(alt(F, 4))
base = from_poly_tau(_polymul(pq3, pq3), [(1, (1, 2, 0)), (-1, (2, 0, 1))])
eta = scale(base, Fr(N * N + 4, 2 * N * N))
ok = True
for _ in range(3):
    X, p = rand_point(N, 4, True, rng)
    lhs = Fr(1, 2) * bracket_explicit(Dsl, 3, Dsl, 3, X, p, N, B, Gi)
    rhs = ev(F223, X, p, N) / (12 * N * N) + d_explicit(eta, X, p, N)
    ok &= (lhs == rhs)
print(f"sl3: (1/2){{D,D}} == F_223/(12N^2) + d eta_3 at 3 random traceless points (explicit bracket and d): {ok}")


# (3) basis change
def piX(X, N):
    t = tr(X)
    return [[X[i][j] - (Fr(t, N) if i == j else 0) for j in range(N)] for i in range(N)]


for N in (3, 4, 5):
    ok1 = ok2 = ok3 = True
    for _ in range(4):
        X, p = rand_point(N, 3, False, rng)
        Xp = [piX(x, N) for x in X]
        lhs = ev(D, Xp, p, N)
        rhs = ev(D, X, p, N) - Fr(2, N) * ev(S1, X, p, N) + Fr(4, N * N) * ev(S0, X, p, N)
        ok1 &= (lhs == rhs)
        br = p[0][0] * p[1][1] - p[0][1] * p[1][0]
        M1 = br * sum(tr(X[a]) * tr(matmul(Xp[b], Xp[c])) for (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))
        ok2 &= (M1 == ev(S1, X, p, N) - Fr(3, N) * ev(S0, X, p, N))
        ok3 &= (ev(jM, X, p, N) == Fr(1, 2) * lhs + M1 / N + ev(S0, X, p, N) / (N * N))
    print(f"gl{N}: pi^*D_sl = D - (2/N)S1 + (4/N^2)S0: {ok1};  M1 = S1 - (3/N)S0: {ok2};  "
          f"j_M = (1/2)pi^*D_sl + M1/N + U/N^2: {ok3}")

# (4) gl2 identity, symbolic
xs = sympy.symbols('a0:4'); ys = sympy.symbols('b0:4'); zs = sympy.symbols('c0:4')
Xm = sympy.Matrix(2, 2, xs); Ym = sympy.Matrix(2, 2, ys); Zm = sympy.Matrix(2, 2, zs)
lhs = (Xm * (Ym * Zm + Zm * Ym)).trace()
rhs = Xm.trace() * (Ym * Zm).trace() + Ym.trace() * (Xm * Zm).trace() + Zm.trace() * (Xm * Ym).trace() - Xm.trace() * Ym.trace() * Zm.trace()
print(f"gl2 symbolic: tr(X{{Y,Z}}) == trX tr(YZ) + trY tr(XZ) + trZ tr(XY) - trX trY trZ: {sympy.expand(lhs - rhs) == 0}")
