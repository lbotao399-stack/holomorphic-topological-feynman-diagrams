"""(c), (f) Functional forms and coefficients.
 (1) dictionary: alt((1/3)(p1a^2 p2b^2 - p1a p1b p2a p2b) f) == [p1p2]^2 f,
                 alt((1/3)(p1a^3 p2b^3 - 3 p1a^2 p1b p2a p2b^2) d) == [p1p2]^3 d   (a, b = components 1, 2)
     i.e. (1/3!) w(A,A,A) for w = [p1p2]^2 f, [p1p2]^3 d equals
          (1/3) int Omega f_abc (d1^2 A^a d2^2 A^b - d1d2 A^a d1d2 A^b) A^c,
          (1/3) int Omega d_abc (d1^3 A^a d2^3 A^b - 3 d1^2d2 A^a d1d2^2 A^b) A^c      (d_i = d/dz_i).
 (2) d([p1p2]^3 d)(x_0..x_3) == sum_{i<j} (-1)^{i+j} [p_k p_l]^3 tr([X_i,X_j]{X_k,X_l})  ({k<l} = complement),
     exact evaluation at random integer matrices (gl_3, sl_3) and momenta.
 (3) coefficients (sympy): c_O = pi nu/48, pi nu = -U/(64 pi^4), a^2 = -12 N^2 eps c_O, i_2 = eps pi nu iota,
     iota = -((3N^2+10)/20)[p1p2]^2 f, i_3 = a eps pi nu xi, xi = -((2N^2+5)/30)[p1p2]^3 d, theta = a pi nu theta',
     theta' = -(N^2/60) d([p1p2]^3 d);  values for eps = U."""
import itertools, random
from common import *
import evalx

fW = {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(-1)}}
dW = {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}}
# raw monomials over 3 slots: (p1a, p1b, p2a, p2b, p3a, p3b)
P2 = {(2, 0, 0, 2, 0, 0): Fr(1, 3), (1, 1, 1, 1, 0, 0): Fr(-1, 3)}
P3 = {(3, 0, 0, 3, 0, 0): Fr(1, 3), (2, 1, 1, 2, 0, 0): Fr(-1)}
w2 = ce.alt_raw(ce.raw_from(P2, fW), 3)
w3 = ce.alt_raw(ce.raw_from(P3, dW), 3)
u = ce.bracket_poly(0, 1, 3)
u3 = ce.poly_mul(ce.poly_mul(u, u), u)
D3 = ce.canon_raw(ce.raw_from(u3, dW), 3)
print("(1) alt((1/3)(p1a^2 p2b^2 - p1a p1b p2a p2b) f) == [p1p2]^2 f:", is_zero(ce.lin((1, w2), (-1, f_bracket_sq()))))
print("    alt((1/3)(p1a^3 p2b^3 - 3 p1a^2 p1b p2a p2b^2) d) == [p1p2]^3 d:", is_zero(ce.lin((1, w3), (-1, D3))))

rng = random.Random(5)
for sl in (False, True):
    N = 3
    c = ce.d_ce(ce.restrict_sl(D3) if sl else D3, 3, sl)
    ok, nz = True, False
    for trial in range(6):
        mats = [evalx.rand_mat(N, 4, traceless=sl, rng=rng) for _ in range(4)]
        moms = evalx.rand_moms(4, 4, rng=rng)
        lhs = evalx.eval_cochain(c, mats, moms, N)
        rhs = Fr(0)
        for i, j in itertools.combinations(range(4), 2):
            k, l = [r for r in range(4) if r not in (i, j)]
            br = moms[k][0] * moms[l][1] - moms[k][1] * moms[l][0]
            A = evalx.matmul(mats[i], mats[j])
            Bm = evalx.matmul(mats[j], mats[i])
            comm = [[A[r][s] - Bm[r][s] for s in range(N)] for r in range(N)]
            A2 = evalx.matmul(mats[k], mats[l])
            B2 = evalx.matmul(mats[l], mats[k])
            acom = [[A2[r][s] + B2[r][s] for s in range(N)] for r in range(N)]
            rhs += (-1) ** (i + j) * br ** 3 * evalx.trace(evalx.matmul(comm, acom))
        ok &= (lhs == rhs)
        nz |= lhs != 0
    print("(2) %s_3: d([p1p2]^3 d) == sum_{i<j} (-1)^{i+j}[p_k p_l]^3 tr([X_i,X_j]{X_k,X_l}) at 6 random points: %s (nonzero values: %s)"
          % ("sl" if sl else "gl", ok, nz))

U, eps, Nn, a = sp.symbols('U epsilon N a')
pinu = -U / (64 * sp.pi ** 4)
cO = pinu / 48
a2 = sp.simplify(-12 * Nn ** 2 * eps * cO)
k2 = sp.simplify(eps * pinu * (-(3 * Nn ** 2 + 10) / sp.Integer(20)))
k3 = sp.simplify(a * eps * pinu * (-(2 * Nn ** 2 + 5) / sp.Integer(30)))
kth = sp.simplify(a * pinu * (-Nn ** 2 / sp.Integer(60)))
print("(3) a^2 = %s ; i_2 = (%s) [p1p2]^2 f ; i_3 = (%s) [p1p2]^3 d ; theta = (%s) d([p1p2]^3 d)" % (a2, k2, k3, kth))
print("    a^2 eta_N + eps zeta, eta_N = ((N^2+4)/(2N^2)) u2f, zeta = -(pi nu N^2/40) u2f:  coefficient %s == %s: %s"
      % (sp.simplify(a2 * (Nn ** 2 + 4) / (2 * Nn ** 2) + eps * (-pinu * Nn ** 2 / 40)), k2,
         sp.simplify(a2 * (Nn ** 2 + 4) / (2 * Nn ** 2) + eps * (-pinu * Nn ** 2 / 40) - k2) == 0))
sub = {eps: U}
print("    eps = U (U^2 = 1):  a^2 = %s ; i_2 coefficient %s ; i_3 coefficient %s ; theta coefficient %s"
      % (sp.simplify(a2.subs(sub).subs(U ** 2, 1)), sp.simplify(k2.subs(sub).subs(U ** 2, 1)),
         sp.simplify(k3.subs(sub).subs(U ** 2, 1)), kth))
for Nv in (3, 4, 5):
    print("    N = %d: a^2 = %s, i_2 coefficient %s, i_3 coefficient / a %s"
          % (Nv, sp.nsimplify(a2.subs(sub).subs(U ** 2, 1).subs(Nn, Nv)), sp.simplify(k2.subs(sub).subs(U ** 2, 1).subs(Nn, Nv)),
             sp.simplify((k3 / a).subs(sub).subs(U ** 2, 1).subs(Nn, Nv))))
