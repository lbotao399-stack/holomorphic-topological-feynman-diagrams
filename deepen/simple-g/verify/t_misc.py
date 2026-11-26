"""(1) Killing / trace ratios r_V for the matrix constructions (K = r_V tr_V), and lambda_V = r_V^2 / 2 from
    F_133 + F_223 = (1/2) beta_K and the rescaling rule.
(2) Rescaling kappa -> c kappa (c = 2, 3): W -> W / c, beta -> c beta, exactly on full tensors (so5, G2, sl3).
(3) sl2: dim C^4_inv: rank of alt(mu * kappa_{pi}) over the 21 bracket products mu and the 3 pairings pi,
    together with beta (6 configurations, full tensors).
(4) O_2 assembly with F_133 = -F_223 + beta/2 (sympy)."""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
import valg, vcoch as vc

# (1)
print("(1) Killing = r_V tr_V:")
for name, N in [('so5', 5), ('so6', 6), ('so7', 7), ('so8', 8), ('sp4', 4), ('sp6', 6), ('G2', 8), ('sl3', 3), ('sl4', 4), ('sl5', 5)]:
    A = valg.build(name)
    tr = np.array([[int(np.trace(a @ b)) for b in A.mats] for a in A.mats], dtype=object)
    K = A.K.astype(object)
    i, j = np.argwhere(tr != 0)[0]
    r = Fr(int(K[i, j]), int(tr[i, j]))
    ok = all(K[a, b] * r.denominator == tr[a, b] * r.numerator for a in range(A.n) for b in range(A.n))
    note = " (tr on O = R + Im O, D(1) = 0, so tr_8 = tr_7)" if name == 'G2' else ""
    print("   %-4s K = %s tr_%d on all entries: %s ; lambda_tr = r^2/2 = %s%s" % (name, r, N if name != 'G2' else 7, ok, r * r / 2, note))

# (2)
print("(2) rescaling:")
for name in ('so5', 'G2', 'sl3'):
    A = valg.build(name)
    S = vc.Setup(A)
    W0 = {G: Fr(1, S.s ** 7) for G in ('133', '223')}
    Wt = {G: S.W_full(G) for G in ('133', '223')}
    rng = np.random.default_rng(1)
    P = vc.config4(rng)
    B0 = vc.beta_tensor(S.TT, P)
    res = []
    for c in (2, 3):
        S2 = vc.Setup(A)
        # rebuild with kappa' = c kappa
        S2.kap = S.kap * c
        S2.f = np.einsum('abd,dc->abc', A.c, S2.kap)
        inv = sp.Matrix(S2.kap.tolist()).inv()
        s2 = int(sp.ilcm(*[sp.Rational(x).q for x in inv]))
        S2.s = s2
        S2.Gi = np.array((inv * s2).tolist(), dtype=np.int64)
        S2.h = vc.ein('yk,aku->ayu', S2.Gi, S2.f)
        t = vc.ein('xk,kwr->xwr', S2.Gi, S2.f)
        t = vc.ein('wk,xkr->xwr', S2.Gi, t)
        S2.fh = vc.ein('rk,xwk->xwr', S2.Gi, t)
        S2.TT = vc.ein('abe,cde->abcd', A.c, S2.f)
        okW = True
        for G in ('133', '223'):
            W2 = S2.W_full(G)
            # W2 / s2^7 == (W / s^7) / c  <=>  W2 * s^7 * c == W * s2^7
            okW &= bool(np.array_equal(W2.astype(object) * (S.s ** 7) * c, Wt[G].astype(object) * (s2 ** 7)))
        okB = bool(np.array_equal(vc.beta_tensor(S2.TT, P), c * B0))
        res.append("c = %d: W_133, W_223 -> W/c: %s, beta -> c beta: %s" % (c, okW, okB))
    print("   %-4s %s" % (name, " ; ".join(res)))

# (3)
A = valg.sl(2)
S = vc.Setup(A)
kap = S.kap
pairings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
rng = np.random.default_rng(4)
while True:
    Ps = [vc.config4(rng) for _ in range(6)]
    if sp.Matrix([vc.quad_basis_eval(P) for P in Ps]).det() != 0:
        break
brs = list(itertools.combinations(range(4), 2))
cols = []
for (pa, pb) in pairings:
    T = np.einsum('ab,cd->abcd', kap, kap)
    T = np.transpose(T, list(np.argsort([pa[0], pa[1], pb[0], pb[1]])))
    for m1, m2 in itertools.combinations_with_replacement(range(6), 2):
        poly = lambda Q, m1=m1, m2=m2: vc.brk(Q[brs[m1][0]], Q[brs[m1][1]]) * vc.brk(Q[brs[m2][0]], Q[brs[m2][1]])
        v = []
        for P in Ps:
            out = np.zeros((3,) * 4, dtype=np.int64)
            for s in itertools.permutations(range(4)):
                k = vc.perm_sign(s) * poly([P[s[0]], P[s[1]], P[s[2]], P[s[3]]])
                if k:
                    out += k * np.transpose(T, list(np.argsort(s)))
            v.append(out.reshape(-1))
        cols.append(np.concatenate(v))
Bv = np.concatenate([vc.beta_tensor(S.TT, P).reshape(-1) for P in Ps])
M = np.stack(cols, axis=1)
print("(3) sl2: rank of the %d cochains alt(mu kappa_pi) = %d ; with beta: %d ; beta nonzero: %s ; "
      "dim (g^(x)4)^g = 3 (Weyl), kappa_pi independent: rank %d"
      % (len(cols), valg.exact_rank(M.tolist()), valg.exact_rank(np.concatenate([M, Bv[:, None]], axis=1).tolist()),
         bool(Bv.any()),
         valg.exact_rank([np.transpose(np.einsum('ab,cd->abcd', kap, kap), list(np.argsort([a[0], a[1], b[0], b[1]]))).reshape(-1).tolist() for a, b in pairings])))

# (4)
Nn, pi = sp.symbols('Nnorm pi')
F2, b = sp.symbols('F223 beta_K')
F1 = -F2 + b / 2
for lab, expr in [('c3 sign', (pi * Nn / 240) * (-3 * F1 + 2 * F2)), ('K_223 flipped', (pi * Nn / 240) * (-3 * F1 - 2 * F2)),
                  ('K_133 flipped', (pi * Nn / 240) * (3 * F1 + 2 * F2))]:
    e = sp.expand(expr)
    print("(4) %-14s O_2 = %s * F223 + (%s) * beta_K" % (lab, sp.factor(e.coeff(F2)), sp.factor(e.coeff(b))))
print("    check (Nnorm/4)(-(pi/20) w-coeff) : sum_G |Aut|^-1 K_G W_G = (1/4)(-pi/20 W133 + pi/30 W223) [12][34] = (pi/240)(-3 W133 + 2 W223) [12][34]:",
      sp.simplify(sp.Rational(1, 4) * (-pi / 20) - pi / 240 * (-3)) == 0 and sp.simplify(sp.Rational(1, 4) * (pi / 30) - pi / 240 * 2) == 0)
