"""Miscellaneous exact checks.
 (1) Killing / tr_V for the matrix constructions (all entries).
 (2) d-tensor tr(X{Y,Z}) nonzero for sl_N (lower bound dim (S^3 g*)^g >= 1).
 (3) d G_kk = 0 symbolically: G_kk = alt([p0p1][p2p3] k01 k23);  dG(x0..x4) = sum_{l<m} P_lm(p) k_lm f_abc,
     computed with f([X_i,X_j],X_k) = sgn f_abc and the 10 coefficient polynomials P_lm (sympy).
 (4) closedness of alt([p0p1][p2p3] W) for a non-symmetric invariant W = tr(X0 X1 X2 X3) on sl3 (random points mod p),
     and a control: the same with the QB coefficients permuted (not closed).
 (5) Vogel-parameter identities (sympy) and comparison with the root-data n4 of vbk.
"""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
import vlie, vpoint, vbk

out = []
# (1)
for key, rr in [('so5', 3), ('so7', 5), ('so8', 6), ('so9', 7), ('so10', 8), ('so11', 9), ('so12', 10), ('sp4', 6), ('sp6', 8), ('sp8', 10),
                ('sp10', 12), ('sl3', 6), ('sl4', 8), ('sl5', 10), ('sl6', 12), ('sl7', 14), ('so6', 4), ('G2', None), ('F4', None)]:
    L = vlie.build(key)
    if L.kappaV is None:
        continue
    ratios = set()
    for a in range(L.d):
        for b in range(L.d):
            kv = L.kappaV[a][b]; K = L.K[a][b]
            if kv != 0:
                ratios.add(K / kv)
            elif K != 0:
                ratios.add('inconsistent')
    out.append('(1) %-10s Killing / tr_V on all entries: %s (expected %s)' % (L.name, sorted(ratios, key=str), rr))
# G2 in so7: Killing / tr_7 on G2
L = vlie.build('G2')
ratios = set(L.K[a][b] / L.kappaV[a][b] for a in range(L.d) for b in range(L.d) if L.kappaV[a][b] != 0)
out.append('(1) G2(phi) Killing / tr_7: %s' % sorted(ratios))
# (2)
rng = np.random.default_rng(1)
for N in range(3, 8):
    X, Y, Z = [rng.integers(-3, 4, (N, N)) for _ in range(3)]
    X, Y, Z = [M - np.trace(M) * np.eye(N, dtype=int) * 0 for M in (X, Y, Z)]
    for M in (X, Y, Z):
        M[N - 1, N - 1] -= np.trace(M)
    dval = np.trace(X @ (Y @ Z + Z @ Y))
    out.append('(2) sl%d: tr(X{Y,Z}) at a random traceless triple = %d (nonzero %s)' % (N, dval, dval != 0))
# (3)
p = [sp.symbols('q%d_1 q%d_2' % (i, i)) for i in range(5)]
Bk = lambda a, b: a[0] * b[1] - a[1] * b[0]
def Gpoly(qs, s):
    return vlie.perm_sign(s) * Bk(qs[s[0]], qs[s[1]]) * Bk(qs[s[2]], qs[s[3]])
# G(y0..y3) = sum_s sgn(s)[q_s0 q_s1][q_s2 q_s3] k(Y_s0,Y_s1) k(Y_s2,Y_s3)
# dG(x0..x4) = sum_{i<j} (-1)^{i+j} G([x_i,x_j], x_k, x_l, x_m)
coef = {}
for i, j in itertools.combinations(range(5), 2):
    rest = [k for k in range(5) if k not in (i, j)]
    labels = ['B'] + rest            # slot 0 = [x_i, x_j]
    qs = [(p[i][0] + p[j][0], p[i][1] + p[j][1])] + [p[k] for k in rest]
    for s in itertools.permutations(range(4)):
        pairs = [(labels[s[0]], labels[s[1]]), (labels[s[2]], labels[s[3]])]
        # the pair containing 'B' gives k([X_i,X_j], X_k) = f(X_i,X_j,X_k); the other pair gives k_lm
        pb = [pr for pr in pairs if 'B' in pr][0]
        po = [pr for pr in pairs if 'B' not in pr][0]
        kk = [x for x in pb if x != 'B'][0]
        trip = (i, j, kk)
        lm = tuple(sorted(po))
        sg = vlie.perm_sign(sorted(range(3), key=lambda t: trip[t]))   # f(X_i,X_j,X_k) = sgn * f(sorted)
        key = lm
        coef[key] = coef.get(key, 0) + (-1) ** (i + j) * sg * Gpoly(qs, s)
sub = {p[4][c]: -(p[0][c] + p[1][c] + p[2][c] + p[3][c]) for c in range(2)}
allzero = all(sp.expand(v.subs(sub, simultaneous=True)) == 0 for v in coef.values())
out.append('(3) d G_kk: %d coefficient polynomials P_lm (k_lm f_abc), all identically 0 after sum p = 0: %s' % (len(coef), allzero))
# (4)
L = vlie.build('sl3'); P_ = 1048573
E = vpoint.Engine(L, P_)
mats = None
# sl3 basis matrices (rebuild): use coordinates of vlie.sl
N = 3
basis = []
for i in range(N):
    for j in range(N):
        if i != j:
            M = np.zeros((N, N), dtype=np.int64); M[i, j] = 1; basis.append(M)
for i in range(N - 1):
    M = np.zeros((N, N), dtype=np.int64); M[i, i] = 1; M[i + 1, i + 1] = -1; basis.append(M)
def mat(X):
    return sum(int(x) * B for x, B in zip(X, basis)) % P_
def Wtr(Xs):
    M = np.eye(N, dtype=np.int64)
    for X in Xs:
        M = M @ mat(X) % P_
    return int(np.trace(M)) % P_
def Fcoch(Ys, Q, perm_mu=None):
    u, v, w = [(int(Q[a][0]) * int(Q[b][1]) - int(Q[a][1]) * int(Q[b][0])) % P_ for a, b in ((0, 1), (0, 2), (1, 2))]
    mono = [u * u, v * v, w * w, u * v, u * w, v * w]
    if perm_mu:
        mono = [mono[k] for k in perm_mu]
    tot = 0
    for s in vpoint.PERMS:
        sg, c = vpoint.CQ[s]
        val = sum(cm * mm for cm, mm in zip(c, mono)) % P_
        tot = (tot + sg * val * Wtr([Ys[s[0]], Ys[s[1]], Ys[s[2]], Ys[s[3]]])) % P_
    return tot
rng = np.random.default_rng(9)
for trial in range(3):
    X = [rng.integers(0, P_, L.d).astype(np.float64) for _ in range(5)]
    Pm = [rng.integers(0, P_, 2) for _ in range(4)]; Pm.append(np.array([(-sum(q[0] for q in Pm)) % P_, (-sum(q[1] for q in Pm)) % P_]))
    res = []
    for perm_mu in (None, (3, 1, 2, 0, 4, 5)):
        tot = 0
        for i, j in itertools.combinations(range(5), 2):
            rest = [k for k in range(5) if k not in (i, j)]
            Ys = [E.bracket(X[i], X[j])] + [X[k] for k in rest]
            Q = [(Pm[i] + Pm[j]) % P_] + [Pm[k] for k in rest]
            tot = (tot + (-1) ** (i + j) * Fcoch(Ys, Q, perm_mu)) % P_
        res.append(tot)
    out.append('(4) sl3, W = tr(X0X1X2X3): d alt([p0p1][p2p3] W) = %d ; control (QB coefficients uu<->uv permuted) = %d' % tuple(res))
# (5)
al, be, ga = sp.symbols('alpha beta gamma')
t = al + be + ga
t2 = al ** 2 + be ** 2 + ga ** 2; t3 = al ** 3 + be ** 3 + ga ** 3
dV = (al - 2 * t) * (be - 2 * t) * (ga - 2 * t) / (al * be * ga)
C4 = (3 * t * t2 - t3) / (16 * t ** 3)
n4V = dV * (C4 - sp.Rational(1, 8))
e1 = sp.simplify(8 * n4V - 3 * dV + 3 * dV * (t + al) * (t + be) * (t + ga) / (2 * t ** 3))
e2 = sp.simplify(8 * n4V - 3 * dV - 3 * dV * (dV - 3) * al * be * ga / (4 * t ** 3))
out.append('(5) 8 n4 - 3d + 3d(t+a)(t+b)(t+c)/(2t^3) = %s ; 8 n4 - 3d - 3d(d-3)abc/(4t^3) = %s' % (e1, e2))
params = {('A', 1): (-2, 2, 2), ('G', 2): (-2, Fr(10, 3), Fr(8, 3)), ('F', 4): (-2, 5, 6), ('E', 6): (-2, 6, 8), ('E', 7): (-2, 8, 12), ('E', 8): (-2, 12, 20)}
for N in range(5, 31):
    params[('B', (N - 1) // 2) if N % 2 else ('D', N // 2)] = (-2, 4, N - 4)
for N in range(4, 31, 2):
    params[('C', N // 2)] = (-2, 1, Fr(N, 2) + 2)
for N in range(3, 21):
    params[('A', N - 1)] = (-2, 2, N)
okall = True; bad = []
for (typ, n), (a_, b_, c_) in sorted(params.items()):
    r = vbk.analyse(typ, n)
    sub = {al: sp.Rational(Fr(a_).numerator, Fr(a_).denominator), be: sp.Rational(Fr(b_).numerator, Fr(b_).denominator),
           ga: sp.Rational(Fr(c_).numerator, Fr(c_).denominator)}
    dv = dV.subs(sub); nv = n4V.subs(sub)
    ok = (dv == r['d']) and (nv == sp.Rational(r['n'][4].numerator, r['n'][4].denominator))
    okall &= ok
    if not ok:
        bad.append((typ, n))
out.append('(5) Vogel formula d, n4 = d(C_4 - 1/8) with (alpha,beta,gamma): sl_N (-2,2,N), so_N (-2,4,N-4), sp_N (-2,1,N/2+2), G2 (-2,10/3,8/3), F4 (-2,5,6), E6 (-2,6,8), E7 (-2,8,12), E8 (-2,12,20): equal to root-data d, n4 for %d algebras: %s %s'
           % (len(params), okall, bad))
print('\n'.join(out))
