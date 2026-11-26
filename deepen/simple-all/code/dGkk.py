"""d G_kk = 0 for every Lie algebra with invariant form: G_kk(y_0..y_3) = sum_s sgn(s) [q_s0 q_s1][q_s2 q_s3] kappa(Y_s0,Y_s1) kappa(Y_s2,Y_s3).
(dG)(x_0..x_4) = sum_{i<j} (-1)^{i+j} G([x_i,x_j], x_rest); kappa([X_i,X_j],X_k) = f(X_i,X_j,X_k) (totally antisymmetric).
dG = sum_{l<m} P_lm(p) kappa(X_l,X_m) f(X_a,X_b,X_c), {a<b<c} = complement of {l,m}: the 10 polynomials P_lm are computed
(p_4 = -(p_0+..+p_3)) and checked to vanish identically; then dG = 0 for every g."""
import itertools
import sympy as sp
P = [[sp.Symbol('p%d_%d' % (i, c)) for c in (1, 2)] for i in range(4)]
P.append([-(P[0][c] + P[1][c] + P[2][c] + P[3][c]) for c in range(2)])
br = lambda u, v: u[0] * v[1] - u[1] * v[0]
def sgn(seq):
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s
coef = {}
for i, j in itertools.combinations(range(5), 2):
    rest = [r for r in range(5) if r not in (i, j)]
    ys = [('B', [P[i][0] + P[j][0], P[i][1] + P[j][1]])] + [(r, P[r]) for r in rest]
    for s in itertools.permutations(range(4)):
        Y = [ys[k] for k in s]
        poly = br(Y[0][1], Y[1][1]) * br(Y[2][1], Y[3][1]) * sgn(s) * (-1) ** (i + j)
        pairs = [(Y[0][0], Y[1][0]), (Y[2][0], Y[3][0])]
        # the pair containing 'B' gives f(X_i,X_j,X_k) = sign * f(sorted); the other pair gives kappa(X_l,X_m)
        for a, b in (pairs, pairs[::-1]):
            if 'B' in a:
                k = a[1] if a[0] == 'B' else a[0]
                lm = tuple(sorted(b))
                trip = [i, j, k]
                srt = sorted(trip)
                coef[lm] = coef.get(lm, 0) + poly * sgn([srt.index(t) for t in trip])
                break
ok = all(sp.expand(v) == 0 for v in coef.values())
print("structures kappa_lm f_abc: %d ; all 10 coefficient polynomials vanish identically: %s" % (len(coef), ok))
