"""Referee check, own code (no project imports).
(A) sl3, kappa = tr: F_223 is not exact in the G-invariant subcomplex of alternating CE cochains of
    g[[z1,z2]] on sum p = 0, bidegree (2,2) (G-Reynolds => not exact at all). Certificate: exact rank
    increase over Q at sample points. Controls: F_133+F_223, beta, sl2 F_223.
(B) Index forms W_133, W_223 vs the trace formulas of flagship Thm 3.5 (sl3).
(C) Deformation: dD = 0, {s0,w} = -dw, class of (1/2){D,D} for sl3 = lambda [F_223] with lambda solved
    exactly; explicit coboundary of flagship Thm 4.6.
(D) Extension: G-invariant weight-(1,1) 3-cochains for so5 (and sl3 as control): closed subspace.
All arithmetic exact (Python ints / Fractions); invariant-tensor dimensions by nullity mod p
(upper bound on the Q-nullity) together with exact independent invariants (lower bound)."""
from fractions import Fraction as Fr
import itertools, random, sys, math
import numpy as np

random.seed(7)
LOG = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); LOG.append(s)

# ------------------------------------------------------------------ algebras (integer matrices)
def E(n, i, j):
    m = [[0] * n for _ in range(n)]; m[i][j] = 1; return np.array(m, dtype=object)

def sl(n):
    B = [E(n, i, j) for i in range(n) for j in range(n) if i != j]
    B += [E(n, i, i) - E(n, i + 1, i + 1) for i in range(n - 1)]
    return B

def so(n):
    return [E(n, i, j) - E(n, j, i) for i in range(n) for j in range(i + 1, n)]

def sp(m):
    n = 2 * m
    J = np.zeros((n, n), dtype=object);
    for i in range(m):
        J[i][m + i] = 1; J[m + i][i] = -1
    # sp = { X : X^T J + J X = 0 } = { J S : S symmetric } (J^{-1} = -J)
    B = []
    for i in range(n):
        for j in range(i, n):
            S = E(n, i, j) + (E(n, j, i) if i != j else 0 * E(n, i, j))
            B.append(J.dot(S))
    return B

def tr(X): return sum(X[i][i] for i in range(X.shape[0]))
def mul(X, Y): return X.dot(Y)
def com(X, Y): return X.dot(Y) - Y.dot(X)

def frmat(A): return [[Fr(x) for x in r] for r in A]

def inv(A):
    n = len(A); M = [list(map(Fr, A[i])) + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [M[r][k] - f * M[c][k] for k in range(2 * n)]
    return [r[n:] for r in M]

def rank(cols):
    """exact rank over Q of a list of column vectors (Fractions)"""
    M = [list(r) for r in zip(*cols)]  # rows
    if not M: return 0
    nr, nc = len(M), len(M[0]); r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [M[i][k] - f * M[r][k] for k in range(nc)]
        r += 1
        if r == nr: break
    return r

def solve_ls(cols, b):
    """exact: find x with sum x_k cols_k = b (assumes consistent); returns x and residual check"""
    nc = len(cols); nr = len(b)
    M = [[cols[k][i] for k in range(nc)] + [b[i]] for i in range(nr)]
    piv = []; r = 0
    for c in range(nc):
        p = next((i for i in range(r, nr) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [M[i][k] - f * M[r][k] for k in range(nc + 1)]
        piv.append(c); r += 1
    consistent = all(M[i][nc] == 0 for i in range(r, nr))
    x = [Fr(0)] * nc
    for i, c in enumerate(piv): x[c] = M[i][nc]
    return x, consistent

class Alg:
    def __init__(self, name, basis):
        self.name = name; self.B = basis; self.dim = len(basis)
        G = [[Fr(tr(mul(a, b))) for b in basis] for a in basis]
        self.Ginv = inv(G)
        # dual basis T^a = sum_b Ginv[a][b] T_b  (tr(T^a T_c) = delta)
        self.Bd = [sum((self.Ginv[a][b] * basis[b] for b in range(self.dim)), 0 * basis[0]) for a in range(self.dim)]
        for a in range(self.dim):
            for c in range(self.dim):
                assert tr(mul(self.Bd[a], basis[c])) == (1 if a == c else 0)
        # coordinates: X = sum_a tr(T^a X) T_a
    def coords(self, X): return [tr(mul(self.Bd[a], X)) for a in range(self.dim)]
    def f(self, X, Y, Z): return tr(mul(com(X, Y), Z))
    def dsym(self, X, Y, Z): return tr(mul(X, mul(Y, Z) + mul(Z, Y)))

# ------------------------------------------------------------------ invariant 3-tensors: nullity mod p
def inv3_nullity(A, p=1000003):
    n = A.dim
    # structure constants [T_a,T_b] = sum_c C[a][b][c] T_c
    C = [[A.coords(com(A.B[a], A.B[b])) for b in range(n)] for a in range(n)]
    def md(x): x = Fr(x); return x.numerator * pow(x.denominator, -1, p) % p
    rows = []
    idx = lambda a, b, c: (a * n + b) * n + c
    for y in range(n):
        for a in range(n):
            for b in range(n):
                for e in range(n):
                    row = {}
                    for c in range(n):
                        for (k, coef) in ((idx(c, b, e), C[y][a][c]), (idx(a, c, e), C[y][b][c]), (idx(a, b, c), C[y][e][c])):
                            if coef != 0: row[k] = (row.get(k, 0) + md(coef)) % p
                    if any(v for v in row.values()): rows.append(row)
    M = np.zeros((len(rows), n ** 3), dtype=np.int64)
    for i, r in enumerate(rows):
        for k, v in r.items(): M[i, k] = v
    # rank mod p
    M = M % p; r = 0; nr, nc = M.shape
    for c in range(nc):
        piv = None
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        piv = r + nz[0]
        M[[r, piv]] = M[[piv, r]]
        invv = pow(int(M[r, c]), -1, p)
        M[r] = (M[r] * invv) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r]) % p) % p
        r += 1
        if r == nr: break
    return n ** 3 - r

# ------------------------------------------------------------------ cochains
def brk(a, b): return a[0] * b[1] - a[1] * b[0]

def eta_factory(A, mono, T):
    """alternating 3-cochain alt(m(q0,q1) T(X0,X1,X2)); mono=(i,j): (q0^1)^i (q1^1)^(2-i)... general exps"""
    e01, e11, e02, e12 = mono   # exponents of q0^1, q1^1, q0^2, q1^2
    def m(q0, q1): return q0[0] ** e01 * q1[0] ** e11 * q0[1] ** e02 * q1[1] ** e12
    def eta(xs):
        tot = 0
        for perm in itertools.permutations(range(3)):
            s = perm_sign(perm)
            q = [xs[perm[k]][1] for k in range(3)]; X = [xs[perm[k]][0] for k in range(3)]
            tot += s * m(q[0], q[1]) * T(X[0], X[1], X[2])
        return tot
    return eta

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def d3(eta):
    """(d eta)(x0..x3) = sum_{i<j} (-1)^{i+j} eta([x_i,x_j], rest)"""
    def de(xs):
        tot = 0
        for i in range(4):
            for j in range(i + 1, 4):
                rest = [xs[k] for k in range(4) if k not in (i, j)]
                xij = (com(xs[i][0], xs[j][0]), (xs[i][1][0] + xs[j][1][0], xs[i][1][1] + xs[j][1][1]))
                tot += (-1) ** (i + j) * eta([xij] + rest)
        return tot
    return de

def d4(om):
    def de(xs):
        tot = 0
        for i in range(5):
            for j in range(i + 1, 5):
                rest = [xs[k] for k in range(5) if k not in (i, j)]
                xij = (com(xs[i][0], xs[j][0]), (xs[i][1][0] + xs[j][1][0], xs[i][1][1] + xs[j][1][1]))
                tot += (-1) ** (i + j) * om([xij] + rest)
        return tot
    return de

def bracket33(A, a, b):
    """{a,b}(x1..x4) = sum_{shuffles} sgn a(x_s1,x_s2,X_e) kappa^{ee'} b(X_e',x_s3,x_s4)"""
    sh = [(i, j, k, l) for (i, j) in itertools.combinations(range(4), 2) for (k, l) in [tuple(sorted(set(range(4)) - {i, j}))]]
    def br_(xs):
        tot = 0
        for (i, j, k, l) in sh:
            s = perm_sign((i, j, k, l))
            Q = (xs[i][1][0] + xs[j][1][0], xs[i][1][1] + xs[j][1][1])
            mQ = (-Q[0], -Q[1])
            for e in range(A.dim):
                tot += s * a([xs[i], xs[j], (A.B[e], mQ)]) * b([(A.Bd[e], Q), xs[k], xs[l]])
        return tot
    return br_

# ---- W index forms (flagship Prop 3.3), primed pairs contracted with kappa^{-1}
def W_tensors(A):
    n = A.dim
    F = np.array([[[A.f(A.B[a], A.B[b], A.B[c]) for c in range(n)] for b in range(n)] for a in range(n)], dtype=object)
    Fr_ = np.vectorize(Fr)
    F = Fr_(F)
    Gi = np.array(A.Ginv, dtype=object)
    # raise index helper: Fup[a', ...] = sum_a Gi[a', a] ... ; we need contractions kappa^{x x'}
    # f_{x y z} with all three raised: Fuuu
    Fuuu = np.einsum('xyz,xX,yY,zZ->XYZ', F, Gi, Gi, Gi)
    # W133 = - f_{xyz} f_{a1 y' u} f_{a2 u' w} f_{a3 z' s} f_{a4 s' r} f_{x' w' r'}
    # contract: G1[a1,y,u] = f_{a1 y' u} with y' -> raised; easier: use Fuuu on N0 (x,y,z raised) and lower others.
    W133 = -np.einsum('XYZ,aYu,uU,bUw,cZs,sS,dSr,Xwr->abcd',
                      Fuuu, F, Gi, F, F, Gi, F,
                      np.einsum('xwr,wW,rR->xWR', F, Gi, Gi), optimize=True)
    # W223 = - f_{xyz} f_{a1 x' u} f_{a2 y' w} f_{a3 z' s} f_{a4 s' r} f_{u' w' r'}
    Fnp = np.einsum('uwr,uU,wW,rR->UWR', F, Gi, Gi, Gi)
    W223 = -np.einsum('XYZ,aXu,bYw,cZs,sS,dSr,uwr->abcd', Fuuu, F, F, F, Gi, F, Fnp, optimize=True)
    out = []
    for W in (W133, W223):
        den = 1
        for v in W.flat: den = den * Fr(v).denominator // math.gcd(den, Fr(v).denominator)
        Wi = np.array([int(Fr(v) * den) for v in W.flat], dtype=object).reshape(W.shape)
        mx = max(abs(int(v)) for v in Wi.flat)
        assert mx * n ** 4 * 3 ** 4 < 9 * 10 ** 18
        out.append((W, np.array(Wi, dtype=np.int64), den))
    return out[0], out[1]

def W_eval(Wt, cs):
    W, Wi, den = Wt
    v = np.einsum('abcd,a,b,c,d->', Wi, *[np.array([int(x) for x in c], dtype=np.int64) for c in cs])
    return Fr(int(v), den)

def F_eval(A, W, xs):
    """F_G = sum_s sgn(s) [p_s1 p_s2][p_s3 p_s4] W(X_s1..X_s4), xs=[(X,p)...]"""
    cs = [A.coords(x[0]) for x in xs]
    tot = 0
    for s in itertools.permutations(range(4)):
        tot += perm_sign(s) * brk(xs[s[0]][1], xs[s[1]][1]) * brk(xs[s[2]][1], xs[s[3]][1]) * \
            W_eval(W, [cs[s[0]], cs[s[1]], cs[s[2]], cs[s[3]]])
    return tot

def beta_eval(A, xs):
    tot = 0
    for i in range(4):
        for j in range(i + 1, 4):
            k, l = sorted(set(range(4)) - {i, j})
            tot += (-1) ** (i + j) * brk(xs[k][1], xs[l][1]) ** 2 * tr(mul(com(xs[i][0], xs[j][0]), com(xs[k][0], xs[l][0])))
    return tot

def rand_pt(A, k, R=3, Rp=6):
    xs = []
    ps = [(random.randint(-Rp, Rp), random.randint(-Rp, Rp)) for _ in range(k - 1)]
    ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
    for i in range(k):
        X = sum((random.randint(-R, R) * b for b in A.B), 0 * A.B[0])
        xs.append((X, ps[i]))
    return xs

monos22 = [(i, 2 - i, j, 2 - j) for i in range(3) for j in range(3)]
monos11 = [(i, 1 - i, j, 1 - j) for i in range(2) for j in range(2)]

def F_trace_sl(xs, N, which):
    """flagship Thm 3.5 sl_N formulas, X's traceless"""
    X = [x[0] for x in xs]
    def t(*ix):
        M = X[ix[0] - 1]
        for i in ix[1:]: M = mul(M, X[i - 1])
        return tr(M)
    if which == '133':
        return (2 * t(1, 2, 3, 4) + (2 + N * N) * t(1, 2, 4, 3) - 4 * t(1, 3, 2, 4) + (2 + N * N) * t(1, 3, 4, 2)
                - 4 * t(1, 4, 2, 3) + 2 * t(1, 4, 3, 2) + 6 * N * t(1, 2) * t(3, 4))
    if which == '223':
        return (-2 * t(1, 2, 3, 4) - 2 * t(1, 2, 4, 3) + 4 * t(1, 3, 2, 4) - 2 * t(1, 3, 4, 2) + 4 * t(1, 4, 2, 3)
                - 2 * t(1, 4, 3, 2) - 4 * N * t(1, 2) * t(3, 4) + 2 * N * t(1, 3) * t(2, 4) + 2 * N * t(1, 4) * t(2, 3))
    if which == 'tot':
        return (-10 * t(1, 2, 3, 4) - (10 + 3 * N * N) * t(1, 2, 4, 3) + 20 * t(1, 3, 2, 4) - (10 + 3 * N * N) * t(1, 3, 4, 2)
                + 20 * t(1, 4, 2, 3) - 10 * t(1, 4, 3, 2) - 26 * N * t(1, 2) * t(3, 4) + 4 * N * t(1, 3) * t(2, 4) + 4 * N * t(1, 4) * t(2, 3))

# ================================================================== run
def run_algebra(A, invT, nsamp, extra=None):
    P('=' * 70); P('algebra', A.name, 'dim', A.dim, 'kappa = tr')
    W133, W223 = W_tensors(A)
    pts = [rand_pt(A, 4) for _ in range(nsamp)]
    etas = [(mono, Tn, eta_factory(A, mono, T)) for mono in monos22 for (Tn, T) in invT]
    dcols = [[Fr(d3(e)(xs)) for xs in pts] for (_, _, e) in etas]
    rd = rank(dcols)
    P('G-inv C^3_(2,2) generators:', len(etas), '| rank of d on them at', nsamp, 'sample points:', rd)
    F133 = [Fr(F_eval(A, W133, xs)) for xs in pts]
    F223 = [Fr(F_eval(A, W223, xs)) for xs in pts]
    beta = [Fr(beta_eval(A, xs)) for xs in pts]
    S = [a + b for a, b in zip(F133, F223)]
    P('rank[d | F_223]          =', rank(dcols + [F223]), '(rank increase => F_223 not exact)')
    P('rank[d | F_133]          =', rank(dcols + [F133]))
    P('rank[d | -3F133+2F223]   =', rank(dcols + [[-3 * a + 2 * b for a, b in zip(F133, F223)]]))
    P('rank[d | F_133+F_223]    =', rank(dcols + [S]))
    P('rank[d | beta]           =', rank(dcols + [beta]))
    # beta = d([q0q1]^2 f)/6 check
    ef = None
    def eqq(xs):
        tot = 0
        for perm in itertools.permutations(range(3)):
            q = [xs[perm[k]][1] for k in range(3)]; X = [xs[perm[k]][0] for k in range(3)]
            tot += perm_sign(perm) * brk(q[0], q[1]) ** 2 * A.f(X[0], X[1], X[2])
        return tot
    dq = [Fr(d3(eqq)(xs)) for xs in pts]
    P('d(alt([q0q1]^2 f)) == 6 beta at all samples:', all(a == 6 * b for a, b in zip(dq, beta)))
    # F_133+F_223 = c beta ?
    nzb = [i for i in range(nsamp) if beta[i] != 0]
    c = S[nzb[0]] / beta[nzb[0]]
    P('F_133+F_223 == c*beta at all samples:', all(S[i] == c * beta[i] for i in range(nsamp)), ' c =', c,
      '(kappa = tr; Killing = r tr gives c_Killing = c/r^2)')
    return dict(W133=W133, W223=W223, pts=pts, dcols=dcols, F133=F133, F223=F223, beta=beta, etas=etas, rd=rd)

# ---------------- sl2 control
sl2 = Alg('sl2', sl(2))
P('inv 3-tensor nullity mod p, sl2:', inv3_nullity(sl2))
R2 = run_algebra(sl2, [('f', sl2.f)], 30)
# ---------------- sl3
sl3 = Alg('sl3', sl(3))
P('inv 3-tensor nullity mod p, sl3:', inv3_nullity(sl3))
R3 = run_algebra(sl3, [('f', sl3.f), ('d', sl3.dsym)], 50)

# (B) trace formulas, sl3
pts = R3['pts']
ok = {w: all(Fr(W_eval(R3['W' + w], [sl3.coords(x[0]) for x in xs])) == F_trace_sl(xs, 3, w) for xs in pts[:10]) for w in ('133', '223')}
okt = all(Fr(-3 * W_eval(R3['W133'], [sl3.coords(x[0]) for x in xs]) + 2 * W_eval(R3['W223'], [sl3.coords(x[0]) for x in xs])) == F_trace_sl(xs, 3, 'tot') for xs in pts[:10])
P('index-form W_133 == Thm 3.5 sl_N formula (N=3, 10 random points):', ok['133'])
P('index-form W_223 == Thm 3.5 sl_N formula (N=3, 10 random points):', ok['223'])
P('index-form -3W133+2W223 == Thm 3.5 W_tot first nine terms (N=3):', okt)
# sl2 formula
pts2 = R2['pts']
def k(xs, i, j): return tr(mul(xs[i - 1][0], xs[j - 1][0]))
ok2 = all(Fr(-3 * W_eval(R2['W133'], [sl2.coords(x[0]) for x in xs]) + 2 * W_eval(R2['W223'], [sl2.coords(x[0]) for x in xs]))
          == -104 * k(xs, 1, 2) * k(xs, 3, 4) + 16 * k(xs, 1, 3) * k(xs, 2, 4) + 40 * k(xs, 1, 4) * k(xs, 2, 3) for xs in pts2[:10])
P('sl2 W_tot == -104 k12k34 + 16 k13k24 + 40 k14k23 (10 points):', ok2)
# sl2: F_133, F_223 vs beta in Killing normalization: Killing(sl2) = 4 tr.
b2 = R2['beta']; i0 = next(i for i in range(len(b2)) if b2[i] != 0)
c133 = R2['F133'][i0] / b2[i0]; c223 = R2['F223'][i0] / b2[i0]
P('sl2 (tr): F_133 = %s beta_tr, F_223 = %s beta_tr, exact at all samples: %s' % (c133, c223,
  all(R2['F133'][i] == c133 * b2[i] and R2['F223'][i] == c223 * b2[i] for i in range(len(b2)))))
P('sl2 converted to Killing (kappa=4tr => F->F/4, beta->4beta): F_133 = %s beta_K, F_223 = %s beta_K' % (c133 / 16, c223 / 16))

# ================================================================== (C) deformation, sl3
A = sl3
def Dco(xs):   # D = [p1p2] tr(1{23}), arity 3, xs=[(X,q)...]
    return brk(xs[0][1], xs[1][1]) * A.dsym(xs[0][0], xs[1][0], xs[2][0])
def s0(xs):
    return tr(mul(xs[0][0], com(xs[1][0], xs[2][0])))
pts3 = [rand_pt(A, 3) for _ in range(20)]
P('D alternating (20 pts, all 6 perms):', all(Dco([xs[p[0]], xs[p[1]], xs[p[2]]]) == perm_sign(p) * Dco(xs) for xs in pts3 for p in itertools.permutations(range(3))))
pts4 = R3['pts']
P('dD == 0 at all', len(pts4), 'samples:', all(d3(Dco)(xs) == 0 for xs in pts4))
# {s0,w} = -dw for w = a random G-inv 3-cochain of bidegree (2,2)
wtest = lambda xs: R3['etas'][4][2](xs) + 3 * R3['etas'][13][2](xs)
bs = bracket33(A, s0, wtest); dw = d3(wtest)
P('{s0,w} == -dw (w = eta_4 + 3 eta_13), 15 samples:', all(Fr(bs(xs)) == -Fr(dw(xs)) for xs in pts4[:15]),
  '| nonzero:', any(dw(xs) != 0 for xs in pts4[:15]))
P('{s0,s0} == 0, 10 samples:', all(Fr(bracket33(A, s0, s0)(xs)) == 0 for xs in pts4[:10]))
bDD = bracket33(A, Dco, Dco)
half = [Fr(bDD(xs)) / 2 for xs in pts4]
dcols = R3['dcols']
P('rank[d | (1/2){D,D}] =', rank(dcols + [half]), '(rank of d alone =', R3['rd'], ')')
x, cons = solve_ls(dcols + [R3['F223']], half)
P('(1/2){D,D} in span(d C^3_inv) + Q F_223 at samples (consistent):', cons, '| coefficient of F_223 =', x[-1],
  '| flagship a^2/(12N^2) at N=3: 1/108 ->', x[-1] == Fr(1, 108))
# explicit coboundary: (1/2){D,D} = F223/(12N^2) + d[ (N^2+4)/(2N^2) [p1p2]^2 tr(1[23]) ]
N = 3
def etaN(xs): return Fr(N * N + 4, 2 * N * N) * brk(xs[0][1], xs[1][1]) ** 2 * A.f(xs[0][0], xs[1][0], xs[2][0])
deN = d3(etaN)
lhs_ok = all(half[i] == R3['F223'][i] / (12 * N * N) + Fr(deN(pts4[i])) for i in range(len(pts4)))
lhs_ok_m = all(half[i] == R3['F223'][i] / (12 * N * N) - Fr(deN(pts4[i])) for i in range(len(pts4)))
P('explicit identity (1/2){D,D} = F_223/108 + d eta_N at all samples:', lhs_ok, '| with -d eta_N:', lhs_ok_m)
# double trace claim: [alt([p1p2][p3p4] tr12 tr34)] = -[F_223]/(6N)
def G4(xs):
    tot = 0
    for s in itertools.permutations(range(4)):
        tot += perm_sign(s) * brk(xs[s[0]][1], xs[s[1]][1]) * brk(xs[s[2]][1], xs[s[3]][1]) * \
            tr(mul(xs[s[0]][0], xs[s[1]][0])) * tr(mul(xs[s[2]][0], xs[s[3]][0]))
    return tot
Gc = [Fr(G4(xs)) for xs in pts4]
x2, cons2 = solve_ls(dcols + [R3['F223']], Gc)
P('alt([p1p2][p3p4]tr12tr34) class = %s [F_223] (consistent %s); flagship -1/(6N) = %s' % (x2[-1], cons2, Fr(-1, 18)))
# closedness of F's (5-point d), a few samples
pts5 = [rand_pt(A, 5) for _ in range(6)]
F223f = lambda xs: F_eval(A, R3['W223'], xs); F133f = lambda xs: F_eval(A, R3['W133'], xs)
P('dF_223 == 0, dF_133 == 0 at 6 five-point samples (sl3):', all(d4(F223f)(xs) == 0 for xs in pts5), all(d4(F133f)(xs) == 0 for xs in pts5))

# ================================================================== (D) weight-(1,1) G-invariant 3-cochains
import sympy as sp_
q = sp_.symbols('a1 a2 b1 b2')   # q0=(a1,a2), q1=(b1,b2), q2=-(q0+q1)
def qv(k):
    if k == 0: return (q[0], q[1])
    if k == 1: return (q[2], q[3])
    return (-q[0] - q[2], -q[1] - q[3])
syms = []
for (e01, e11, e02, e12) in monos11:
    tot = 0
    for perm in itertools.permutations(range(3)):
        q0, q1 = qv(perm[0]), qv(perm[1])
        tot += q0[0] ** e01 * q1[0] ** e11 * q0[1] ** e02 * q1[1] ** e12
    syms.append(sp_.Poly(sp_.expand(tot), *q))
mons = sorted(set(m for p in syms for m in p.monoms()))
Msym = sp_.Matrix([[p.as_dict().get(m, 0) for m in mons] for p in syms])
P('rank of S3-symmetrizations of the 4 weight-(1,1) monomials (exact, sympy):', Msym.rank(),
  '=> alt(m f) spans a space of dim <= this for any g with (T^3)^g = C f')
for (nm, basis) in (('sl3', sl(3)), ('so5', so(5))):
    B_ = Alg(nm, basis)
    nul = inv3_nullity(B_)
    invT = [('f', B_.f)] + ([('d', B_.dsym)] if nm.startswith('sl') else [])
    ptsD3 = [rand_pt(B_, 3, R=2) for _ in range(25)]
    ptsD4 = [rand_pt(B_, 4, R=2) for _ in range(25)]
    etas = [eta_factory(B_, mono, T) for mono in monos11 for (_, T) in invT]
    ecols = [[Fr(e(xs)) for xs in ptsD3] for e in etas]
    dc = [[Fr(d3(e)(xs)) for xs in ptsD4] for e in etas]
    # closed combinations: kernel of d restricted to span(etas)
    r_e = rank(ecols); r_d = rank(dc)
    P('%s: nullity mod p of invariance on T^3 = %d; generators %d; dim span C^3_inv,(1,1) >= %d; rank d >= %d'
      % (nm, nul, len(etas), r_e, r_d))
    if nm.startswith('so') or nm.startswith('sp'):
        P('   => if dim span == rank d, no nonzero closed G-inv weight-(1,1) 3-cochain:', r_e == r_d)

open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/notes/flagship/referee/ce_check.log', 'w').write('\n'.join(LOG) + '\n')

# ================================================================== (E) other weights for so5; E_f for so7, sp6
monos10 = [(1, 0, 0, 0), (0, 1, 0, 0)]
monos20 = [(2, 0, 0, 0), (1, 1, 0, 0), (0, 2, 0, 0)]
for (nm, basis) in (('so5', so(5)),):
    B_ = Alg(nm, basis)
    for (wn, ms) in (('(1,0)', monos10), ('(2,0)', monos20)):
        ptsD3 = [rand_pt(B_, 3, R=2) for _ in range(25)]
        ptsD4 = [rand_pt(B_, 4, R=2) for _ in range(25)]
        etas = [eta_factory(B_, mono, B_.f) for mono in ms]
        ecols = [[Fr(e(xs)) for xs in ptsD3] for e in etas]
        dc = [[Fr(d3(e)(xs)) for xs in ptsD4] for e in etas]
        P('%s weight %s: G-inv generators alt(m f): %d; sample rank of cochains %d; sample rank of d %d'
          % (nm, wn, len(etas), rank(ecols), rank(dc)))
for (nm, basis) in (('so7', so(7)), ('sp6', sp(3))):
    B_ = Alg(nm, basis)
    def Ef(xs, B_=B_):
        tot = 0
        for perm in itertools.permutations(range(3)):
            q = [xs[perm[k]][1] for k in range(3)]; X = [xs[perm[k]][0] for k in range(3)]
            tot += perm_sign(perm) * q[0][0] * q[0][1] * B_.f(X[0], X[1], X[2])
        return tot
    ptsD4 = [rand_pt(B_, 4, R=2) for _ in range(4)]
    vals = [d3(Ef)(xs) for xs in ptsD4]
    P('%s: d(alt(q0^1 q0^2 f)) at 4 sample points: %s  (nonzero => the unique G-inv weight-(1,1) 3-cochain is not closed, given dim (T^3)^g = 1 from flagship Prop 3.10)' % (nm, vals))
open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/notes/flagship/referee/ce_check.log', 'w').write('\n'.join(LOG) + '\n')
