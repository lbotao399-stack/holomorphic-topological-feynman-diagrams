# referee2: so9 (new in run 3), kappa = tr_9 (Killing = 7 tr_9).  Own construction, no project imports.
# (1) F133 + F223 = (49/2) beta_tr ; (2) F223 = 3 beta_tr - (15/2) G_kk  (= table (x,a) = (3/49, -15/686) in Killing normalization);
# (3) rank certificate: rank[beta | F223] = 2, rank[beta | G_kk] = 2 by exact nonzero 2x2 minors  =>  F223, G_kk not in C beta;
# (4) Prop. 5.7 table entries for so9: G-invariant C^3_(a,b) = Sym(m) f, dims and ranks of d at (1,1), (2,0), (0,2), (2,2), (1,2).
import sys, random, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 99)
out = []
def P(*s):
    line = ' '.join(str(x) for x in s); print(line, flush=True); out.append(line)

n = 9
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
dim = len(pairs)
T = np.zeros((dim, n, n), dtype=np.int64)
for a, (i, j) in enumerate(pairs):
    T[a, i, j] = 1; T[a, j, i] = -1
gram = np.einsum('aij,bji->ab', T, T)
assert (gram == -2*np.eye(dim, dtype=np.int64)).all()      # kappa = tr: kappa_ab = -2 delta_ab, kappa^{-1} = -1/2 delta
fT = np.einsum('aij,bjk,cki->abc', T, T, T) - np.einsum('bij,ajk,cki->abc', T, T, T)   # tr([Ta,Tb]Tc)
P("so9: dim", dim, "; nonzero f entries", int((fT != 0).sum()), "; Jacobi check:",
  bool((np.einsum('abx,xcd->abcd', fT, fT) + np.einsum('bcx,xad->abcd', fT, fT) + np.einsum('cax,xbd->abcd', fT, fT) == 0).all()))

def Mx(X):  # M[b,c] = tr([X,T_b] T_c)
    return np.einsum('ij,bjk,cki->bc', X, T, T) - np.einsum('bij,jk,cki->bc', T, X, T)

def W(graph, Xs):
    M = [Mx(X) for X in Xs]
    if graph == '223':
        PA, sA = M[0], Fr(1, 4); PB, sB = M[1], Fr(1, 4); PC, sC = M[2] @ M[3], Fr(-1, 8)
    else:
        PA, sA = np.eye(dim, dtype=np.int64), Fr(-1, 2); PB, sB = M[0] @ M[1], Fr(-1, 8); PC, sC = M[2] @ M[3], Fr(-1, 8)
    Q = np.einsum('XYZ,xX,yY,zZ->xyz', fT, PA, PB, PC, optimize=True)
    S = int((fT * Q).sum())
    Qf = np.einsum('XYZ,xX,yY,zZ->xyz', fT.astype(float), PA.astype(float), PB.astype(float), PC.astype(float), optimize=True)
    assert abs(float((fT*Qf).sum()) - S) <= 1e-6*max(1, abs(S)), "overflow"
    return -sA*sB*sC*S

def br(a, b): return a[0]*b[1] - a[1]*b[0]
def tr(m): return int(np.trace(m))
def comm(a, b): return a @ b - b @ a
def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def FG(graph, xs):
    s = Fr(0)
    for perm in itertools.permutations(range(4)):
        pp = [xs[i][0] for i in perm]; XX = [xs[i][1] for i in perm]
        s += perm_sign(perm) * br(pp[0], pp[1]) * br(pp[2], pp[3]) * W(graph, XX)
    return s

def Gkk(xs):
    s = 0
    for perm in itertools.permutations(range(4)):
        pp = [xs[i][0] for i in perm]; XX = [xs[i][1] for i in perm]
        s += perm_sign(perm) * br(pp[0], pp[1]) * br(pp[2], pp[3]) * tr(XX[0] @ XX[1]) * tr(XX[2] @ XX[3])
    return s

def dco(om, xs):
    k = len(xs); s = 0
    for i in range(k):
        for j in range(i+1, k):
            rest = [xs[m] for m in range(k) if m not in (i, j)]
            s += (-1)**(i+j) * om(((xs[i][0][0]+xs[j][0][0], xs[i][0][1]+xs[j][0][1]), comm(xs[i][1], xs[j][1])), *rest)
    return s

def u2f(x0, x1, x2): return br(x0[0], x1[0])**2 * tr(x0[1] @ comm(x1[1], x2[1]))
def beta(xs): return dco(u2f, xs)

def rmat():
    A = np.zeros((n, n), dtype=np.int64)
    for i, j in pairs:
        v = random.randint(-2, 2); A[i, j] = v; A[j, i] = -v
    return A
def rpoint(k):
    ps = [(random.randint(-4, 4), random.randint(-4, 4)) for _ in range(k-1)]
    ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
    return [(p, rmat()) for p in ps]

P("== (1), (2) at 4 random points (exact)")
vals = []
for t in range(4):
    x = rpoint(4)
    f133, f223, b, g = FG('133', x), FG('223', x), beta(x), Gkk(x)
    vals.append((f223, b, g))
    P(f" point {t}: F133+F223 == (49/2) beta: {f133+f223 == Fr(49, 2)*b} ; F223 == 3 beta - (15/2) G_kk: {f223 == 3*b - Fr(15, 2)*g} ;"
      f" F223 = {f223}, beta = {b}, G_kk = {g}")
P("== (3) rank certificates (exact 2x2 minors over the 4 points)")
m1 = vals[0][0]*vals[1][1] - vals[0][1]*vals[1][0]
m2 = vals[0][2]*vals[1][1] - vals[0][1]*vals[1][2]
P(" minor det[[F223(P0), beta(P0)], [F223(P1), beta(P1)]] =", m1, " -> rank[beta|F223] = 2:", m1 != 0)
P(" minor det[[G_kk(P0), beta(P0)], [G_kk(P1), beta(P1)]] =", m2, " -> rank[beta|G_kk] = 2:", m2 != 0)
P(" Killing-normalized a = -15/686 = (3/16) v1v2v3/t^3 with Vogel (v1,v2,v3,t) for so_N = (-2, 4, N-4, N-2):",
  Fr(3, 16)*Fr(-2*4*(n-4), (n-2)**3) == Fr(-15, 686), "; a(so_N) = -3(N-4)/(2(N-2)^3) at N=9:", Fr(-3*(n-4), 2*(n-2)**3))

P("== (4) Prop. 5.7 for so9: C^3_(a,b) = Sym(m) f; rank of d by exact minors at random points")
q01, q11, q02, q12 = sp.symbols('q01 q11 q02 q12')   # q_i^c : component c of momentum i  (q0 = (q01,q02)? no: q0 = (q01, q02), q1 = (q11, q12))
Q = {0: (q01, q02), 1: (q11, q12), 2: (-q01-q11, -q02-q12)}
def sym_space(a, b):
    mons = []
    for i0 in range(a+1):
        for j0 in range(b+1):
            mons.append(Q[0][0]**i0 * Q[1][0]**(a-i0) * Q[0][1]**j0 * Q[1][1]**(b-j0))
    syms = []
    for m in mons:
        s = 0
        for perm in itertools.permutations(range(3)):
            s += m.subs({q01: sp.Symbol('A1'), q02: sp.Symbol('A2'), q11: sp.Symbol('B1'), q12: sp.Symbol('B2')}, simultaneous=True).subs(
                {sp.Symbol('A1'): Q[perm[0]][0], sp.Symbol('A2'): Q[perm[0]][1], sp.Symbol('B1'): Q[perm[1]][0], sp.Symbol('B2'): Q[perm[1]][1]}, simultaneous=True)
        syms.append(sp.Poly(sp.expand(s), q01, q11, q02, q12))
    # basis via row reduction of coefficient vectors
    allm = sorted({mm for p in syms for mm in p.as_dict()})
    Mt = sp.Matrix([[p.as_dict().get(mm, 0) for mm in allm] for p in syms])
    rref, piv = Mt.T.rref()
    return [syms[i] for i in piv]
for (a, b) in [(1, 1), (2, 0), (0, 2), (2, 2), (1, 2)]:
    B = sym_space(a, b)
    def cochain(poly):
        def c(x0, x1, x2):
            v = poly.eval({q01: x0[0][0], q02: x0[0][1], q11: x1[0][0], q12: x1[0][1]})
            return int(v) * tr(x0[1] @ comm(x1[1], x2[1]))
        return c
    rows = []
    for t in range(len(B) + 2):
        x = rpoint(4)
        rows.append([dco(cochain(pb), x) for pb in B])
    rk = sp.Matrix(rows).rank()
    # alternation sanity of the cochains
    x = rpoint(3)
    alt_ok = all(cochain(pb)(*x) == -cochain(pb)(x[1], x[0], x[2]) == -cochain(pb)(x[0], x[2], x[1]) for pb in B)
    P(f" ({a},{b}): dim C^3 = {len(B)}  alternating: {alt_ok}  rank d_3 >= {rk} (exact rank of the value matrix) -> Z^3 = 0: {rk == len(B)}")

with open(sys.argv[2] if len(sys.argv) > 2 else 'r2_so9.log', 'w') as fh:
    fh.write('\n'.join(out) + '\n')
