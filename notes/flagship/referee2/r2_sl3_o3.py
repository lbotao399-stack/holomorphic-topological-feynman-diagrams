# referee2: sl3, kappa = tr.  Own colour factors (index networks of flagship Prop. 4.4 / eq. after Prop. 6.2), own CE differential,
# own shuffle bracket, own evaluation of o', theta'(d), {D,iota}, d o_3 and the order-hbar^3 primitive.  Kernels K^(v) taken from the
# table of flagship Prop. 6.4 (three of them recomputed in r2_kernels.py).  No project imports.  Exact integers / Fractions.
import sys, random, itertools
from fractions import Fraction as Fr
import numpy as np

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
out = []
def P(*s):
    line = ' '.join(str(x) for x in s); print(line, flush=True); out.append(line)

N = 3
def E(i, j):
    m = np.zeros((3, 3), dtype=object); m[:] = 0; m[i, j] = 1; return m
basis = [E(0, 1), E(0, 2), E(1, 0), E(1, 2), E(2, 0), E(2, 1)]
H1 = E(0, 0) - E(1, 1); H2 = E(1, 1) - E(2, 2)
basis += [H1, H2]
dim = 8
def tr(m): return sum(m[i, i] for i in range(3))
def mm(a, b): return a.dot(b)
def comm(a, b): return mm(a, b) - mm(b, a)
kap = np.array([[tr(mm(a, b)) for b in basis] for a in basis], dtype=object)
import sympy
kinv = sympy.Matrix(kap.tolist()).inv()
kinv3 = np.array([[int(3*kinv[i, j]) for j in range(dim)] for i in range(dim)], dtype=object)
assert all(3*kinv[i, j] == int(3*kinv[i, j]) for i in range(dim) for j in range(dim))
kinv = {(i, j): Fr(int(3*kinv[i, j]), 3) for i in range(dim) for j in range(dim)}
fT = np.array([[[tr(mm(comm(a, b), c)) for c in basis] for b in basis] for a in basis], dtype=object)
dT = np.array([[[tr(mm(a, mm(b, c) + mm(c, b))) for c in basis] for b in basis] for a in basis], dtype=object)

def coords(X):
    # X = sum_{i!=j} X_ij E_ij + h1 H1 + h2 H2 ; h1 = X00, h2 = -X22
    v = [X[0, 1], X[0, 2], X[1, 0], X[1, 2], X[2, 0], X[2, 1], X[0, 0], -X[2, 2]]
    return np.array(v, dtype=object)

# ---- 4-tensors of the theta networks, local order L' (flagship Prop. 4.4, Prop. 6.4 table) ----
def ein(spec, *ops):
    return np.einsum(spec, *[np.array(o, dtype=np.int64) for o in ops], optimize=True)

def W_tensor(graph, v, T):
    """W^{(v,T)}_G = sgn(L->L') prod_{u != v} f(L'_u) T(L'_v) prod_e kappa^{-1}, sgn = -1; v = None: all f (W_G).
    Returns integer tensor of 3^7 * W, slots (x1,x2,x3,x4)."""
    ten = {u: (T if u == v else fT) for u in range(6)}
    k = kinv3
    if graph == '133':
        # A: e0 (0->5); B: 0->1->2->5 (e1,e2,e3); C: 0->3->4->5 (e4,e5,e6)
        PA = np.array(k, dtype=np.int64)
        PB = ein('xb,pbc,cd,qde,ey->xpqy', k, ten[1], k, ten[2], k)
        PC = ein('xb,pbc,cd,qde,ey->xpqy', k, ten[3], k, ten[4], k)
        W = ein('xyz,xX,ypqY,zrsZ,XYZ->pqrs', ten[0], PA, PB, PC, ten[5])
    else:
        # A: 0->1->5 (e0,e1); B: 0->2->5 (e2,e3); C: 0->3->4->5 (e4,e5,e6)
        PA = ein('xb,pbc,cy->xpy', k, ten[1], k)
        PB = ein('xb,pbc,cy->xpy', k, ten[2], k)
        PC = ein('xb,pbc,cd,qde,ey->xpqy', k, ten[3], k, ten[4], k)
        W = ein('xyz,xpX,yqY,zrsZ,XYZ->pqrs', ten[0], PA, PB, PC, ten[5])
    Wf = None
    return -np.array(W, dtype=object)

# overflow guard: recompute one tensor in float64 and compare
def W_tensor_float(graph):
    k = np.array(kinv3, float); f = np.array(fT, float)
    if graph == '133':
        PB = np.einsum('xb,pbc,cd,qde,ey->xpqy', k, f, k, f, k)
        W = np.einsum('xyz,xX,ypqY,zrsZ,XYZ->pqrs', f, k, PB, PB, f, optimize=True)
    return -W

SC = 3**7
WG = {g: W_tensor(g, None, fT) for g in ('133', '223')}
assert np.allclose(np.array(WG['133'], float), W_tensor_float('133')), "int64 overflow"
P("max |3^7 W| entries:", max(abs(int(x)) for x in WG['133'].ravel()), max(abs(int(x)) for x in WG['223'].ravel()))
WV = {(g, v): W_tensor(g, v, dT) for g in ('133', '223') for v in range(6)}

def Wev(Wt, Xs):
    c = [coords(X) for X in Xs]
    r = np.tensordot(Wt, c[3], axes=([3], [0]))
    r = np.tensordot(r, c[2], axes=([2], [0]))
    r = np.tensordot(r, c[1], axes=([1], [0]))
    return Fr(int(np.dot(r, c[0])), SC)

# ---- cochains: points are lists of (p, X), p a 2-tuple, sum p = 0 ----
def br(a, b): return a[0]*b[1] - a[1]*b[0]
def add2(a, b): return (a[0]+b[0], a[1]+b[1])
def neg2(a): return (-a[0], -a[1])

def Dc(x1, x2, x3):   return br(x1[0], x2[0]) * tr(mm(x1[1], mm(x2[1], x3[1]) + mm(x3[1], x2[1])))
def U2F(x1, x2, x3):  return br(x1[0], x2[0])**2 * tr(mm(x1[1], comm(x2[1], x3[1])))
def D3c(x1, x2, x3):  return br(x1[0], x2[0])**3 * tr(mm(x1[1], mm(x2[1], x3[1]) + mm(x3[1], x2[1])))

def dco(om):
    def r(*xs):
        k = len(xs); s = 0
        for i in range(k):
            for j in range(i+1, k):
                rest = [xs[m] for m in range(k) if m not in (i, j)]
                s += (-1)**(i+j) * om((add2(xs[i][0], xs[j][0]), comm(xs[i][1], xs[j][1])), *rest)
        return s
    return r

def perm_sign(p):
    s = 1; p = list(p)
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s

def bracket(a, b):   # arity 3, 3 -> 4 ; flagship Sec. 5.3
    def r(*xs):
        s = Fr(0)
        for I in itertools.combinations(range(4), 2):
            J = [m for m in range(4) if m not in I]
            sg = perm_sign(list(I) + J)
            q = neg2(add2(xs[I[0]][0], xs[I[1]][0]))
            va = [a(xs[I[0]], xs[I[1]], (q, T)) for T in basis]
            vb = [b((neg2(q), T), xs[J[0]], xs[J[1]]) for T in basis]
            s += sg * sum(Fr(va[e]) * kinv[e, f] * vb[f] for e in range(dim) for f in range(dim) if kinv[e, f] != 0)
        return s
    return r

def lin(*terms):  # linear combination of cochains
    def r(*xs): return sum(Fr(c) * f(*xs) for c, f in terms)
    return r
def scale(c, f): return lambda *xs: Fr(c) * f(*xs)

def FG(Wt):
    def r(*xs):
        s = Fr(0)
        for perm in itertools.permutations(range(4)):
            sg = perm_sign(perm)
            pp = [xs[i][0] for i in perm]; XX = [xs[i][1] for i in perm]
            s += sg * br(pp[0], pp[1]) * br(pp[2], pp[3]) * Wev(Wt, XX)
        return s
    return r

def b_of(pp):
    b = {}
    for i in range(4):
        for j in range(4): b[(i+1, j+1)] = br(pp[i], pp[j])
    return b

KT = {('133', 0): lambda b: Fr(1, 252) * b[1,2]*(b[1,3]+b[2,3])*(b[1,2]-2*b[1,3]-b[2,3]),
      ('133', 1): lambda b: Fr(1, 420) * b[1,2]*(b[1,3]+b[2,3])*(7*b[1,2]-b[1,3]),
      ('133', 2): lambda b: Fr(1, 420) * b[1,2]*(b[1,3]+b[2,3])*(8*b[1,2]-b[2,3]),
      ('133', 3): lambda b: Fr(1, 420) * b[1,2]*(b[1,3]+b[2,3])*(8*b[1,3]+7*b[2,3]),
      ('133', 4): lambda b: Fr(-1, 420) * b[1,2]*(b[1,3]+b[2,3])*(b[1,2]-7*b[1,3]-8*b[2,3]),
      ('133', 5): lambda b: Fr(-1, 252) * b[1,2]*(b[1,3]+b[2,3])*(2*b[1,2]-b[1,3]-2*b[2,3]),
      ('223', 0): lambda b: Fr(1, 840) * b[1,2]*(b[1,3]+b[2,3])*(5*b[1,2]-2*b[1,3]+2*b[2,3]),
      ('223', 1): lambda b: Fr(1, 420) * b[1,2]*(b[1,3]+b[2,3])*(b[1,2]+2*b[1,3]),
      ('223', 2): lambda b: Fr(-1, 420) * b[1,2]*(b[1,3]+b[2,3])*(b[1,2]-2*b[2,3]),
      ('223', 3): lambda b: Fr(-1, 84) * b[1,2]*(b[1,3]+b[2,3])**2,
      ('223', 4): lambda b: Fr(-1, 84) * b[1,2]*(b[1,3]+b[2,3])**2,
      ('223', 5): lambda b: Fr(1, 840) * b[1,2]*(b[1,3]+b[2,3])*(9*b[1,2]+2*b[1,3]-2*b[2,3])}

def thetap(*xs):  # theta'(d) = sum_G |Aut|^{-1} sum_v alt(K^(v)/pi W^(v,d))
    s = Fr(0)
    for perm in itertools.permutations(range(4)):
        sg = perm_sign(perm)
        pp = [xs[i][0] for i in perm]; XX = [xs[i][1] for i in perm]
        b = b_of(pp)
        for key, K in KT.items():
            kv = K(b)
            if kv != 0: s += sg * Fr(1, 4) * kv * Wev(WV[key], XX)
    return s

def rmat():
    while True:
        X = np.array([[random.randint(-4, 4) for _ in range(3)] for _ in range(3)], dtype=object)
        X[2, 2] = -X[0, 0] - X[1, 1]
        if any(X[i, j] != 0 for i in range(3) for j in range(3)): return X
def rpoint(k):
    ps = [(random.randint(-5, 5), random.randint(-5, 5)) for _ in range(k-1)]
    ps.append((-sum(p[0] for p in ps), -sum(p[1] for p in ps)))
    return [(p, rmat()) for p in ps]

def trw(X, order): # tr(X_{o1} X_{o2} X_{o3} X_{o4})
    m = X[order[0]-1]
    for o in order[1:]: m = mm(m, X[o-1])
    return tr(m)

# -------- colour-factor claims (Prop. 6.5, Thm. 4.6 at N = 3) --------
P("== colour factors at 6 random X-tuples")
okc = True; ok223 = True; okth = True
for t in range(6):
    Xs = [rmat() for _ in range(4)]
    for v in range(6):
        w133 = Wev(WV[('133', v)], Xs); w223 = Wev(WV[('223', v)], Xs)
        tgt = N**2 * (trw(Xs, (1, 3, 4, 2)) - trw(Xs, (1, 2, 4, 3))) * (1 if v < 3 else -1)
        okc &= (w133 == tgt); ok223 &= (w223 == 0)
    W133 = Wev(WG['133'], Xs); W223 = Wev(WG['223'], Xs)
    t2 = lambda i, j: tr(mm(Xs[i-1], Xs[j-1]))
    W133c = 2*trw(Xs,(1,2,3,4)) + (2+N**2)*trw(Xs,(1,2,4,3)) - 4*trw(Xs,(1,3,2,4)) + (2+N**2)*trw(Xs,(1,3,4,2)) - 4*trw(Xs,(1,4,2,3)) + 2*trw(Xs,(1,4,3,2)) + 6*N*t2(1,2)*t2(3,4)
    W223c = -2*trw(Xs,(1,2,3,4)) - 2*trw(Xs,(1,2,4,3)) + 4*trw(Xs,(1,3,2,4)) - 2*trw(Xs,(1,3,4,2)) + 4*trw(Xs,(1,4,2,3)) - 2*trw(Xs,(1,4,3,2)) - 4*N*t2(1,2)*t2(3,4) + 2*N*t2(1,3)*t2(2,4) + 2*N*t2(1,4)*t2(2,3)
    okth &= (W133 == W133c and W223 == W223c)
P("W_133^(v,d) = +-N^2(tr(1342)-tr(1243)) (+ for v=0,1,2):", okc, "; W_223^(v,d) = 0 all v:", ok223,
  "; Thm 4.6 sl_N trace formulas for W_133, W_223 at N=3:", okth)

# -------- alternation sanity --------
x = rpoint(3)
P("D, u2f, D3 alternating (one transposition):", Dc(*x) == -Dc(x[1], x[0], x[2]), U2F(*x) == -U2F(x[1], x[0], x[2]), D3c(*x) == -D3c(x[0], x[2], x[1]))

F133 = FG(WG['133']); F223 = FG(WG['223'])
beta = dco(U2F)
DD = bracket(Dc, Dc)
op = lin((Fr(-3, 240), F133), (Fr(2, 240), F223))
iota = scale(Fr(-(3*N**2+10), 20), U2F)
Diota = bracket(Dc, iota)
dD3 = dco(D3c)
o3n = lin((1, Diota), (1, thetap))            # o_3 / (a pi nu)
xi = scale(Fr(-(2*N**2+5), 30), D3c)

checks = {
 'F133+F223 = 2N^2 beta':                      (lin((1, F133), (1, F223)), scale(2*N**2, beta)),
 "o' = F223/48 - (N^2/40) beta":               (op, lin((Fr(1, 48), F223), (Fr(-N**2, 40), beta))),
 '(1/2){D,D} = F223/(12N^2) + d(((N^2+4)/(2N^2)) u2f)': (scale(Fr(1, 2), DD), lin((Fr(1, 12*N**2), F223), (Fr(N**2+4, 2*N**2), beta))),
 "o' - (N^2/8){D,D} = d iota, iota = -((3N^2+10)/20) u2f": (lin((1, op), (Fr(-N**2, 8), DD)), dco(iota)),
 "theta'(d) = -(N^2/60) dD3":                  (thetap, scale(Fr(-N**2, 60), dD3)),
 '{D,u2f} = (1/3) dD3':                        (bracket(Dc, U2F), scale(Fr(1, 3), dD3)),
 '{D,iota} + theta\' = d xi, xi = -((2N^2+5)/30) D3': (o3n, dco(xi)),
 "{D,iota} - theta' = d(-((N^2+5)/30) D3)":    (lin((1, Diota), (-1, thetap)), scale(Fr(-(N**2+5), 30), dD3)),
}
P("== cochain identities at 6 random 4-point configurations (sum p = 0), exact")
for name, (lhs, rhs) in checks.items():
    res = []
    for t in range(6):
        x = rpoint(4)
        l, r = lhs(*x), rhs(*x)
        res.append((l == r, l))
    P(f"{name}: {all(a for a, _ in res)}   (lhs values: {[str(v) for _, v in res[:3]]} ...)")

P("== controls (must fail)")
x = rpoint(4)
P("theta'(d) = +(N^2/60) dD3 :", thetap(*x) == scale(Fr(N**2, 60), dD3)(*x), "; {D,iota} - theta' = d xi :",
  lin((1, Diota), (-1, thetap))(*x) == dco(xi)(*x), "; {D,iota} - theta' = d(((N^2+5)/30) D3):", lin((1, Diota), (-1, thetap))(*x) == scale(Fr(N**2+5, 30), dD3)(*x))

P("== closedness: d o_3 = 0 (5-point configurations), d o' = 0; control d(random 4-cochain) != 0")
for t in range(3):
    y = rpoint(5)
    P(" d(o_3/(a pi nu)) =", dco(o3n)(*y), "; d o' =", dco(op)(*y), "; d theta' =", dco(thetap)(*y),
      "; control d(F223 * [p0p1]-weighted) =", dco(lambda *xs: br(xs[0][0], xs[1][0]) * F223(*xs))(*y))

P("== order-hbar^3 coefficients (N = 3)")
P(" i_2 = eps pi nu iota = eps U ((3N^2+10)/1280) pi^-4 u2f :", Fr(3*N**2+10, 1280), "  (pi nu = -U/(64 pi^4); -(-(3N^2+10)/20)/64 =", Fr(3*N**2+10, 20*64), ")")
P(" i_3 = a eps pi nu xi = a eps U ((2N^2+5)/1920) pi^-4 D3 :", Fr(2*N**2+5, 1920), "  (= -(-(2N^2+5)/30)/64 =", Fr(2*N**2+5, 30*64), ")")

with open(sys.argv[2] if len(sys.argv) > 2 else 'r2_sl3_o3.log', 'w') as fh:
    fh.write('\n'.join(out) + '\n')
