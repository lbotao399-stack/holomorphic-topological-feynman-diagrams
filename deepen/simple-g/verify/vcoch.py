"""Verifier's cochain machinery (independent of the explorer's core.py).

Form used: kappa = K / g, K = Killing form, g = gcd of the entries of K (so kappa is an integer matrix).
f[a,b,c] = kappa([e_a,e_b],e_c);  kappa^{-1} = Gi / s  (Gi integer).
Index forms (anomaly-structure claim c3; primed pairs contracted with kappa^{-1}):
  W_133(A1..A4) = - f_xyz f_{A1 y'u} f_{A2 u'w} f_{A3 z's} f_{A4 s'r} f_{x'w'r'}
  W_223(A1..A4) = - f_xyz f_{A1 x'u} f_{A2 y'w} f_{A3 z's} f_{A4 s'r} f_{u'w'r'}
Implementation:  h[A][y,u] = sum_y' Gi[y,y'] f[A,y',u]   (= s * kappa^{yy'} f_{A y' u}),
  fh[x,w,r] = Gi Gi Gi f (all raised, times s^3),   PP[A1,A2][y,w] = (h[A1] h[A2])[y,w],
  s^7 W_133 = - sum f[x,y,z] PP[A1,A2][y,w] PP[A3,A4][z,r] fh[x,w,r]
  s^7 W_223 = - sum f[x,y,z] h[A1][x,u] h[A2][y,w] PP[A3,A4][z,r] fh[u,w,r]
Cochains are evaluated at numerical momentum configurations P = (p0,..,p3), p3 = -(p0+p1+p2):
  F_G(P)[a0..a3] = sum_{s in S4} sgn(s) [P_s0 P_s1][P_s2 P_s3] W_G[a_s0, a_s1, a_s2, a_s3]
  beta(P)[a0..a3] = sum_{i<j} (-1)^{i+j} [P_k P_l]^2 kappa([e_ai,e_aj],[e_ak,e_al])
Exact integer arithmetic: every product is done in float64 only if the abs-value bound of all partial sums
is < 2^52 (exact), else in int64 if the bound is < 2^62, else in Python integers (object arrays).
"""
import itertools
import math
import numpy as np

LIM52 = 2.0 ** 52
LIM62 = 2.0 ** 62
STATS = {'float': 0, 'int64': 0, 'object': 0, 'maxbound': 0.0}


def perm_sign(s):
    s = list(s)
    sg = 1
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            if s[i] > s[j]:
                sg = -sg
    return sg


def mm(A, B):
    """exact integer matrix product A (m x k) @ B (k x n)."""
    if A.dtype == object or B.dtype == object:
        STATS['object'] += 1
        return A.astype(object).dot(B.astype(object))
    Af = A.astype(np.float64)
    Bf = B.astype(np.float64)
    with np.errstate(all='ignore'):          # spurious FP flags of the macOS BLAS; results checked below
        bound = np.abs(Af) @ np.abs(Bf)
    assert np.isfinite(bound).all()
    mx = float(bound.max()) if bound.size else 0.0
    STATS['maxbound'] = max(STATS['maxbound'], mx)
    if mx < LIM52:
        STATS['float'] += 1
        with np.errstate(all='ignore'):
            C = Af @ Bf
        assert np.isfinite(C).all() and float(np.abs(C).max(initial=0.0)) <= mx
        Ci = np.rint(C)
        assert np.array_equal(Ci, C), "non-integer float product"
        return Ci.astype(np.int64)
    if mx < LIM62:
        STATS['int64'] += 1
        return A.astype(np.int64) @ B.astype(np.int64)
    STATS['object'] += 1
    return A.astype(object).dot(B.astype(object))


def ein(spec, A, B):
    """exact two-operand einsum via mm (contracted indices = those absent from the output)."""
    ins, out = spec.split('->')
    ia, ib = ins.split(',')
    con = [x for x in ia if x in ib and x not in out]
    fa = [x for x in ia if x not in con]
    fb = [x for x in ib if x not in con]
    assert all(x in out for x in fa + fb) and sorted(out) == sorted(fa + fb)
    dims = {}
    for x, d in zip(ia, A.shape):
        dims[x] = d
    for x, d in zip(ib, B.shape):
        dims[x] = d
    At = np.transpose(A, [ia.index(x) for x in fa + con]).reshape(
        int(np.prod([dims[x] for x in fa])) if fa else 1, int(np.prod([dims[x] for x in con])) if con else 1)
    Bt = np.transpose(B, [ib.index(x) for x in con + fb]).reshape(
        int(np.prod([dims[x] for x in con])) if con else 1, int(np.prod([dims[x] for x in fb])) if fb else 1)
    C = mm(np.ascontiguousarray(At), np.ascontiguousarray(Bt))
    C = C.reshape([dims[x] for x in fa + fb])
    return np.transpose(C, [(fa + fb).index(x) for x in out])


class Setup:
    def __init__(self, A):
        self.A = A
        n = A.n
        K = A.K
        g = 0
        for x in K.reshape(-1):
            g = math.gcd(g, abs(int(x)))
        self.g = g                          # kappa = K / g
        kap = K // g
        assert np.array_equal(kap * g, K)
        self.kap = kap
        self.f = np.einsum('abd,dc->abc', A.c, kap)
        import sympy as sp
        inv = sp.Matrix(kap.tolist()).inv()
        s = int(sp.ilcm(*[sp.Rational(x).q for x in inv]))
        self.s = s
        self.Gi = np.array((inv * s).tolist(), dtype=np.int64)
        assert np.array_equal(self.Gi @ kap, s * np.eye(n, dtype=np.int64))
        f, Gi = self.f, self.Gi
        self.h = ein('yk,aku->ayu', Gi, f)                       # h[A,y,u]
        t = ein('xk,kwr->xwr', Gi, f)
        t = ein('wk,xkr->xwr', Gi, t)
        self.fh = ein('rk,xwk->xwr', Gi, t)
        # kappa([e_a,e_b],[e_c,e_d]) = sum_e c[a,b,e] f[c,d,e]
        self.TT = ein('abe,cde->abcd', A.c, self.f)

    # ---------------------------------------------------------------- full 4-tensors (times s^7)
    def W_full(self, graph):
        f, h, fh = self.f, self.h, self.fh
        PP = ein('ayu,buw->abyw', h, h)                          # PP[A1,A2,y,w]
        if graph == '133':
            G4 = ein('xyz,xwr->ywzr', f, fh)
            T = ein('abyw,ywzr->abzr', PP, G4)
        elif graph == '223':
            R = ein('xyz,axu->azyu', f, h)                       # R[A1,z,y,u]
            M = ein('byw,uwr->yubr', h, fh)                      # M[y,u,A2,r]
            T = ein('azyu,yubr->abzr', R, M)                     # T[A1,A2,z,r]
        else:
            raise ValueError(graph)
        W = ein('abzr,cdzr->abcd', T, PP)
        return -W

    # ---------------------------------------------------------------- pointwise (vectors), Python ints
    def W_point(self, graph, X):
        """s^7 W_graph(X1..X4) for integer coordinate vectors X[0..3] (exact, object arithmetic)."""
        f = self.f.astype(object)
        h = self.h.astype(object)
        fh = self.fh.astype(object)
        H = [np.tensordot(np.array(x, dtype=object), h, axes=([0], [0])) for x in X]   # n x n each
        P34 = H[2].dot(H[3])
        if graph == '133':
            P12 = H[0].dot(H[1])
            # sum f[x,y,z] P12[y,w] P34[z,r] fh[x,w,r]
            T1 = np.tensordot(f, P12, axes=([1], [0]))          # [x,z,w]
            T2 = np.tensordot(T1, P34, axes=([1], [0]))         # [x,w,r]
            val = (T2 * fh).sum()
        else:
            # sum f[x,y,z] H0[x,u] H1[y,w] P34[z,r] fh[u,w,r]
            T1 = np.tensordot(f, H[0], axes=([0], [0]))         # [y,z,u]
            T2 = np.tensordot(T1, H[1], axes=([0], [0]))        # [z,u,w]
            T3 = np.tensordot(T2, P34, axes=([0], [0]))         # [u,w,r]
            val = (T3 * fh).sum()
        return -val


# ------------------------------------------------------------------------------- momenta
def brk(p, q):
    return int(p[0]) * int(q[1]) - int(p[1]) * int(q[0])


def config4(rng, R=9):
    p = [rng.integers(-R, R + 1, 2) for _ in range(3)]
    p.append(-(p[0] + p[1] + p[2]))
    return [tuple(int(t) for t in x) for x in p]


def quad_basis_eval(P):
    """(u^2, v^2, w^2, uv, uw, vw) at P, u = [p0p1], v = [p0p2], w = [p1p2]."""
    u, v, w = brk(P[0], P[1]), brk(P[0], P[2]), brk(P[1], P[2])
    return [u * u, v * v, w * w, u * v, u * w, v * w]


def alt_tensor(W, P):
    """sum_s sgn(s) [P_s0 P_s1][P_s2 P_s3] W[a_s0..a_s3] as a full tensor."""
    out = None
    for s in itertools.permutations(range(4)):
        k = perm_sign(s) * brk(P[s[0]], P[s[1]]) * brk(P[s[2]], P[s[3]])
        if k == 0:
            continue
        # value at (a0..a3) = W[a_s0, a_s1, a_s2, a_s3]  ->  np.transpose(W, argsort(s))
        V = np.transpose(W, list(np.argsort(s)))
        term = k * V
        out = term if out is None else out + term
    return np.zeros_like(W) if out is None else out


def beta_tensor(TT, P):
    out = None
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        cf = (-1) ** (i + j) * brk(P[k], P[l]) ** 2
        if cf == 0:
            continue
        V = np.transpose(TT, list(np.argsort([i, j, k, l])))
        out = cf * V if out is None else out + cf * V
    return np.zeros_like(TT) if out is None else out


def alt_point(Wfun, X, P):
    tot = 0
    for s in itertools.permutations(range(4)):
        k = perm_sign(s) * brk(P[s[0]], P[s[1]]) * brk(P[s[2]], P[s[3]])
        if k:
            tot += k * Wfun([X[s[0]], X[s[1]], X[s[2]], X[s[3]]])
    return tot


def beta_point(S, X, P):
    A = S.A
    tot = 0
    co = A.c.astype(object)
    kap = S.kap.astype(object)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        bij = np.tensordot(np.tensordot(np.array(X[i], dtype=object), co, axes=([0], [0])),
                           np.array(X[j], dtype=object), axes=([0], [0]))
        bkl = np.tensordot(np.tensordot(np.array(X[k], dtype=object), co, axes=([0], [0])),
                           np.array(X[l], dtype=object), axes=([0], [0]))
        tot += (-1) ** (i + j) * brk(P[k], P[l]) ** 2 * int(bij.dot(kap).dot(bkl))
    return tot
