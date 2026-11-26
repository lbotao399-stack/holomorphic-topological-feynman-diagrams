"""Two-loop obstruction cochains of CS on R x C^2 for a Lie algebra L (algebras.LA), evaluated exactly.

Colour weights (anomaly-structure claim c3 index forms; f_abc = kappa([X_a,X_b],X_c), primed pairs
contracted with kappa^{-1}):
  W_133(A1..A4) = - f_xyz f_{A1 y'u} f_{A2 u'w} f_{A3 z's} f_{A4 s'r} f_{x'w'r'}
  W_223(A1..A4) = - f_xyz f_{A1 x'u} f_{A2 y'w} f_{A3 z's} f_{A4 s'r} f_{u'w'r'}
With G = s kappa^{-1} (integer), h[a][b,c] = sum_k G[b,k] f[a,k,c], fhat = G G G f (all indices raised):
  W_133 = -s^{-7} sum f[x,y,z] (h_A1 h_A2)[y,w] (h_A3 h_A4)[z,r] fhat[x,w,r]
  W_223 = -s^{-7} sum f[x,y,z] h_A1[x,u] h_A2[y,w] (h_A3 h_A4)[z,r] fhat[u,w,r]
Integer arrays below are s^7 W.

Cochains (functions of x_i = X_i e^{p_i z}, sum p = 0, 0-based slots 0..3):
  F_G      = sum_sigma sgn(sigma) [p_s0 p_s1][p_s2 p_s3] W_G(X_s0, X_s1, X_s2, X_s3)
  beta     = d( [q1 q2]^2 f(Y1,Y2,Y3) ) = sum_{i<j} (-1)^{i+j} [p_k p_l]^2 kappa([X_i,X_j],[X_k,X_l])
Momentum polynomials are reduced exactly to the basis (u^2, v^2, w^2, uv, uw, vw),
u = [p0 p1], v = [p0 p2], w = [p1 p2], p3 = -(p0+p1+p2).
Evaluation set: X_0 = b_theta, X_1, X_2, X_3 = all basis triples.  For g simple an invariant tensor T
vanishes iff T(b_theta, ., ., .) = 0 (X -> T(X,.,.,.) is g-equivariant, its kernel is an ideal), so this
evaluation is injective on (quadratics) x (g^{*(x)4})^g.
"""
import itertools
import math
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from algebras import perm_sign

LIM = 2.0 ** 62

# ----------------------------------------------------------------------------- momentum algebra
QB = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']


def br(a, b):
    """[p_a p_b] as coefficient vector over (u, v, w), 0 <= a, b <= 3, p3 = -(p0+p1+p2)."""
    if a == b:
        return np.zeros(3, dtype=np.int64)
    if a > b:
        return -br(b, a)
    base = {(0, 1): 0, (0, 2): 1, (1, 2): 2}
    if b < 3:
        v = np.zeros(3, dtype=np.int64); v[base[(a, b)]] = 1
        return v
    return -sum(br(a, c) for c in range(3))


def quad(l1, l2):
    """product of two linear forms over (u,v,w) -> coefficients over QB."""
    out = np.zeros(6, dtype=np.int64)
    pos = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3, (0, 2): 4, (1, 2): 5}
    for i in range(3):
        for j in range(3):
            if l1[i] and l2[j]:
                out[pos[(min(i, j), max(i, j))]] += l1[i] * l2[j]
    return out


def check_momentum_algebra():
    """sympy: the reduction rules are identities and the 6 quadratics are independent."""
    P = [[sp.Symbol('p%d_%d' % (i, c)) for c in range(2)] for i in range(3)]
    P.append([-(P[0][c] + P[1][c] + P[2][c]) for c in range(2)])
    B = lambda a, b: P[a][0] * P[b][1] - P[a][1] * P[b][0]
    u, v, w = B(0, 1), B(0, 2), B(1, 2)
    lin = [u, v, w]
    ok = True
    for a in range(4):
        for b in range(4):
            if sp.expand(B(a, b) - sum(int(c) * x for c, x in zip(br(a, b), lin))) != 0:
                ok = False
    quads = [u * u, v * v, w * w, u * v, u * w, v * w]
    gens = [P[i][c] for i in range(3) for c in range(2)]
    rows = []
    allm = set()
    polys = [sp.Poly(sp.expand(q), *gens) for q in quads]
    for p in polys:
        allm |= set(p.as_dict())
    allm = sorted(allm)
    M = sp.Matrix([[p.as_dict().get(m, 0) for m in allm] for p in polys])
    return ok, M.rank()


# ----------------------------------------------------------------------------- integer contractions
LIM53 = 2.0 ** 52


def ct(spec, A, B):
    """exact integer einsum of two integer operands.  bound = einsum(|A|,|B|) bounds every partial sum.
    bound < 2^52: float64 (BLAS) is exact (all partial sums and products are integers below 2^53);
    bound < 2^62: int64; otherwise OverflowError."""
    Af, Bf = A.astype(np.float64), B.astype(np.float64)
    bound = np.einsum(spec, np.abs(Af), np.abs(Bf), optimize=True)
    mx = bound.max() if bound.size else 0.0
    if mx < LIM53:
        return np.rint(np.einsum(spec, Af, Bf, optimize=True)).astype(np.int64)
    if mx >= LIM:
        raise OverflowError("int64 bound exceeded: %.3e" % mx)
    return np.einsum(spec, A, B, optimize=True)


class Tensors:
    def __init__(self, L):
        self.L = L
        f, G = L.fl, L.G
        self.f = f
        self.h = ct('bk,akc->abc', G, f)
        t = ct('xk,kwr->xwr', G, f)
        t = ct('wk,xkr->xwr', G, t)
        self.fhat = ct('rk,xwk->xwr', G, t)
        self.s7 = L.s ** 7

    def H(self, slot):
        """slot: None (free), basis index (fixed) or integer coordinate vector (fixed)."""
        if slot is None:
            return self.h
        if isinstance(slot, np.ndarray):
            return ct('a,abc->bc', slot.astype(np.int64), self.h)[None]
        return self.h[slot][None]

    def W(self, graph, slots):
        """s^7 W_graph with slots[i] = basis index or None (free); free axes in slot order.
        Loops over the slot-0 index when slot 0 is free (memory <= dim^4)."""
        if slots[0] is None:
            out = [self.W(graph, [a] + list(slots[1:])) for a in range(self.L.dim)]
            return np.stack(out, axis=0)
        H1, H2, H3, H4 = [self.H(x) for x in slots]
        f, fhat = self.f, self.fhat
        MC = ct('ayu,buw->abyw', H3, H4)              # (n3, n4, z, r)
        if graph == '133':
            MB = ct('ayu,buw->abyw', H1, H2)          # (n1, n2, y, w)
            Q = ct('xyz,abyw->abxzw', f, MB)
            R = ct('abxzw,xwr->abzr', Q, fhat)
        elif graph == '223':
            P = ct('xyz,axu->ayzu', f, H1)
            P2 = ct('ayzu,byw->abzuw', P, H2)
            R = ct('abzuw,uwr->abzr', P2, fhat)
        else:
            raise ValueError(graph)
        W = -ct('abzr,cdzr->abcd', R, MC)
        sq = tuple(i for i, x in enumerate(slots) if x is not None)
        sq = tuple(i for i in sq if W.shape[i] == 1)
        return W.reshape(W.shape).squeeze(axis=sq)


def F_eval(T, graph, theta):
    """F_graph(b_theta, ., ., .) as int64 array (6, dim, dim, dim) over QB (times s^7)."""
    d = T.L.dim
    Wj = []
    for j in range(4):
        slots = [None] * 4
        slots[j] = theta
        Wj.append(T.W(graph, slots))
    out = np.zeros((6, d, d, d), dtype=np.int64)
    for sig in itertools.permutations(range(4)):
        sg = perm_sign(sig)
        c = quad(br(sig[0], sig[1]), br(sig[2], sig[3]))
        j = sig.index(0)
        ks = [k for k in range(4) if k != j]
        perm = [next(m for m in range(3) if sig[ks[m]] == t + 1) for t in range(3)]
        V = np.transpose(Wj[j], perm)
        for mu in range(6):
            if c[mu]:
                out[mu] += sg * int(c[mu]) * V
    return out


def tree_full(L):
    """TT[a,b,c,d] = kappa([b_a,b_b],[b_c,b_d])."""
    return ct('abg,cdg->abcd', L.fl, L.fu)


def beta_eval(L, theta, TT=None):
    if TT is None:
        TT = tree_full(L)
    d = L.dim
    out = np.zeros((6, d, d, d), dtype=np.int64)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = quad(br(k, l), br(k, l))
        # Tvar[v0..v3] = TT[v_i, v_j, v_k, v_l]
        order = [i, j, k, l]
        Tvar = np.transpose(TT, list(np.argsort(order)))
        V = Tvar[theta]
        for mu in range(6):
            if c[mu]:
                out[mu] += (-1) ** (i + j) * int(c[mu]) * V
    return out


# ----------------------------------------------------------------------------- exact linear algebra
def to_rows(cols):
    """stack flattened integer columns, keep rows where some column is nonzero; Python ints."""
    A = np.stack([c.reshape(-1) for c in cols], axis=1)
    keep = np.any(A != 0, axis=1)
    return [[int(x) for x in row] for row in A[keep]]


def rank_exact(rows):
    """exact rank over Q (fraction-free Gaussian elimination on Python ints)."""
    M = [list(r) for r in rows if any(r)]
    if not M:
        return 0
    ncol = len(M[0])
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pr = M[r]
        new = M[:r + 1]
        for i in range(r + 1, len(M)):
            if M[i][c]:
                a, b = pr[c], M[i][c]
                row = [a * x - b * y for x, y in zip(M[i], pr)]
                g = 0
                for x in row:
                    g = math.gcd(g, abs(x))
                if g > 1:
                    row = [x // g for x in row]
                if any(row):
                    new.append(row)
            else:
                new.append(M[i])
        M = new
        r += 1
        if r == len(M):
            break
    return r


def in_span(rows, k_target):
    """column k_target in span of the other columns? returns (bool, coefficients as Fractions or None).
    Exact: solve on a pivot subset, then verify on every row."""
    ncol = len(rows[0])
    others = [c for c in range(ncol) if c != k_target]
    rA = rank_exact([[r[c] for c in others] for r in rows])
    rAb = rank_exact(rows)
    if rAb > rA:
        return False, None
    # find independent rows
    sel = []
    for i in range(len(rows)):
        cand = sel + [i]
        if rank_exact([[rows[x][c] for c in others] for x in cand]) == len(cand):
            sel = cand
        if len(sel) == rA:
            break
    As = sp.Matrix([[rows[x][c] for c in others] for x in sel])
    bs = sp.Matrix([rows[x][k_target] for x in sel])
    if rA == len(others):
        sol = As.solve(bs)
    else:
        sol, params = As.gauss_jordan_solve(bs)
        sol = sol.subs({q: 0 for q in params})
    sol = [sp.Rational(x) for x in sol]
    for r in rows:
        if sum(s * r[c] for s, c in zip(sol, others)) != r[k_target]:
            raise AssertionError("verification failed")
    return True, sol
