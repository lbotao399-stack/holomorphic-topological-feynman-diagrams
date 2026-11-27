"""Quadratic Lie superalgebras g_tot for holomorphic CS/BF theories on C^H.

Conventions (checked by `check()`):
  basis T_A, parity p[A];  [T_A, T_B] = sum_C f[A,B,C] T_C;  w[A,B] = omega(T_A, T_B), parity p_omega
  [X,Y] = -(-1)^{|X||Y|}[Y,X];  super Jacobi;  omega(X,Y) = (-1)^{|X||Y|} omega(Y,X);  omega([X,Y],Z) = omega(X,[Y,Z])
  letters phi^A (coordinates of g_tot[1]) have parity pi[A] = p[A] + 1;  Phi = sum_A T_A phi^A is odd.
  Vertex  I = (1/6) omega(Phi, [Phi, Phi]) = (1/6) sum V[A,B,C] phi^A phi^B phi^C  (right-coefficient Koszul rule
  omega(X a, Y b) = (-1)^{|a||Y|} omega(X,Y) a b,  [X a, Y b] = (-1)^{|a||Y|} [X,Y] a b).
  Propagator numerators G = Gint / Gden  (Gint integer), fixed in `propagator()` by Q0^2 = 0.
All structure constants are small integers (numpy int64).

Constructors: matrix Lie algebras (so(N), sl(N), sp(2N)), reps; semidirect g |x PiR; L_U = g + PiU + g^v (6d hyper form);
tensor with Lambda[eps_1..eps_k] (odd) ; BF double T*_p L with pairing parity p.
"""
import itertools
from fractions import Fraction as Fr

import numpy as np


# ============================================================================ matrix Lie algebras
def _solve_coords(basis, M):
    """coordinates of matrix M in the span of basis (exact, Fractions)."""
    n = len(basis)
    rows = []
    flat = [b.flatten() for b in basis]
    target = M.flatten()
    # least squares exact: pick independent coordinates
    import sympy as sp
    A = sp.Matrix([[int(flat[j][k]) for j in range(n)] for k in range(len(target))])
    bvec = sp.Matrix([int(x) for x in target])
    sol, params = A.gauss_jordan_solve(bvec)
    if params.shape[0]:
        raise ValueError('basis not independent')
    return [Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in sol]


class MatrixLie:
    """Lie algebra spanned by integer matrices; f integer; kappa(X,Y) = scale * tr(XY)."""

    def __init__(self, name, mats, scale=1):
        self.name = name
        self.mats = [np.array(m, dtype=np.int64) for m in mats]
        n = self.dim = len(mats)
        f = np.zeros((n, n, n), dtype=object)
        for a in range(n):
            for b in range(n):
                C = self.mats[a] @ self.mats[b] - self.mats[b] @ self.mats[a]
                if not C.any():
                    continue
                co = _solve_coords(self.mats, C)
                for c in range(n):
                    f[a, b, c] = co[c]
        if any(x.denominator != 1 for x in f.flatten() if isinstance(x, Fr)):
            raise ValueError('non-integer structure constants')
        self.f = np.array([[[int(x) for x in r2] for r2 in r1] for r1 in f], dtype=np.int64)
        k = np.zeros((n, n), dtype=object)
        for a in range(n):
            for b in range(n):
                k[a, b] = Fr(scale) * int(np.trace(self.mats[a] @ self.mats[b]))
        if any(x.denominator != 1 for x in k.flatten()):
            raise ValueError('non-integer Killing form')
        self.kappa = np.array([[int(x) for x in r] for r in k], dtype=np.int64)


def E(N, i, j):
    m = np.zeros((N, N), dtype=np.int64)
    m[i, j] = 1
    return m


def so(N):
    mats = [E(N, i, j) - E(N, j, i) for i in range(N) for j in range(i + 1, N)]
    g = MatrixLie('so%d' % N, mats, scale=Fr(-1, 2))
    return g


def sl(N):
    mats = [E(N, i, j) for i in range(N) for j in range(N) if i != j]
    mats += [E(N, i, i) - E(N, i + 1, i + 1) for i in range(N - 1)]
    return MatrixLie('sl%d' % N, mats)


def gl(N):
    return MatrixLie('gl%d' % N, [E(N, i, j) for i in range(N) for j in range(N)])


def sp(N2):
    """sp(N2) = {X : X^T J + J X = 0}, J = [[0,1],[-1,0]] blocks; integer basis."""
    n = N2 // 2
    J = np.zeros((N2, N2), dtype=np.int64)
    J[:n, n:] = np.eye(n, dtype=np.int64)
    J[n:, :n] = -np.eye(n, dtype=np.int64)
    mats = []
    for i in range(n):
        for j in range(n):
            mats.append(E(N2, i, j) - E(N2, n + j, n + i))        # A block
    for i in range(n):
        for j in range(i, n):
            mats.append(E(N2, i, n + j) + (E(N2, j, n + i) if i != j else 0))   # B symmetric
            mats.append(E(N2, n + i, j) + (E(N2, n + j, i) if i != j else 0))   # C symmetric
    for m in mats:
        assert not (m.T @ J + J @ m).any()
    g = MatrixLie('sp%d' % N2, mats)
    g.J = J
    return g


def su2_eps():
    """su(2) with f = eps_abc, kappa = delta (adjoint-only theories)."""
    class G:
        pass
    g = G()
    g.name, g.dim = 'su2', 3
    f = np.zeros((3, 3, 3), dtype=np.int64)
    for a, b, c in itertools.permutations(range(3)):
        f[a, b, c] = 1 if (a, b, c) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1
    g.f, g.kappa = f, np.eye(3, dtype=np.int64)
    g.mats = None
    return g


# ============================================================================ superalgebra container
class SuperLie:
    def __init__(self, name, par, f, w=None, pomega=None, names=None, blocks=None):
        self.name = name
        self.par = np.array(par, dtype=np.int64) % 2
        self.N = len(self.par)
        self.f = np.array(f, dtype=np.int64)
        self.w = None if w is None else np.array(w, dtype=np.int64)
        self.pomega = pomega
        self.names = names or ['T%d' % i for i in range(self.N)]
        self.blocks = blocks or {}          # block name -> list of basis indices
        self.pi = (self.par + 1) % 2         # letter parities

    def sgn(self, x):
        return -1 if x % 2 else 1

    # ------------------------------------------------------------------ checks
    def check(self, quadratic=True):
        N, p, f = self.N, self.par, self.f
        # parity of brackets
        for A, B, C in zip(*np.nonzero(f)):
            assert (p[A] + p[B] - p[C]) % 2 == 0, 'bracket parity'
        # graded antisymmetry
        for A in range(N):
            for B in range(N):
                s = -self.sgn(p[A] * p[B])
                assert (f[A, B] == s * f[B, A]).all(), 'antisymmetry %d %d' % (A, B)
        # super Jacobi: [A,[B,C]] = [[A,B],C] + (-1)^{|A||B|}[B,[A,C]]
        fA = f.astype(object)
        lhs = np.einsum('bcd,ade->abce', fA, fA)                 # [A,[B,C]] -> e
        t1 = np.einsum('abd,dce->abce', fA, fA)                  # [[A,B],C]
        t2 = np.einsum('acd,bde->abce', fA, fA)                  # [B,[A,C]]
        sg = np.array([[self.sgn(p[a] * p[b]) for b in range(N)] for a in range(N)], dtype=object)
        rhs = t1 + sg[:, :, None, None] * t2
        assert (lhs == rhs).all(), 'Jacobi'
        if not quadratic:
            return True
        w = self.w
        assert w is not None
        for A, B in zip(*np.nonzero(w)):
            assert (p[A] + p[B] - self.pomega) % 2 == 0, 'pairing parity'
        for A in range(N):
            for B in range(N):
                assert w[A, B] == self.sgn(p[A] * p[B]) * w[B, A], 'pairing symmetry'
        # invariance omega([A,B],C) = omega(A,[B,C])
        L = np.einsum('abd,dc->abc', fA, w.astype(object))
        R = np.einsum('bcd,ad->abc', fA, w.astype(object))
        assert (L == R).all(), 'invariance'
        assert _rank_mod(w) == N, 'degenerate pairing'
        return True

    # ------------------------------------------------------------------ vertex and propagator
    def vertex(self):
        """V[A,B,C]: I = (1/6) sum V phi^A phi^B phi^C (graded symmetric; checked)."""
        p, pi, f, w = self.par, self.pi, self.f, self.w
        N = self.N
        fw = np.einsum('bcd,ad->abc', f, w)          # omega(T_A, [T_B, T_C])
        V = np.zeros((N, N, N), dtype=np.int64)
        for A, B, C in zip(*np.nonzero(fw)):
            s = self.sgn(pi[B] * p[C] + pi[A] * (p[B] + p[C]))
            V[A, B, C] = s * fw[A, B, C]
        # graded symmetry under letter swaps
        for A, B, C in zip(*np.nonzero(V)):
            assert V[B, A, C] == self.sgn(pi[A] * pi[B]) * V[A, B, C], 'vertex symmetry 12'
            assert V[A, C, B] == self.sgn(pi[B] * pi[C]) * V[A, B, C], 'vertex symmetry 23'
        self.V = V
        return V

    def propagator(self, variant=0):
        """G = Gint/Gden with Gint integer: inverse pairing, index convention by `variant`."""
        w = self.w
        winv = np.array(_inverse_fr(w), dtype=object)
        den = 1
        for x in winv.flatten():
            den = den * x.denominator // np.gcd(den, x.denominator)
        Gint = np.array([[int(x * den) for x in r] for r in winv], dtype=np.int64)
        if variant == 0:
            G = Gint.T.copy()                   # G[A,B] = (w^{-1})[B,A]
        else:
            G = Gint.copy()
        self.Gint, self.Gden = G, int(den)
        return G, int(den)


def _rank_mod(M, p=2147483647):
    A = [[int(x) % p for x in r] for r in M]
    rk, nr, nc = 0, len(A), len(A[0])
    for c in range(nc):
        piv = next((i for i in range(rk, nr) if A[i][c]), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        inv = pow(A[rk][c], p - 2, p)
        A[rk] = [x * inv % p for x in A[rk]]
        for i in range(nr):
            if i != rk and A[i][c]:
                fac = A[i][c]
                A[i] = [(x - fac * y) % p for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def _inverse_fr(w):
    import sympy as sp
    M = sp.Matrix(w.tolist()).inv()
    return [[Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in row] for row in M.tolist()]


# ============================================================================ constructions
def lie_even(g, with_pairing=True):
    """ordinary Lie algebra as an even superalgebra (pairing = kappa, parity 0)."""
    return SuperLie(g.name, [0] * g.dim, g.f, g.kappa if with_pairing else None, 0,
                    names=['%s%d' % ('X', a) for a in range(g.dim)], blocks={'g': list(range(g.dim))})


def semidirect(g, rho, rname='r', rpar=1):
    """L = g |x R (R abelian, parity rpar), rho: list over a of matrices, X_a r_i = sum_j rho_a[j,i] r_j."""
    n, d = g.dim, rho[0].shape[0]
    N = n + d
    f = np.zeros((N, N, N), dtype=np.int64)
    f[:n, :n, :n] = g.f
    for a in range(n):
        for i in range(d):
            for j in range(d):
                if rho[a][j, i]:
                    f[a, n + i, n + j] = rho[a][j, i]
                    f[n + i, a, n + j] = -rho[a][j, i]
    par = [0] * n + [rpar] * d
    names = ['X%d' % a for a in range(n)] + ['%s%d' % (rname, i) for i in range(d)]
    return SuperLie(g.name + '+' + rname, par, f, None, None, names,
                    blocks={'g': list(range(n)), rname: list(range(n, N))})


def hyper_LU(g, rho, Om, uname='u'):
    """L_U = g + PiU + g^v : [g,g], [g,U] = rho, [g,g^v] coadjoint, [U,U] -> g^v (moment map);
    omega(X_a, X^b) = delta, omega(u_i,u_j) = Om[i,j] (antisymmetric, g-invariant).  Even pairing."""
    n, d = g.dim, Om.shape[0]
    N = 2 * n + d
    iU = list(range(n, n + d))
    iD = list(range(n + d, N))
    f = np.zeros((N, N, N), dtype=np.int64)
    f[:n, :n, :n] = g.f
    for a in range(n):
        for i in range(d):
            for j in range(d):
                if rho[a][j, i]:
                    f[a, iU[i], iU[j]] = rho[a][j, i]
                    f[iU[i], a, iU[j]] = -rho[a][j, i]
        for b in range(n):
            for c in range(n):
                if g.f[c, a, b]:                      # [X_a, X^b] = sum_c f[c,a,b] X^c
                    f[a, iD[b], iD[c]] = g.f[c, a, b]
                    f[iD[b], a, iD[c]] = -g.f[c, a, b]
    for b in range(n):
        M = -(Om @ rho[b])                             # [u_i,u_j] = sum_b M_b[i,j] X^b
        assert (M == M.T).all()
        for i in range(d):
            for j in range(d):
                if M[i, j]:
                    f[iU[i], iU[j], iD[b]] = M[i, j]
    w = np.zeros((N, N), dtype=np.int64)
    for a in range(n):
        w[a, iD[a]] = 1
        w[iD[a], a] = 1
    w[n:n + d, n:n + d] = Om
    par = [0] * n + [1] * d + [0] * n
    names = ['X%d' % a for a in range(n)] + ['%s%d' % (uname, i) for i in range(d)] + ['Xv%d' % a for a in range(n)]
    return SuperLie(g.name + '+U', par, f, w, 0, names,
                    blocks={'g': list(range(n)), uname: iU, 'gv': iD})


def tensor_lambda(L, k, gen='e'):
    """L (x) Lambda[eps_1..eps_k] (odd generators), no pairing.  Basis (i, mask)."""
    masks = list(range(1 << k))
    def mpar(m):
        return bin(m).count('1') % 2
    def mul(m1, m2):
        if m1 & m2:
            return 0, 0
        # sign of reordering eps's of m1 m2 into increasing order
        s = 0
        for a in range(k):
            if m1 >> a & 1:
                s += bin(m2 & ((1 << a) - 1)).count('1')
        return (-1) ** s, m1 | m2
    idx = {(i, m): t for t, (m, i) in enumerate(itertools.product(masks, range(L.N)))}
    N = len(idx)
    par = [0] * N
    names = [''] * N
    for (i, m), t in idx.items():
        par[t] = (L.par[i] + mpar(m)) % 2
        names[t] = L.names[i] + ''.join('%s%d' % (gen, a + 1) for a in range(k) if m >> a & 1)
    f = np.zeros((N, N, N), dtype=np.int64)
    for (i, m), t in idx.items():
        for (j, m2), u in idx.items():
            s, mm = mul(m, m2)
            if s == 0:
                continue
            s *= (-1) ** (mpar(m) * L.par[j])
            for kk in np.nonzero(L.f[i, j])[0]:
                f[t, u, idx[(kk, mm)]] += s * L.f[i, j, kk]
    blocks = {}
    for bname, ids in L.blocks.items():
        for m in masks:
            tag = bname + ''.join('%s%d' % (gen, a + 1) for a in range(k) if m >> a & 1)
            blocks[tag] = [idx[(i, m)] for i in ids]
    return SuperLie(L.name + '(x)L%d' % k, par, f, None, None, names, blocks)


def double(L, pomega, suffix='^v'):
    """T*_p L = L + L^v, [e_i, e^j] = sum_l f[l,i,j] e^l, omega(e_i, e^j) = delta, parity(e^i) = parity(e_i) + p."""
    n = L.N
    N = 2 * n
    f = np.zeros((N, N, N), dtype=np.int64)
    f[:n, :n, :n] = L.f
    par = list(L.par) + [(x + pomega) % 2 for x in L.par]
    for i in range(n):
        for j in range(n):
            for l in range(n):
                if L.f[l, i, j]:
                    f[i, n + j, n + l] = L.f[l, i, j]
                    s = -(-1) ** (par[n + j] * par[i])
                    f[n + j, i, n + l] = s * L.f[l, i, j]
    w = np.zeros((N, N), dtype=np.int64)
    for i in range(n):
        w[i, n + i] = 1
        w[n + i, i] = (-1) ** (par[i] * par[n + i])
    names = list(L.names) + [x + suffix for x in L.names]
    blocks = dict(L.blocks)
    for b, ids in L.blocks.items():
        blocks[b + suffix] = [n + i for i in ids]
    return SuperLie('T*' + L.name, par, f, w, pomega, names, blocks)


def berezin_tensor(g, k):
    """g (x) Lambda[th_1..th_k] with pairing omega(x a, y b) = kappa(x,y) int a b (6d/LLC superfield form).
    Basis order: mask-major as in tensor_lambda; int th_1...th_k = 1."""
    L = lie_even(g, with_pairing=False)
    T = tensor_lambda(L, k, gen='t')
    n = g.dim
    top = (1 << k) - 1
    masks = list(range(1 << k))
    def mul(m1, m2):
        if m1 & m2:
            return 0, 0
        s = 0
        for a in range(k):
            if m1 >> a & 1:
                s += bin(m2 & ((1 << a) - 1)).count('1')
        return (-1) ** s, m1 | m2
    w = np.zeros((T.N, T.N), dtype=np.int64)
    for m1 in masks:
        for m2 in masks:
            s, mm = mul(m1, m2)
            if s and mm == top:
                for a in range(n):
                    for b in range(n):
                        if g.kappa[a, b]:
                            w[m1 * n + a, m2 * n + b] = s * g.kappa[a, b]
    T.w, T.pomega = w, k % 2
    T.name = '%s(x)Lambda%d' % (g.name, k)
    return T


# ============================================================================ representations
def fundamental(g):
    return [m.copy() for m in g.mats]


def adjoint_rep(g):
    # (ad X_a)[c, b] = f[a,b,c]
    return [g.f[a].T.copy() for a in range(g.dim)]


def dual_rep(rho):
    return [-m.T.copy() for m in rho]


def direct_sum(*reps):
    n = len(reps[0])
    out = []
    for a in range(n):
        blocks = [r[a] for r in reps]
        d = sum(b.shape[0] for b in blocks)
        M = np.zeros((d, d), dtype=np.int64)
        o = 0
        for b in blocks:
            M[o:o + b.shape[0], o:o + b.shape[0]] = b
            o += b.shape[0]
        out.append(M)
    return out


def cotangent_hyper(rho):
    """U = R + R^*, rho_U = rho (+) -rho^T, Om = [[0, 1], [-1, 0]]."""
    d = rho[0].shape[0]
    U = direct_sum(rho, dual_rep(rho))
    Om = np.zeros((2 * d, 2 * d), dtype=np.int64)
    Om[:d, d:] = np.eye(d, dtype=np.int64)
    Om[d:, :d] = -np.eye(d, dtype=np.int64)
    for M in U:
        assert not (Om @ M + M.T @ Om).any()
    return U, Om


if __name__ == '__main__':
    for g in (so(3), so(5), sl(2), sl(3), sp(4), su2_eps()):
        print(g.name, g.dim, 'f integer', g.f.dtype, 'kappa det', round(np.linalg.det(g.kappa.astype(float))))
