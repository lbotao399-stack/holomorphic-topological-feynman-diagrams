"""Own Lie-algebra layer for VERIFY.md (does not import the explorer's lie.py / cochain.py).

Alg(mats, form): basis matrices T_a (integer numpy arrays), invariant form kappa_ab = form(T_a, T_b)
(default tr(T_a T_b)).  fu[a,b,c]: [T_a, T_b] = sum_c fu[a,b,c] T_c (exact solve).
f[a,b,c] = kappa([T_a,T_b], T_c).  kinv = kappa^{-1} as Fractions; kinv_int = s * kinv integer.
Colour weights from an explicit list of vertex index-words and edge pairs (einsum, object or int64).
"""
from fractions import Fraction as Fr
import itertools
import numpy as np


def _solve_coords(basis_flat, target_flat):
    """basis_flat: dim x m integer, find exact c with c @ basis_flat = target (unique)."""
    dim, m = basis_flat.shape
    # pick dim independent columns greedily
    A = [[Fr(int(x)) for x in row] for row in basis_flat.T.tolist()]  # m x dim
    rows = []
    M = []
    for r in range(m):
        cand = M + [A[r]]
        if _rank(cand) > len(M):
            M = cand
            rows.append(r)
        if len(M) == dim:
            break
    assert len(M) == dim
    return rows


def _rank(A):
    A = [list(r) for r in A]
    if not A:
        return 0
    n, m = len(A), len(A[0])
    r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        for i in range(r + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
    return r


def inv_frac(A):
    n = len(A)
    M = [[Fr(x) for x in A[i]] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


class Alg:
    def __init__(self, name, mats, form=None):
        self.name = name
        self.mats = [np.array(m, dtype=np.int64) for m in mats]
        self.N = self.mats[0].shape[0]
        self.dim = d = len(mats)
        form = form or (lambda X, Y: int(np.trace(X @ Y)))
        self.form = form
        flat = np.array([m.reshape(-1) for m in self.mats], dtype=np.int64)  # d x N^2
        cols = _solve_coords(flat, None)
        Bsub = [[Fr(int(flat[a, c])) for c in cols] for a in range(d)]  # d x d
        Binv = inv_frac(Bsub)
        self._cols, self._Binv = cols, Binv
        fu = np.zeros((d, d, d), dtype=object)
        for a in range(d):
            for b in range(d):
                cm = self.mats[a] @ self.mats[b] - self.mats[b] @ self.mats[a]
                fu[a, b] = self.coords(cm)
        self.fu = fu
        assert all(Fr(x).denominator == 1 for x in fu.reshape(-1)), "non-integer structure constants"
        self.fu_int = np.array([[[int(x) for x in r2] for r2 in r1] for r1 in fu], dtype=np.int64)
        self.kappa = [[Fr(form(self.mats[a], self.mats[b])) for b in range(d)] for a in range(d)]
        self.kinv = inv_frac(self.kappa)
        den = 1
        for r in self.kinv:
            for x in r:
                den = den * x.denominator // np.gcd(den, x.denominator)
        self.kinv_scale = den
        self.kinv_int = np.array([[int(x * den) for x in r] for r in self.kinv], dtype=np.int64)
        kap = np.array([[int(x) for x in r] for r in self.kappa], dtype=np.int64)
        self.f = np.einsum('abe,ec->abc', self.fu_int, kap)
        # checks: invariance of kappa, antisymmetry of f, Jacobi
        assert np.array_equal(self.f, -np.transpose(self.f, (1, 0, 2)))
        assert np.array_equal(self.f, np.transpose(self.f, (1, 2, 0)))
        J = (np.einsum('abe,ecf->abcf', self.fu_int, self.fu_int) + np.einsum('bce,eaf->abcf', self.fu_int, self.fu_int)
             + np.einsum('cae,ebf->abcf', self.fu_int, self.fu_int))
        assert not J.any()

    def coords(self, m):
        """exact coordinates of a matrix in the span of the basis (asserts membership)."""
        flat = np.array(m).reshape(-1)
        c = [sum(Fr(int(flat[self._cols[j]])) * self._Binv[j][a] for j in range(self.dim)) for a in range(self.dim)]
        rec = sum((int(c[a]) if Fr(c[a]).denominator == 1 else 0) * self.mats[a] for a in range(self.dim))
        assert all(Fr(x).denominator == 1 for x in c) and np.array_equal(rec, np.array(m)), "not in span / non-integer"
        return [int(x) for x in c]

    def mat(self, c):
        return sum(int(c[a]) * self.mats[a] for a in range(self.dim))

    def bracket(self, c1, c2):
        return [int(x) for x in np.einsum('a,b,abc->c', np.array(c1, dtype=np.int64), np.array(c2, dtype=np.int64), self.fu_int)]

    def center_dim(self):
        M = self.fu_int.transpose(0, 2, 1).reshape(self.dim, -1)
        return self.dim - _rank([[Fr(int(x)) for x in r] for r in M.tolist()])


def E(N, i, j):
    m = np.zeros((N, N), dtype=np.int64)
    m[i, j] = 1
    return m


def gl(N, twist=None):
    """gl_N; twist: optional integer unimodular matrix U (basis T_a = sum_b U_ab E_b)."""
    base = [E(N, i, j) for i in range(N) for j in range(N)]
    if twist is not None:
        base = [sum(int(twist[a][b]) * base[b] for b in range(len(base))) for a in range(len(base))]
    return Alg("gl%d" % N, base)


def sl(N):
    """sl_N with basis E_ij (i != j) and E_kk - E_NN (differs from the explorer's H_k)."""
    base = [E(N, i, j) for i in range(N) for j in range(N) if i != j]
    base += [E(N, k, k) - E(N, N - 1, N - 1) for k in range(N - 1)]
    return Alg("sl%d" % N, base)


def so(n):
    base = [E(n, i, j) - E(n, j, i) for i in range(n) for j in range(i + 1, n)]
    return Alg("so%d" % n, base)


def block_sum(algs_mats, sizes, name, form):
    """block-diagonal embedding of several matrix algebras (list of lists of matrices)."""
    Ntot = sum(sizes)
    mats = []
    off = 0
    for ms, s in zip(algs_mats, sizes):
        for m in ms:
            M = np.zeros((Ntot, Ntot), dtype=np.int64)
            M[off:off + s, off:off + s] = m
            mats.append(M)
        off += s
    return Alg(name, mats, form)


def unimodular(d, rng, steps=None):
    U = [[int(i == j) for j in range(d)] for i in range(d)]
    for _ in range(steps or 2 * d):
        i, j = rng.sample(range(d), 2)
        s = rng.choice([-1, 1])
        U[i] = [x + s * y for x, y in zip(U[i], U[j])]
    return U


# ---------------------------------------------------------------- colour weights
def contract(alg, vertices, edges_pairs, ext_letters, dtype=object):
    """vertices: list of 3-letter words (f_{word}); edges_pairs: list of (l, l') contracted with
    kinv_int; ext_letters: output order. Returns integer tensor = kinv_scale^{#edges} * W."""
    ops, terms = [], []
    f = alg.f.astype(dtype)
    ki = alg.kinv_int.astype(dtype)
    for w in vertices:
        terms.append(w)
        ops.append(f)
    for (l, lp) in edges_pairs:
        terms.append(l + lp)
        ops.append(ki)
    return np.einsum(','.join(terms) + '->' + ''.join(ext_letters), *ops, optimize='greedy'), alg.kinv_scale ** len(edges_pairs)


def W133_index(alg, dtype=object):
    # W_133 = - f_xyz f_{a1 y'u} f_{a2 u'w} f_{a3 z's} f_{a4 s'r} f_{x'w'r'}  (c3)
    # letters: x y z u w s r ; primes X Y Z U W S R ; externals a b c d
    V = ['xyz', 'aYu', 'bUw', 'cZs', 'dSr', 'XWR']
    Ep = [('x', 'X'), ('y', 'Y'), ('z', 'Z'), ('u', 'U'), ('w', 'W'), ('s', 'S'), ('r', 'R')]
    T, sc = contract(alg, V, Ep, 'abcd', dtype)
    return -T, sc


def W223_index(alg, dtype=object):
    # W_223 = - f_xyz f_{a1 x'u} f_{a2 y'w} f_{a3 z's} f_{a4 s'r} f_{u'w'r'}
    V = ['xyz', 'aXu', 'bYw', 'cZs', 'dSr', 'UWR']
    Ep = [('x', 'X'), ('y', 'Y'), ('z', 'Z'), ('u', 'U'), ('w', 'W'), ('s', 'S'), ('r', 'R')]
    T, sc = contract(alg, V, Ep, 'abcd', dtype)
    return -T, sc


def orientation_rule(alg, nv, edges, ext, dtype=object):
    """W = sgn(L -> L') * contraction (rule stated in the explorer's c3 / cochain.py docstring),
    implemented here from the statement."""
    L = [('x', v) for v in ext]
    for e in range(len(edges)):
        L += [('e', e, 0), ('e', e, 1)]
    Lp = []
    words = []
    for v in range(nv):
        hs = [('x', v)] if v in ext else []
        hs += sorted([('e', e, s) for e, ab in enumerate(edges) for s in (0, 1) if ab[s] == v], key=lambda h: (h[1], h[2]))
        assert len(hs) == 3
        Lp += hs
        words.append(hs)
    pos = {h: i for i, h in enumerate(L)}
    perm = [pos[h] for h in Lp]
    inv = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
    sgn = -1 if inv % 2 else 1
    letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    let = {h: letters[i] for i, h in enumerate(L)}
    V = [''.join(let[h] for h in hs) for hs in words]
    Ep = [(let[('e', e, 0)], let[('e', e, 1)]) for e in range(len(edges))]
    T, sc = contract(alg, V, Ep, ''.join(let[('x', v)] for v in ext), dtype)
    return sgn * T, sc, sgn


def perm_sign(p):
    p = list(p)
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def alt_part(W):
    k = W.ndim
    out = None
    for s in itertools.permutations(range(k)):
        term = perm_sign(s) * np.transpose(W, s)
        out = term if out is None else out + term
    return out
