"""Verifier's own Lie algebra constructions (independent of the explorer's algebras.py).

Data: Alg.n, Alg.c[a,b,d] (int64):  [e_a, e_b] = sum_d c[a,b,d] e_d.
Killing form  K[a,b] = tr(ad e_a ad e_b) = sum_{d,e} c[a,d,e] c[b,e,d]   (int, exact).
Lowered       fK[a,b,c] = K([e_a,e_b], e_c) = sum_d c[a,b,d] K[d,c].
Constructions:
  so(N):   e = E_ij - E_ji (i<j), N x N matrices.
  sp(2m):  e = J S, S in {E_ii, E_ij + E_ji}, J = [[0,1],[-1,0]] (m x m blocks)   (X^T J + J X = 0).
  sl(N):   e = E_ij (i != j), E_ii - E_NN (i < N).
  G2:      Der(O) on the Cayley-Dickson octonions O = H + H l, nullspace of the derivation system (exact).
  F4:      so(9) + S_16 with octonionic real gamma matrices on O + O:
           g_i = [[0, L_i^T], [L_i, 0]] (i = 0..7, L_i = left mult. by e_i), g_8 = diag(1_8, -1_8);
           basis Y_ij = 2(E_ij - E_ji) (i<j), s_alpha;  [Y_ij, s] = g_i g_j s,
           [s, t] = sum_{i<j} (s^T g_i g_j t) Y_ij.  Jacobi is checked exactly, not assumed.
"""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr


# ------------------------------------------------------------------------------------------- helpers
def coords_solver(mats):
    """exact coordinate map for a list of integer matrices (linearly independent)."""
    B = np.array([m.reshape(-1) for m in mats], dtype=object)       # n x N^2
    n = B.shape[0]
    # positions: choose columns of B (entries) giving an invertible n x n minor
    Mb = sp.Matrix(B.tolist())
    _, pivpos = Mb.rref()
    pivpos = list(pivpos)
    assert len(pivpos) == n, "basis matrices not independent"
    S = Mb[:, pivpos]               # n x n
    Sinv = S.inv()

    def coords(m):
        v = sp.Matrix([[m.reshape(-1)[p] for p in pivpos]])
        x = v * Sinv                # 1 x n
        x = [sp.Rational(t) for t in x]
        # verify
        R = np.zeros(m.size, dtype=object)
        for a in range(n):
            if x[a] != 0:
                R = R + B[a] * x[a]
        assert all(sp.Rational(R[i]) == int(m.reshape(-1)[i]) for i in range(m.size)), "not in span"
        return x
    return coords


class Alg:
    def __init__(self, name, c, info=""):
        self.name = name
        self.c = np.asarray(c, dtype=np.int64)
        self.n = self.c.shape[0]
        self.info = info
        self._derived()

    def _derived(self):
        c = self.c
        n = self.n
        self.ad = np.transpose(c, (0, 2, 1)).copy()      # ad[a][d,b] = c[a,b,d]  (matrix of ad e_a)
        K = np.einsum('ade,bed->ab', c, c)
        self.K = K
        self.fK = np.einsum('abd,dc->abc', c, K)
        Ks = sp.Matrix(K.tolist())
        self.detK = Ks.det()
        inv = Ks.inv()
        s = int(sp.ilcm(*[sp.Rational(x).q for x in inv]))
        self.s = s
        self.Gi = np.array((inv * s).tolist(), dtype=np.int64)   # K^{-1} = Gi / s
        assert np.array_equal(self.Gi @ K, s * np.eye(n, dtype=np.int64))

    def checks(self):
        c, n = self.c, self.n
        out = {}
        out['antisym'] = bool(np.array_equal(c, -c.transpose(1, 0, 2)))
        ok = True
        for a in range(n):
            # [[e_a,e_b],e_c] + [[e_b,e_c],e_a] + [[e_c,e_a],e_b] = 0
            t1 = np.einsum('bd,dce->bce', c[a], c)
            t2 = np.einsum('bcd,de->bce', c, c[:, a, :])
            t3 = np.einsum('cd,dbe->bce', c[:, a, :], c)
            if (t1 + t2 + t3).any():
                ok = False
                break
        out['jacobi'] = ok
        out['K_invariant'] = bool(np.array_equal(self.fK, -self.fK.transpose(0, 2, 1)))
        out['detK!=0'] = self.detK != 0
        return out

    def bracket_vec(self, x, y):
        """[x, y] for coordinate vectors (Python ints / object)."""
        return np.einsum('a,b,abd->d', x, y, self.c.astype(object))


def from_matrices(name, mats, info=""):
    coords = coords_solver(mats)
    n = len(mats)
    c = [[[None] * n for _ in range(n)] for _ in range(n)]
    den = 1
    cq = np.zeros((n, n, n), dtype=object)
    for a in range(n):
        for b in range(n):
            m = mats[a] @ mats[b] - mats[b] @ mats[a]
            x = coords(m)
            for d in range(n):
                cq[a, b, d] = x[d]
                den = sp.ilcm(den, sp.Rational(x[d]).q)
    den = int(den)
    # rescale basis by den so that structure constants are integers: [den X, den Y] = den (den [X,Y])
    ci = np.array([[[int(sp.Rational(cq[a, b, d]) * den) for d in range(n)] for b in range(n)]
                   for a in range(n)], dtype=np.int64)
    A = Alg(name, ci, info + ("; basis rescaled by %d" % den if den != 1 else ""))
    A.mats = [m * den for m in mats]
    return A


def E(N, i, j):
    m = np.zeros((N, N), dtype=np.int64)
    m[i, j] = 1
    return m


def so(N):
    mats = [E(N, i, j) - E(N, j, i) for i in range(N) for j in range(i + 1, N)]
    return from_matrices("so%d" % N, mats, "E_ij - E_ji in gl(%d)" % N)


def sp_(m):
    N = 2 * m
    J = np.zeros((N, N), dtype=np.int64)
    J[:m, m:] = np.eye(m, dtype=np.int64)
    J[m:, :m] = -np.eye(m, dtype=np.int64)
    mats = []
    for i in range(N):
        for j in range(i, N):
            S = E(N, i, j) + E(N, j, i) if i != j else E(N, i, i)
            X = J @ S
            assert not (X.T @ J + J @ X).any()
            mats.append(X)
    return from_matrices("sp%d" % N, mats, "J S, S symmetric, in gl(%d)" % N)


def sl(N):
    mats = [E(N, i, j) for i in range(N) for j in range(N) if i != j]
    mats += [E(N, i, i) - E(N, N - 1, N - 1) for i in range(N - 1)]
    return from_matrices("sl%d" % N, mats, "E_ij, E_ii - E_NN in gl(%d)" % N)


# ------------------------------------------------------------------------------------------- octonions
def quat_table():
    # basis 1, i, j, k ; table[a][b] = (sign, index)
    t = {}
    names = ['1', 'i', 'j', 'k']
    mult = {('i', 'j'): (1, 'k'), ('j', 'k'): (1, 'i'), ('k', 'i'): (1, 'j'),
            ('j', 'i'): (-1, 'k'), ('k', 'j'): (-1, 'i'), ('i', 'k'): (-1, 'j'),
            ('i', 'i'): (-1, '1'), ('j', 'j'): (-1, '1'), ('k', 'k'): (-1, '1')}
    M = np.zeros((4, 4, 4), dtype=np.int64)
    for a in range(4):
        for b in range(4):
            if a == 0:
                M[a, b, b] = 1
            elif b == 0:
                M[a, b, a] = 1
            else:
                s, r = mult[(names[a], names[b])]
                M[a, b, names.index(r)] = s
    return M


def octonion_table():
    """Cayley-Dickson: (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c)); O[x,y,z] = coefficient of e_z in e_x e_y."""
    Q = quat_table()
    conj = np.diag([1, -1, -1, -1])

    def qm(x, y):
        return np.einsum('a,b,abc->c', x, y, Q)

    O = np.zeros((8, 8, 8), dtype=np.int64)
    for x in range(8):
        for y in range(8):
            ex = np.zeros(8, dtype=np.int64); ex[x] = 1
            ey = np.zeros(8, dtype=np.int64); ey[y] = 1
            a, b = ex[:4], ex[4:]
            c, d = ey[:4], ey[4:]
            first = qm(a, c) - qm(conj @ d, b)
            second = qm(d, a) + qm(b, conj @ c)
            O[x, y] = np.concatenate([first, second])
    return O


def octonion_checks(O):
    """alternativity (x x) y = x (x y), y (x x) = (y x) x on all basis pairs and random integer pairs;
    norm multiplicativity |xy|^2 = |x|^2 |y|^2 at random integer points."""
    rng = np.random.default_rng(1)

    def m(x, y):
        return np.einsum('a,b,abc->c', x, y, O)
    ok = True
    for _ in range(30):
        x = rng.integers(-5, 6, 8); y = rng.integers(-5, 6, 8)
        ok &= np.array_equal(m(m(x, x), y), m(x, m(x, y)))
        ok &= np.array_equal(m(y, m(x, x)), m(m(y, x), x))
        ok &= int(m(x, y) @ m(x, y)) == int(x @ x) * int(y @ y)
    return bool(ok)


def g2():
    O = octonion_table()
    # D in gl(8) with D(e_x e_y) = D(e_x) e_y + e_x D(e_y);  D[r, x] = coefficient of e_r in D(e_x)
    n = 8
    rows = []
    for x in range(n):
        for y in range(n):
            for z in range(n):
                # component z of D(e_x e_y) - D(e_x) e_y - e_x D(e_y)
                row = [0] * 64
                for w in range(n):
                    if O[x, y, w]:
                        row[z * n + w] += int(O[x, y, w])          # D[z, w] * O[x,y,w]
                for r in range(n):
                    if O[r, y, z]:
                        row[r * n + x] -= int(O[r, y, z])          # D[r, x] O[r, y, z]
                    if O[x, r, z]:
                        row[r * n + y] -= int(O[x, r, z])          # D[r, y] O[x, r, z]
                if any(row):
                    rows.append(row)
    M = sp.Matrix(rows)
    ns = M.nullspace()
    mats = []
    for v in ns:
        den = sp.ilcm(*[sp.Rational(t).q for t in v])
        w = [int(sp.Rational(t) * den) for t in v]
        g = 0
        import math
        for t in w:
            g = math.gcd(g, abs(t))
        w = [t // g for t in w]
        mats.append(np.array(w, dtype=np.int64).reshape(8, 8))
    A = from_matrices("G2", mats, "Der(O), Cayley-Dickson octonions, %d x %d matrices on O" % (8, 8))
    A.octonion_ok = octonion_checks(O)
    A.n_der = len(ns)
    return A


def gammas9():
    O = octonion_table()
    L = []
    for i in range(8):
        # L_i[z, y] = coefficient of e_z in e_i e_y
        L.append(O[i].T.copy())
    g = []
    Z = np.zeros((8, 8), dtype=np.int64)
    I8 = np.eye(8, dtype=np.int64)
    for i in range(8):
        g.append(np.block([[Z, L[i].T], [L[i], Z]]))
    g.append(np.block([[I8, Z], [Z, -I8]]))
    return g


def f4():
    g = gammas9()
    I16 = np.eye(16, dtype=np.int64)
    cliff = all(np.array_equal(g[i] @ g[j] + g[j] @ g[i], 2 * (i == j) * I16) for i in range(9) for j in range(9))
    sym = all(np.array_equal(x, x.T) for x in g)
    pairs = [(i, j) for i in range(9) for j in range(i + 1, 9)]
    n = 36 + 16
    c = np.zeros((n, n, n), dtype=np.int64)
    # so(9) part: Y_ij = 2 (E_ij - E_ji); [Y_ij, Y_kl] = 2 * (Y-coordinates of commutator of (E_ij - E_ji))
    mats = [E(9, i, j) - E(9, j, i) for (i, j) in pairs]
    pos = {p: a for a, p in enumerate(pairs)}
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            m = mats[a] @ mats[b] - mats[b] @ mats[a]
            for (p, q), d in pos.items():
                if m[p, q]:
                    c[a, b, d] += 2 * int(m[p, q])
    rho = [g[i] @ g[j] for (i, j) in pairs]          # rho(Y_ij) = g_i g_j  (= 2 * (1/2) g_i g_j)
    for a in range(36):
        for al in range(16):
            v = rho[a][:, al]
            for be in range(16):
                if v[be]:
                    c[a, 36 + al, 36 + be] += int(v[be])
                    c[36 + al, a, 36 + be] -= int(v[be])
    for al in range(16):
        for be in range(16):
            for a in range(36):
                val = int(rho[a][al, be])     # s_al^T g_i g_j s_be
                if val:
                    c[36 + al, 36 + be, a] += val
    A = Alg("F4", c, "so(9) + S_16, octonionic gamma matrices; Y_ij = 2(E_ij - E_ji)")
    A.clifford_ok = cliff and sym
    # rho is a representation: [rho(Y_a), rho(Y_b)] = sum_d c[a,b,d] rho(Y_d)
    A.rho_ok = all(np.array_equal(rho[a] @ rho[b] - rho[b] @ rho[a],
                                  sum(c[a, b, d] * rho[d] for d in range(36))) for a in range(36) for b in range(36))
    return A


# ------------------------------------------------------------------------------------------- structure
def rank_mod_p(M, p):
    """rank of an integer matrix modulo prime p (dense Gaussian elimination, int64)."""
    A = np.array(M, dtype=np.int64) % p
    r = 0
    rows, cols = A.shape
    for col in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, col]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, col]), p - 2, p)
        A[r] = (A[r] * inv) % p
        nz = np.nonzero(A[:, col])[0]
        for i in nz:
            if i != r:
                A[i] = (A[i] - A[i, col] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r


def exact_rank(M):
    """exact rank over Q of an integer matrix (sympy DomainMatrix over ZZ -> QQ)."""
    from sympy.polys.matrices import DomainMatrix
    from sympy import ZZ
    rows = [[int(x) for x in r] for r in M]
    if not rows:
        return 0
    D = DomainMatrix([[ZZ(x) for x in r] for r in rows], (len(rows), len(rows[0])), ZZ)
    return D.convert_to(sp.QQ).rank()


def invariant_forms_nullity(A, p):
    """dim of {B symmetric bilinear: B([a,b],c) + B(b,[a,c]) = 0} modulo p (upper bound for Q-dim)."""
    n = A.n
    idx = {}
    for i in range(n):
        for j in range(i, n):
            idx[(i, j)] = len(idx)
    rows = []
    c = A.c
    for a in range(n):
        for b in range(n):
            for cc in range(b, n):
                row = {}
                for d in np.nonzero(c[a, b])[0]:
                    k = idx[tuple(sorted((d, cc)))]
                    row[k] = row.get(k, 0) + int(c[a, b, d])
                for d in np.nonzero(c[a, cc])[0]:
                    k = idx[tuple(sorted((b, d)))]
                    row[k] = row.get(k, 0) + int(c[a, cc, d])
                row = {k: v for k, v in row.items() if v}
                if row:
                    rows.append(row)
    M = np.zeros((len(rows), len(idx)), dtype=np.int64)
    for r, row in enumerate(rows):
        for k, v in row.items():
            M[r, k] = v
    # remove duplicate rows
    M = np.unique(M, axis=0)
    return len(idx) - rank_mod_p(M, p)


def rank_of_algebra(A, seed=0):
    """dim ker ad_X for a random integer X (>= rank; = rank for regular X), exact."""
    rng = np.random.default_rng(seed)
    x = rng.integers(-7, 8, A.n)
    adx = np.einsum('a,adb->db', x, A.ad)
    return A.n - exact_rank(adx.tolist())


def build(name):
    if name.startswith('so'):
        return so(int(name[2:]))
    if name.startswith('sp'):
        return sp_(int(name[2:]) // 2)
    if name.startswith('sl'):
        return sl(int(name[2:]))
    if name == 'G2':
        return g2()
    if name == 'F4':
        return f4()
    raise ValueError(name)
