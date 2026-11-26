"""Simple Lie algebras with exact integer (or rational, rescaled) structure constants.

Data of an algebra L:
  L.dim, L.fu[a,b,c]  : [b_a, b_b] = sum_c fu[a,b,c] b_c   (int64)
  L.kappa[a,b]        : invariant form (int64), L.form = description
  L.cartan            : basis indices spanning a Cartan subalgebra (diagonal action checked) or None
  L.wt[a]             : weight of b_a (tuple of ints, eigenvalues of ad(b_h), h in L.cartan) or None
  L.theta             : basis index of a highest root vector (generic functional), or 0

Constructions
  kac(A, M)           : simply-laced, Kac 'Infinite dimensional Lie algebras' Sec. 7.8 construction:
                        g = h + sum_alpha C E_alpha, [h, E_a] = (h|a) E_a, [E_a, E_-a] = -a,
                        [E_a, E_b] = eps(a,b) E_{a+b} (a+b root), eps(a,b) = (-1)^{a^T M b},
                        M_ii = 1, M_ij = 1 for an oriented edge i -> j; form (h|h') = A, (E_a|E_b) = -delta_{a,-b}.
                        Jacobi and invariance are checked exactly, not assumed.
  fold(L, sigma)      : fixed subalgebra of the diagram automorphism E_alpha -> E_{sigma alpha}
                        (an automorphism iff fu is invariant under the induced basis permutation: checked);
                        basis = orbit sums, form = restriction.
  so_matrix(N), sp_matrix(n) : X^T J + J X = 0, J antidiagonal (so) / antidiagonal with signs (sp),
                        basis E_ij -+ E_{j'i'}, form tr_V(XY) in the defining representation.
  g2_octonion()       : stabilizer of the 3-form phi in so(7) (orthonormal basis), form tr_7(XY).
"""
import itertools
import numpy as np
import sympy as sp
from fractions import Fraction as Fr


def perm_sign(seq):
    seq = list(seq); s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


class LA:
    def __init__(self, name, fu, kappa, form, cartan=None):
        self.name = name
        self.fu = np.asarray(fu, dtype=np.int64)
        self.kappa = np.asarray(kappa, dtype=np.int64)
        self.dim = self.fu.shape[0]
        self.form = form
        self.cartan = cartan
        self.wt = None
        self.theta = 0
        self._check()
        if cartan is not None:
            self._weights()

    # exact checks: antisymmetry, Jacobi, invariance and nondegeneracy of kappa
    def _check(self):
        fu, k, d = self.fu, self.kappa, self.dim
        assert np.array_equal(fu, -fu.transpose(1, 0, 2)), "antisymmetry"
        for a in range(d):
            # J[b,c,f] = [[a,b],c] + [[b,c],a] + [[c,a],b]
            t1 = np.einsum('be,ecf->bcf', fu[a], fu)
            # [[b,c],a] = sum_e fu[b,c,e] fu[e,a,f]
            t2 = np.einsum('bce,ef->bcf', fu, fu[:, a, :])
            # [[c,a],b] = sum_e fu[c,a,e] fu[e,b,f]
            t3 = np.einsum('ce,ebf->bcf', fu[:, a, :], fu)
            assert not (t1 + t2 + t3).any(), "Jacobi fails"
        assert np.array_equal(k, k.T), "kappa symmetric"
        fl = np.einsum('abe,ec->abc', fu, k)
        assert np.array_equal(fl, -fl.transpose(0, 2, 1)), "kappa invariant"
        self.fl = fl
        M = sp.Matrix(k.tolist())
        det = M.det()
        assert det != 0, "kappa degenerate"
        inv = M.inv()
        den = sp.ilcm(*[sp.Rational(x).q for x in inv])
        self.s = int(den)
        self.G = np.array((inv * den).tolist(), dtype=np.int64)       # G = s * kappa^{-1}
        assert np.array_equal(self.G @ k, self.s * np.eye(self.dim, dtype=np.int64))

    def _weights(self):
        d = self.dim
        wt = []
        for a in range(d):
            w = []
            for h in self.cartan:
                row = self.fu[h, a]
                nz = np.nonzero(row)[0]
                assert all(c == a for c in nz), "basis element %d is not a weight vector" % a
                w.append(int(row[a]) if len(nz) else 0)
            wt.append(tuple(w))
        self.wt = wt
        gen = np.array([1000003 ** i for i in range(len(self.cartan))][::-1], dtype=object)
        val = [sum(int(x) * int(g) for x, g in zip(w, gen)) for w in wt]
        self.theta = int(np.argmax(np.array(val, dtype=object)))

    def bracket_vec(self, x, y):
        return np.einsum('a,b,abc->c', x, y, self.fu)


# ----------------------------------------------------------------------------- Cartan matrices
def cartan_A(n):
    A = 2 * np.eye(n, dtype=int)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = -1
    return A


def cartan_D(n):
    A = 2 * np.eye(n, dtype=int)
    for i in range(n - 2):
        A[i, i + 1] = A[i + 1, i] = -1
    A[n - 3, n - 1] = A[n - 1, n - 3] = -1
    return A


def cartan_E6():
    # Bourbaki: 1-3-4-5-6, 2-4 (0-based: 0-2-3-4-5, 1-3)
    A = 2 * np.eye(6, dtype=int)
    for i, j in [(0, 2), (2, 3), (3, 4), (4, 5), (1, 3)]:
        A[i, j] = A[j, i] = -1
    return A


def orient(A, arrows):
    """M_ii = 1, M_ij = 1 for i -> j; arrows must cover every edge exactly once."""
    n = len(A)
    M = np.eye(n, dtype=int)
    edges = {(min(i, j), max(i, j)) for i in range(n) for j in range(n) if i != j and A[i, j]}
    got = set()
    for i, j in arrows:
        assert A[i, j] == -1
        M[i, j] = 1
        got.add((min(i, j), max(i, j)))
    assert got == edges
    return M


def positive_roots_simply_laced(A):
    n = len(A)
    simple = [tuple(int(i == j) for j in range(n)) for i in range(n)]
    roots = set(simple)
    frontier = list(simple)
    while frontier:
        new = []
        for r in frontier:
            for i in range(n):
                s = list(r); s[i] += 1
                v = np.array(s)
                if v @ A @ v == 2 and tuple(s) not in roots:
                    roots.add(tuple(s)); new.append(tuple(s))
        frontier = new
    return sorted(roots, key=lambda r: (sum(r), r))


def kac(name, A, M):
    A = np.array(A); n = len(A)
    pos = positive_roots_simply_laced(A)
    roots = pos + [tuple(-x for x in r) for r in pos]
    idx = {}
    # basis: h_1..h_n, then E_alpha
    for k, r in enumerate(roots):
        idx[r] = n + k
    d = n + len(roots)
    fu = np.zeros((d, d, d), dtype=np.int64)
    kap = np.zeros((d, d), dtype=np.int64)
    kap[:n, :n] = A
    for r in roots:
        a = idx[r]
        rv = np.array(r)
        Ar = A @ rv
        for i in range(n):
            fu[i, a, a] = Ar[i]
            fu[a, i, a] = -Ar[i]
        neg = tuple(-x for x in r)
        kap[a, idx[neg]] = -1
        for i in range(n):
            fu[a, idx[neg], i] = -rv[i]
        for s in roots:
            t = tuple(x + y for x, y in zip(r, s))
            if t in idx:
                e = (-1) ** int(rv @ M @ np.array(s))
                fu[a, idx[s], idx[t]] = e
    L = LA(name, fu, kap, "Kac form (h|h') = Cartan matrix, (E_a|E_-a) = -1", cartan=list(range(n)))
    L.roots = roots
    L.ridx = idx
    L.rank_ = n
    return L


def fold(name, L, sigma):
    """sigma: permutation of the simple roots (list); returns the fixed subalgebra."""
    n = L.rank_
    def act(r):
        out = [0] * n
        for i in range(n):
            out[sigma[i]] = r[i]
        return tuple(out)
    perm = list(range(L.dim))
    for i in range(n):
        perm[i] = sigma[i]
    for r in L.roots:
        perm[L.ridx[r]] = L.ridx[act(r)]
    P = np.array(perm)
    # automorphism check: fu[P a, P b, P c] = fu[a, b, c]
    assert np.array_equal(L.fu[np.ix_(P, P, P)], L.fu), "not an automorphism"
    assert np.array_equal(L.kappa[np.ix_(P, P)], L.kappa)
    # orbits
    seen, orbits = set(), []
    for a in range(L.dim):
        if a in seen:
            continue
        o, b = [], a
        while b not in o:
            o.append(b); b = perm[b]
        seen |= set(o)
        orbits.append(sorted(o))
    cart = [k for k, o in enumerate(orbits) if o[0] < n]
    m = len(orbits)
    B = np.zeros((m, L.dim), dtype=np.int64)
    for k, o in enumerate(orbits):
        B[k, o] = 1
    rep = [o[0] for o in orbits]
    fu = np.zeros((m, m, m), dtype=np.int64)
    for I in range(m):
        for J in range(m):
            v = np.einsum('a,b,abc->c', B[I], B[J], L.fu)
            c = v[rep]
            assert np.array_equal(c @ B, v), "bracket not in fixed subalgebra"
            fu[I, J] = c
    kap = B @ L.kappa @ B.T
    F = LA(name, fu, kap, "restriction of the %s Kac form to the fixed subalgebra" % L.name, cartan=cart)
    F.orbits = orbits
    return F


# ----------------------------------------------------------------------------- matrix algebras
def _from_matrices(name, mats, reps, form):
    """mats: basis matrices; reps[k] = (i, j) entry with mats[k][i,j] = 1 and mats[l][i,j] = 0 (l != k)."""
    d = len(mats)
    R = np.array([[int(m[i, j]) for (i, j) in reps] for m in mats])
    assert np.array_equal(R, np.eye(d, dtype=int)), "coordinate entries"
    fu = np.zeros((d, d, d), dtype=np.int64)
    for a in range(d):
        for b in range(d):
            c = mats[a] @ mats[b] - mats[b] @ mats[a]
            co = np.array([c[i, j] for (i, j) in reps], dtype=np.int64)
            assert np.array_equal(sum(co[k] * mats[k] for k in range(d)), c), "not closed"
            fu[a, b] = co
    kap = np.array([[int(np.trace(mats[a] @ mats[b])) for b in range(d)] for a in range(d)], dtype=np.int64)
    cart = [k for k in range(d) if np.array_equal(mats[k], np.diag(np.diag(mats[k])))]
    L = LA(name, fu, kap, form, cartan=cart)
    L.mats = mats
    return L


def E(N, i, j):
    m = np.zeros((N, N), dtype=np.int64); m[i, j] = 1
    return m


def so_matrix(N):
    pr = lambda k: N - 1 - k
    mats, reps = [], []
    for i in range(N):
        for j in range(N):
            if i + j < N - 1:
                mats.append(E(N, i, j) - E(N, pr(j), pr(i))); reps.append((i, j))
    assert len(mats) == N * (N - 1) // 2
    return _from_matrices("so%d" % N, mats, reps, "tr_%d(XY), defining (vector) representation of so(%d)" % (N, N))


def sp_matrix(n):
    N = 2 * n
    pr = lambda k: N - 1 - k
    s = lambda k: 1 if k < n else -1
    mats, reps = [], []
    for i in range(N):
        for j in range(N):
            if i + j < N - 1:
                mats.append(E(N, i, j) + s(pr(j)) * s(i) * E(N, pr(j), pr(i))); reps.append((i, j))
            elif i + j == N - 1:
                mats.append(E(N, i, j)); reps.append((i, j))
    J = sum(s(k) * E(N, k, pr(k)) for k in range(N))
    for m in mats:
        assert not (m.T @ J + J @ m).any()
    assert len(mats) == n * (2 * n + 1)
    return _from_matrices("sp%d" % N, mats, reps, "tr_%d(XY), defining representation of sp(%d)" % (N, N))


def g2_octonion():
    """G2 = {X in so(7) : X.phi = 0}, phi = e123 + e145 + e167 + e246 - e257 - e347 - e356."""
    terms = [((0, 1, 2), 1), ((0, 3, 4), 1), ((0, 5, 6), 1), ((1, 3, 5), 1), ((1, 4, 6), -1),
             ((2, 3, 6), -1), ((2, 4, 5), -1)]
    phi = np.zeros((7, 7, 7), dtype=np.int64)
    for (i, j, k), c in terms:
        for p in itertools.permutations(range(3)):
            idx = (i, j, k)
            q = tuple(idx[t] for t in p)
            sgn = perm_sign(p)
            phi[q] = c * sgn
    pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    L = [E(7, i, j) - E(7, j, i) for (i, j) in pairs]
    # (X.phi)_{ijk} = -(X_li phi_ljk + X_lj phi_ilk + X_lk phi_ijl)
    rows = []
    for X in L:
        Xp = -(np.einsum('li,ljk->ijk', X, phi) + np.einsum('lj,ilk->ijk', X, phi) + np.einsum('lk,ijl->ijk', X, phi))
        rows.append(Xp.reshape(-1))
    Amat = sp.Matrix(np.array(rows).T.tolist())
    ns = Amat.nullspace()
    assert len(ns) == 14, "stabilizer dimension %d" % len(ns)
    basis = []
    for v in ns:
        den = sp.ilcm(*[sp.Rational(x).q for x in v])
        vi = [int(x * den) for x in v]
        g = int(np.gcd.reduce(np.abs(vi)))
        vi = [x // g for x in vi]
        basis.append(sum(vi[k] * L[k] for k in range(21)))
    C = sp.Matrix([[int(m[i, j]) for (i, j) in pairs] for m in basis]).T      # 21 x 14
    # left inverse on a set of 14 independent rows
    rowsel = []
    for r in range(21):
        if sp.Matrix([C.row(x) for x in rowsel + [r]]).rank() == len(rowsel) + 1:
            rowsel.append(r)
        if len(rowsel) == 14:
            break
    Csel = sp.Matrix([C.row(x) for x in rowsel])
    Cinv = Csel.inv()
    d = 14
    fuq = [[[None] * d for _ in range(d)] for _ in range(d)]
    for a in range(d):
        for b in range(d):
            c = basis[a] @ basis[b] - basis[b] @ basis[a]
            vec = sp.Matrix([int(c[pairs[x][0], pairs[x][1]]) for x in rowsel])
            co = Cinv * vec
            assert all(sum(co[k] * int(basis[k][i, j]) for k in range(d)) == int(c[i, j]) for i in range(7) for j in range(7))
            for k in range(d):
                fuq[a][b][k] = sp.Rational(co[k])
    den = sp.ilcm(*[x.q for a in fuq for b in a for x in b])
    # rescale the basis X -> den X so that structure constants are integers: [dX, dY] = d^2[X,Y] = d (d [X,Y])
    fu = np.array([[[int(x * den) for x in b] for b in a] for a in fuq], dtype=np.int64)
    mats = [den * m for m in basis]
    kap = np.array([[int(np.trace(mats[a] @ mats[b])) for b in range(d)] for a in range(d)], dtype=np.int64)
    Lg = LA("G2(oct)", fu, kap, "tr_7(XY), 7-dim representation of G2 inside so(7)", cartan=None)
    Lg.mats = mats
    return Lg


# ----------------------------------------------------------------------------- catalogue
def build(key):
    if key == 'A1':
        return kac('sl2(Kac A1)', cartan_A(1), orient(cartan_A(1), []))
    if key == 'A2':
        return kac('sl3(Kac A2)', cartan_A(2), orient(cartan_A(2), [(0, 1)]))
    if key == 'A3':
        return kac('sl4(Kac A3)', cartan_A(3), orient(cartan_A(3), [(0, 1), (2, 1)]))
    if key == 'A4':
        return kac('sl5(Kac A4)', cartan_A(4), orient(cartan_A(4), [(0, 1), (1, 2), (2, 3)]))
    if key == 'A5':
        return kac('sl6(Kac A5)', cartan_A(5), orient(cartan_A(5), [(0, 1), (1, 2), (4, 3), (3, 2)]))
    if key == 'D4':
        return kac('so8(Kac D4)', cartan_D(4), orient(cartan_D(4), [(0, 1), (2, 1), (3, 1)]))
    if key == 'E6':
        return kac('E6(Kac)', cartan_E6(), orient(cartan_E6(), [(0, 2), (2, 3), (5, 4), (4, 3), (1, 3)]))
    if key == 'C2fold':
        return fold('sp4(fold A3)', build('A3'), [2, 1, 0])
    if key == 'C3fold':
        return fold('sp6(fold A5)', build('A5'), [4, 3, 2, 1, 0])
    if key == 'B3fold':
        return fold('so7(fold D4)', build('D4'), [0, 1, 3, 2])
    if key == 'G2fold':
        return fold('G2(fold D4)', build('D4'), [2, 1, 3, 0])
    if key == 'F4fold':
        return fold('F4(fold E6)', build('E6'), [5, 1, 4, 3, 2, 0])
    if key.startswith('so'):
        return so_matrix(int(key[2:]))
    if key.startswith('sp'):
        return sp_matrix(int(key[2:]) // 2)
    if key == 'G2oct':
        return g2_octonion()
    if key.startswith('run1-sl'):
        # run-1 basis and trace form (explore/anomaly-structure/code/lie.py, imported read-only)
        import sys as _s
        _s.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
        import lie as _lie
        N = int(key[7:])
        g = _lie.sl(N)
        return LA('sl%d(run-1 basis)' % N, g.fu, g.kappa, 'tr_%d(XY), run-1 lie.sl(%d)' % (N, N),
                  cartan=list(range(N * (N - 1), g.dim)))
    raise KeyError(key)


def cartan_matrix_of(L):
    """Cartan matrix of the root system of L from the weights (exact), simple roots w.r.t. a generic functional."""
    r = len(L.cartan)
    kh = sp.Matrix([[int(L.kappa[a, b]) for b in L.cartan] for a in L.cartan])
    khinv = kh.inv()
    roots = sorted({w for w in L.wt if any(w)})
    gen = [sp.Integer(1000003) ** i for i in range(r)][::-1]
    # weight w is the functional h_k -> w_k; as element of h* its Gram is w^T kh^{-1} w
    val = lambda w: sum(sp.Integer(x) * g for x, g in zip(w, gen))
    pos = [w for w in roots if val(w) > 0]
    possum = {tuple(a + b for a, b in zip(x, y)) for x in pos for y in pos}
    simple = [w for w in pos if w not in possum]
    ip = lambda x, y: (sp.Matrix([x]) * khinv * sp.Matrix(y))[0]
    C = [[2 * ip(a, b) / ip(b, b) for b in simple] for a in simple]
    return sp.Matrix(C), len(roots), len(pos)
