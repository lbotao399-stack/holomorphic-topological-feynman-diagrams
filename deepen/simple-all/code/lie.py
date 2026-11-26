"""Sparse exact simple Lie algebras.

SLA data (all exact integers):
  fu : COO arrays (A, B, C, V) over all ordered pairs:  [e_a, e_b] = sum V e_c   (antisymmetry checked)
  kap: sparse symmetric integer invariant form kappa_ab (checked: f_low(a,b,c) = kappa([e_a,e_b],e_c)
       antisymmetric in (b,c), i.e. ad-invariance), G = s kappa^{-1} integer (checked G kappa = s 1)
  wt : weights (eigenvalues of ad h_i, h_i in the Cartan basis elements), theta = highest root vector
Jacobi: T1(a,b,c,f) = sum_e fu(a,b,e) fu(e,c,f);  T1(a,b,c,.) + T1(b,c,a,.) + T1(c,a,b,.) = 0 on all entries.

Constructions
  kac(A, M)   simply-laced, Kac 'Infinite dimensional Lie algebras' Sec. 7.8 (as in simple-g/code/algebras.kac):
              [h,E_a] = (h|a)E_a, [E_a,E_-a] = -a, [E_a,E_b] = eps(a,b)E_{a+b}, eps(a,b) = (-1)^{a^T M b},
              (h_i|h_j) = A_ij, (E_a|E_-a) = -1.  Jacobi and invariance are checked, not assumed.
  fold(L, sigma)  fixed points of the diagram automorphism E_a -> E_{sigma a} (automorphism checked), basis
              orbit sums, form = restriction.
  from_dense(L)   conversion of a simple-g/code/algebras.LA object (dense arrays).
"""
import itertools, math, time
import numpy as np
import scipy.sparse as sps
import sympy as sp
from fractions import Fraction as Fr


class SLA:
    def __init__(self, name, d, A, B, C, V, kap_entries, form, cartan):
        self.name, self.dim, self.form = name, d, form
        A, B, C, V = [np.asarray(x, dtype=np.int64) for x in (A, B, C, V)]
        key = (A * d + B) * d + C
        o = np.argsort(key, kind='stable')
        key, A, B, C, V = key[o], A[o], B[o], C[o], V[o]
        u, st = np.unique(key, return_index=True)
        Vs = np.add.reduceat(V, st) if len(V) else V
        keep = Vs != 0
        self.A, self.B, self.C, self.V = A[st][keep], B[st][keep], C[st][keep], Vs[keep]
        self.cartan = list(cartan)
        ki, kj, kv = [np.asarray(x, dtype=np.int64) for x in zip(*kap_entries)]
        self.kap = sps.coo_matrix((kv, (ki, kj)), shape=(d, d)).tocsr()
        self.kap.sum_duplicates()
        self._mats()
        self._check()
        self._weights()

    # ------------------------------------------------------------------ matrices
    def _mats(self):
        d = self.dim
        self.FU = sps.csr_matrix((self.V, (self.A * d + self.B, self.C)), shape=(d * d, d))      # [(a,b), c]
        self.FUr = sps.csr_matrix((self.V, (self.A, self.B * d + self.C)), shape=(d, d * d))     # [a, (b,c)]
        FL = (self.FU @ self.kap).tocoo()                                                       # f_low[(a,b), c]
        self.FL = FL.tocsr()

    def ad(self, a):
        """matrix of ad e_a: column b -> [e_a, e_b]."""
        d = self.dim
        m = self.A == a
        return sps.csr_matrix((self.V[m], (self.C[m], self.B[m])), shape=(d, d))

    def _check(self):
        d = self.dim
        # antisymmetry
        k1 = (self.A * d + self.B) * d + self.C
        k2 = (self.B * d + self.A) * d + self.C
        o1, o2 = np.argsort(k1), np.argsort(k2)
        assert np.array_equal(k1[o1], k2[o2]) and np.array_equal(self.V[o1], -self.V[o2]), "antisymmetry"
        # Jacobi
        T1 = (self.FU @ self.FUr).tocoo()                 # rows (a,b), cols (c,f): coefficient of e_f in [[a,b],c]
        a, b = T1.row // d, T1.row % d
        c, f = T1.col // d, T1.col % d
        v = T1.data.astype(np.int64)
        keys, vals = [], []
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            keys.append(((x * d + y) * d + z) * d + f); vals.append(v)
        keys, vals = np.concatenate(keys), np.concatenate(vals)
        o = np.argsort(keys)
        keys, vals = keys[o], vals[o]
        u, st = np.unique(keys, return_index=True)
        s = np.add.reduceat(vals, st) if len(vals) else vals
        assert not np.any(s), "Jacobi fails"
        self.jacobi_terms = len(T1.data)
        # kappa symmetric, invariant
        K = self.kap
        assert (K - K.T).count_nonzero() == 0, "kappa symmetric"
        FL = self.FL.tocoo()
        a, b, c = FL.row // d, FL.row % d, FL.col
        k1 = (a * d + b) * d + c
        k2 = (a * d + c) * d + b
        o1, o2 = np.argsort(k1), np.argsort(k2)
        assert np.array_equal(k1[o1], k2[o2]) and np.array_equal(FL.data[o1], -FL.data[o2]), "kappa not invariant"
        self.f_low = (a, b, c, FL.data.astype(np.int64))

    def _weights(self):
        d = self.dim
        wt = np.zeros((d, len(self.cartan)), dtype=np.int64)
        for k, h in enumerate(self.cartan):
            M = self.ad(h).tocoo()
            assert np.all(M.row == M.col), "Cartan element %d not diagonal" % h
            wt[M.row, k] = M.data
        self.wt = wt
        self._inverse_form()
        # highest root vector: maximize a generic functional
        gen = np.array([1000003 ** i for i in range(len(self.cartan))][::-1], dtype=object)
        val = [int(np.dot(np.array(w, dtype=object), gen)) for w in wt]
        self.theta = int(np.argmax(np.array(val, dtype=object)))

    def _inverse_form(self):
        """G = s kappa^{-1}: kappa pairs weight w with -w (invariance); invert block by block exactly."""
        d = self.dim
        groups = {}
        for a in range(d):
            groups.setdefault(tuple(self.wt[a]), []).append(a)
        Kd = self.kap.tolil()
        blocks = []
        den = 1
        done = set()
        for w, idx in groups.items():
            if w in done:
                continue
            mw = tuple(-x for x in w)
            jdx = groups[mw]
            done.add(w); done.add(mw)
            M = sp.Matrix([[int(Kd[i, j]) for j in jdx] for i in idx])        # rows weight w, cols weight -w
            assert M.shape[0] == M.shape[1] and M.det() != 0, "kappa degenerate"
            Mi = M.inv()                                                      # Mi[j,i]: (kappa^{-1})_{jdx, idx}
            blocks.append((idx, jdx, Mi))
            for x in Mi:
                den = sp.ilcm(den, sp.Rational(x).q)
        self.s = int(den)
        gi, gj, gv = [], [], []
        for idx, jdx, Mi in blocks:
            for r, j in enumerate(jdx):
                for c, i in enumerate(idx):
                    x = sp.Rational(Mi[r, c]) * self.s
                    if x != 0:
                        assert x.q == 1
                        gi.append(j); gj.append(i); gv.append(int(x))
                        if idx is not jdx and idx != jdx:
                            gi.append(i); gj.append(j); gv.append(int(x))
        self.G = sps.csr_matrix((np.array(gv, dtype=np.int64), (gi, gj)), shape=(d, d))
        P = (self.G @ self.kap).tocoo()
        assert np.all(P.row == P.col) and np.all(P.data == self.s) and len(P.data) == d, "G kappa != s 1"
        assert (self.G - self.G.T).count_nonzero() == 0

    def killing_ratio(self):
        """r with tr(ad a ad b) = r kappa_ab, checked on all entries."""
        d = self.dim
        # K[a,b] = sum_{c,e} fu[a,c,e] fu[b,e,c]
        M1 = sps.csr_matrix((self.V, (self.A, self.B * d + self.C)), shape=(d, d * d))    # [a,(c,e)]
        M2 = sps.csr_matrix((self.V, (self.A, self.C * d + self.B)), shape=(d, d * d))    # [b,(e,c)] -> col (c,e)
        K = (M1 @ M2.T).tocsr()
        K.sum_duplicates()
        Kc, kc = K.tocoo(), self.kap.tocoo()
        i, j = kc.row[0], kc.col[0]
        r = Fr(int(K[i, j]), int(kc.data[0]))
        D = (K * r.denominator - self.kap * r.numerator)
        D.eliminate_zeros()
        return r, D.count_nonzero() == 0


# ------------------------------------------------------------------ Cartan matrices
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


def cartan_E(n):
    """Bourbaki E_n (n = 6,7,8): 1-3-4-5-..-n chain, 2 attached to 4 (0-based: 0-2-3-..-(n-1), 1-3)."""
    A = 2 * np.eye(n, dtype=int)
    edges = [(0, 2)] + [(k, k + 1) for k in range(2, n - 1)] + [(1, 3)]
    for i, j in edges:
        A[i, j] = A[j, i] = -1
    return A


def orient_all(A, center=0):
    """M_ii = 1, M_ij = 1 for an edge i -> j; every edge of the (tree) diagram is oriented towards `center`
    (so every diagram automorphism fixing `center` preserves eps)."""
    n = len(A)
    dist = {center: 0}
    fr = [center]
    while fr:
        new = []
        for i in fr:
            for j in range(n):
                if A[i, j] and j != i and j not in dist:
                    dist[j] = dist[i] + 1; new.append(j)
        fr = new
    M = np.eye(n, dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and A[i, j] and dist[j] < dist[i]:
                M[i, j] = 1
    assert sum(M[i, j] + M[j, i] for i in range(n) for j in range(i + 1, n) if A[i, j]) == \
        sum(1 for i in range(n) for j in range(i + 1, n) if A[i, j])
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


def kac(name, A, M=None, center=0):
    A = np.array(A); n = len(A)
    if M is None:
        M = orient_all(A, center)
    pos = positive_roots_simply_laced(A)
    roots = pos + [tuple(-x for x in r) for r in pos]
    idx = {r: n + k for k, r in enumerate(roots)}
    d = n + len(roots)
    AA, BB, CC, VV = [], [], [], []
    kap = [(i, j, int(A[i, j])) for i in range(n) for j in range(n) if A[i, j]]
    for r in roots:
        a = idx[r]
        rv = np.array(r)
        Ar = A @ rv
        for i in range(n):
            if Ar[i]:
                AA += [i, a]; BB += [a, i]; CC += [a, a]; VV += [int(Ar[i]), -int(Ar[i])]
        neg = tuple(-x for x in r)
        kap.append((a, idx[neg], -1))
        for i in range(n):
            if rv[i]:
                AA.append(a); BB.append(idx[neg]); CC.append(i); VV.append(-int(rv[i]))
        for s in roots:
            t = tuple(x + y for x, y in zip(r, s))
            if t in idx:
                AA.append(a); BB.append(idx[s]); CC.append(idx[t]); VV.append((-1) ** int(rv @ M @ np.array(s)))
    L = SLA(name, d, AA, BB, CC, VV, kap, "Kac form (h_i|h_j) = A_ij, (E_a|E_-a) = -1", list(range(n)))
    L.roots, L.ridx, L.rank_, L.cartan_matrix = roots, idx, n, A
    hts = [sum(r) for r in pos]
    L.theta = idx[pos[int(np.argmax(hts))]]
    return L


def fold(name, L, sigma):
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
    d = L.dim
    # automorphism check on structure constants and form
    k1 = (L.A * d + L.B) * d + L.C
    k2 = (P[L.A] * d + P[L.B]) * d + P[L.C]
    o1, o2 = np.argsort(k1), np.argsort(k2)
    assert np.array_equal(k1[o1], k2[o2]) and np.array_equal(L.V[o1], L.V[o2]), "not an automorphism"
    kc = L.kap.tocoo()
    assert all(L.kap[P[i], P[j]] == v for i, j, v in zip(kc.row, kc.col, kc.data))
    seen, orbits = set(), []
    for a in range(d):
        if a in seen:
            continue
        o, b = [], a
        while b not in o:
            o.append(b); b = perm[b]
        seen |= set(o)
        orbits.append(sorted(o))
    orb_of = np.zeros(d, dtype=np.int64)
    for k, o in enumerate(orbits):
        orb_of[o] = k
    m = len(orbits)
    # [B_I, B_J] = sum_{a in I, b in J} [e_a, e_b]; coefficient of B_K = coefficient of rep(K)
    acc = {}
    for a, b, c, v in zip(L.A, L.B, L.C, L.V):
        key = (orb_of[a], orb_of[b], int(c))
        acc[key] = acc.get(key, 0) + int(v)
    out = {}
    for (I, J, c), v in acc.items():
        out.setdefault((I, J, orb_of[c]), {})[c] = v
    AA, BB, CC, VV = [], [], [], []
    for (I, J, K), dct in out.items():
        vals = set(dct.get(c, 0) for c in orbits[K])
        assert len(vals) == 1, "bracket not in the fixed subalgebra"
        v = vals.pop()
        if v:
            AA.append(I); BB.append(J); CC.append(K); VV.append(v)
    kap = {}
    for i, j, v in zip(kc.row, kc.col, kc.data):
        key = (orb_of[i], orb_of[j])
        kap[key] = kap.get(key, 0) + int(v)
    kap = [(i, j, v) for (i, j), v in kap.items() if v]
    cart = [k for k, o in enumerate(orbits) if o[0] < n]
    F = SLA(name, m, AA, BB, CC, VV, kap, "restriction of the %s Kac form to the fixed subalgebra" % L.name, cart)
    F.orbits = orbits
    return F


def from_dense(L):
    fu = L.fu
    A, B, C = np.nonzero(fu)
    V = fu[A, B, C]
    ki, kj = np.nonzero(L.kappa)
    kap = [(int(i), int(j), int(L.kappa[i, j])) for i, j in zip(ki, kj)]
    S = SLA(L.name, L.dim, A, B, C, V, kap, L.form, L.cartan)
    return S


SG = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code'


def build(key):
    """keys: A<n>, D<n>, E6, E7, E8 (Kac); B<n> = fold D<n+1>, C<n> = fold A<2n-1>, F4 = fold E6, G2 = fold D4;
    'sg:<key>' = simple-g/code/algebras.build(key) converted (dense constructions)."""
    if key.startswith('sg:'):
        import sys
        if SG not in sys.path:
            sys.path.insert(0, SG)
        import algebras as al
        return from_dense(al.build(key[3:]))
    t, n = key[0], key[1:]
    if key in ('E6', 'E7', 'E8'):
        return kac(key + '(Kac)', cartan_E(int(n)), center=3)
    if key == 'F4':
        return fold('F4(fold E6)', build('E6'), [5, 1, 4, 3, 2, 0])
    if key == 'G2':
        return fold('G2(fold D4)', build('D4'), [2, 1, 3, 0])
    n = int(n)
    if t == 'A':
        return kac('sl%d(Kac A%d)' % (n + 1, n), cartan_A(n), center=(n - 1) // 2)
    if t == 'D':
        return kac('so%d(Kac D%d)' % (2 * n, n), cartan_D(n), center=n - 3 if n >= 4 else 0)
    if t == 'B':
        # D_{n+1}: nodes 0..n; the two end nodes n-1, n are exchanged
        sig = list(range(n + 1)); sig[n - 1], sig[n] = n, n - 1
        return fold('so%d(fold D%d)' % (2 * n + 1, n + 1), build('D%d' % (n + 1)), sig)
    if t == 'C':
        # A_{2n-1}: i <-> 2n-2-i
        m = 2 * n - 1
        return fold('sp%d(fold A%d)' % (2 * n, m), build('A%d' % m), [m - 1 - i for i in range(m)])
    raise KeyError(key)


# ------------------------------------------------------------------ matrix algebras with kappa = tr_V (sparse)
def _mat_alg(name, N, basis, coords, form):
    """basis: list of dicts {(i,j): v}; coords(m) -> dict k -> coefficient (exact, verified by reconstruction)."""
    d = len(basis)
    rows = {}
    for k, m in enumerate(basis):
        for (i, j), v in m.items():
            rows.setdefault(i, []).append((k, j, v))          # entries in row i
    AA, BB, CC, VV = [], [], [], []
    for a in range(d):
        ma = basis[a]
        for b in range(a + 1, d):
            mb = basis[b]
            c = {}
            for (i, j), v in ma.items():
                for (j2, l), w in mb.items():
                    if j2 == j:
                        c[(i, l)] = c.get((i, l), 0) + v * w
            for (i, j), v in mb.items():
                for (j2, l), w in ma.items():
                    if j2 == j:
                        c[(i, l)] = c.get((i, l), 0) - v * w
            c = {k: v for k, v in c.items() if v}
            if not c:
                continue
            co = coords(c)
            rec = {}
            for k, x in co.items():
                for e, v in basis[k].items():
                    rec[e] = rec.get(e, 0) + x * v
            rec = {k: v for k, v in rec.items() if v}
            assert rec == c, "bracket not in span"
            for k, x in co.items():
                AA += [a, b]; BB += [b, a]; CC += [k, k]; VV += [x, -x]
    kap = []
    for a in range(d):
        for b in range(d):
            t = 0
            for (i, j), v in basis[a].items():
                w = basis[b].get((j, i))
                if w:
                    t += v * w
            if t:
                kap.append((a, b, t))
    cart = [k for k, m in enumerate(basis) if all(i == j for (i, j) in m)]
    L = SLA(name, d, AA, BB, CC, VV, kap, form, cart)
    L.N = N
    L.basis = basis
    return L


def so_tr(N):
    pr = lambda k: N - 1 - k
    basis, reps = [], []
    for i in range(N):
        for j in range(N):
            if i + j < N - 1:
                m = {(i, j): 1}
                m[(pr(j), pr(i))] = m.get((pr(j), pr(i)), 0) - 1
                basis.append({k: v for k, v in m.items() if v}); reps.append((i, j))
    idx = {r: k for k, r in enumerate(reps)}
    def coords(c):
        return {idx[e]: v for e, v in c.items() if e in idx}
    return _mat_alg('so%d(tr)' % N, N, basis, coords, 'tr_%d(XY), so(%d) w.r.t. antidiagonal J' % (N, N))


def sp_tr(n):
    N = 2 * n
    pr = lambda k: N - 1 - k
    s = lambda k: 1 if k < n else -1
    basis, reps = [], []
    for i in range(N):
        for j in range(N):
            if i + j < N - 1:
                m = {(i, j): 1}
                m[(pr(j), pr(i))] = m.get((pr(j), pr(i)), 0) + s(pr(j)) * s(i)
                basis.append({k: v for k, v in m.items() if v}); reps.append((i, j))
            elif i + j == N - 1:
                basis.append({(i, j): 1}); reps.append((i, j))
    idx = {r: k for k, r in enumerate(reps)}
    def coords(c):
        return {idx[e]: v for e, v in c.items() if e in idx}
    return _mat_alg('sp%d(tr)' % N, N, basis, coords, 'tr_%d(XY), sp(%d) w.r.t. antidiagonal symplectic J' % (N, N))


def sl_tr(N):
    basis = []
    idx = {}
    for i in range(N):
        for j in range(N):
            if i != j:
                idx[(i, j)] = len(basis); basis.append({(i, j): 1})
    h0 = len(basis)
    for i in range(N - 1):
        basis.append({(i, i): 1, (i + 1, i + 1): -1})
    def coords(c):
        out = {idx[e]: v for e, v in c.items() if e[0] != e[1]}
        cum = 0
        for i in range(N - 1):
            cum += c.get((i, i), 0)
            if cum:
                out[h0 + i] = cum
        return out
    return _mat_alg('sl%d(tr)' % N, N, basis, coords, 'tr_%d(XY), sl(%d)' % (N, N))
