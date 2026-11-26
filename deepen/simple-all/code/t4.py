"""Exact sparse invariant 4-tensors, two-loop colour weights, cochains and probes.

A T4 holds integer COO entries t(a0,a1,a2,a3) (all indices lower); its value is t / s^sp.
  perm(T, sig):   T'(x0,x1,x2,x3) = T(x_sig0, x_sig1, x_sig2, x_sig3)
  compose(A, B):  (A o B)(a,b,c,e) = sum A(a,b,x,y) kappa^{xx'} kappa^{yy'} B(x',y',c,e)        (sp_A + sp_B + 2)
  pair(A, B):     <A,B> = sum A(a,b,c,e) kappa^{aa'} kappa^{bb'} kappa^{cc'} kappa^{ee'} B(a',b',c',e')  (Fraction)
Integer sparse products: |A| |B| is evaluated in float64 first; the int64 product is used only if that bound
is < 2^61 (then every partial sum is an exact int64), otherwise the product is computed modulo primes
p < 2^21 (int64, no overflow: |entry| < 2^20, partial sums < 2^62) and reconstructed by CRT.

Building blocks (L = lie.SLA):
  f(a,b,c)       = kappa([e_a,e_b],e_c)
  Om(a,b,c,e)    = sum_{x,x'} f(x,a,c) kappa^{xx'} f(x',b,e)          (rung: rails a->c, b->e)
  ff(a,b,c,e)    = kappa([e_a,e_b],[e_c,e_e])
  kk_{ab|ce}     = kappa_ab kappa_ce
  W133 = - Om o Om o Om,   W223 = (Om2 o sig) o Om,  Om2 = Om o Om, sig = (2,1,0,3)
  (index forms of simple-g/code/core.py; equality is checked against that code, logs/check_W.log).
"""
import itertools, math
import numpy as np
import scipy.sparse as sps
from fractions import Fraction as Fr

LIM = 2.0 ** 61
PR = [2097143, 2097133, 2097131, 2097097, 2097091, 2097083, 2097047, 2097041, 2097031, 2097023,
      2097013, 2096993, 2096971, 2096959, 2096957, 2096947]


def perm_sign(seq):
    seq = list(seq); s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def _crt_lists(res, primes):
    """res: list of int64 arrays (same positions), symmetric residues -> exact Python ints (object array)."""
    M = 1
    x = np.zeros(len(res[0]), dtype=object)
    for r, p in zip(res, primes):
        r = r.astype(object) % p
        # x = x + M * t, t = (r - x) * M^{-1} mod p
        inv = pow(M % p, p - 2, p)
        t = ((r - x) % p) * inv % p
        x = x + M * t
        M *= p
    half = M // 2
    x = np.where(x > half, x - M, x)
    return x, M


def smm(A, B):
    """exact product of two integer csr matrices; returns (rows, cols, values as int64 or object)."""
    Af = A.astype(np.float64); Af.data = np.abs(Af.data)
    Bf = B.astype(np.float64); Bf.data = np.abs(Bf.data)
    bound = (Af @ Bf)
    mx = bound.data.max() if bound.nnz else 0.0
    if mx < LIM:
        C = (A @ B).tocoo()
        C.sum_duplicates()
        m = C.data != 0
        return C.row[m], C.col[m], C.data[m].astype(np.int64)
    need = 2.0 * mx * 4
    primes, prod = [], 1.0
    for p in PR:
        primes.append(p); prod *= p
        if prod > need:
            break
    assert prod > need
    bc = bound.tocoo()
    bc.sum_duplicates()
    rows, cols = bc.row, bc.col
    order = np.lexsort((cols, rows))
    rows, cols = rows[order], cols[order]
    res = []
    for p in primes:
        Ap = A.copy(); Ap.data = ((Ap.data % p + p // 2) % p - p // 2).astype(np.int64)
        Bp = B.copy(); Bp.data = ((Bp.data % p + p // 2) % p - p // 2).astype(np.int64)
        C = (Ap @ Bp).tocsr()
        C.sum_duplicates()
        vals = np.asarray(C[rows, cols]).reshape(-1).astype(np.int64)
        res.append(vals % p)
    x, M = _crt_lists(res, primes)
    m = x != 0
    return rows[m], cols[m], x[m]


class T4:
    def __init__(self, d, idx, val, sp=0):
        self.d = d
        self.idx = np.asarray(idx, dtype=np.int64).reshape(-1, 4)
        self.val = np.asarray(val)
        self.sp = sp
        self._canon()

    def _canon(self):
        d = self.d
        k = self.keys()
        o = np.argsort(k, kind='stable')
        k = k[o]
        u, st = np.unique(k, return_index=True)
        v = self.val[o]
        if v.dtype == object:
            vs = np.array([sum(v[a:b]) for a, b in zip(st, list(st[1:]) + [len(v)])], dtype=object) if len(v) else v
        else:
            vs = np.add.reduceat(v, st) if len(v) else v
        idx = self.idx[o][st]
        m = vs != 0
        self.idx, self.val = idx[m], vs[m]
        if self.val.dtype == object and len(self.val) and max(abs(int(x)) for x in self.val) < 2 ** 62:
            self.val = self.val.astype(np.int64)

    def keys(self):
        d = self.d
        i = self.idx
        return ((i[:, 0] * d + i[:, 1]) * d + i[:, 2]) * d + i[:, 3]

    def mat(self, rows=(0, 1), cols=(2, 3)):
        d = self.d
        r = self.idx[:, rows[0]] * d + self.idx[:, rows[1]]
        c = self.idx[:, cols[0]] * d + self.idx[:, cols[1]]
        assert self.val.dtype != object, "object values: use pair_obj"
        return sps.csr_matrix((self.val, (r, c)), shape=(d * d, d * d))

    def nnz(self):
        return len(self.val)


def from_mat(d, rows, cols, vals, legs, sp):
    """entries M[(x_l0,x_l1),(x_l2,x_l3)] with legs = positions (l0,l1,l2,l3) of the tensor slots."""
    idx = np.zeros((len(vals), 4), dtype=np.int64)
    idx[:, legs[0]] = rows // d
    idx[:, legs[1]] = rows % d
    idx[:, legs[2]] = cols // d
    idx[:, legs[3]] = cols % d
    return T4(d, idx, vals, sp)


def perm(T, sig):
    """T'(x0..x3) = T(x_sig0..x_sig3): an entry of T at y has y_k = x_{sig k}."""
    idx = np.zeros_like(T.idx)
    for k in range(4):
        idx[:, sig[k]] = T.idx[:, k]
    return T4(T.d, idx, T.val.copy(), T.sp)


def lin(terms, d=None, sp=None):
    """sum c_i T_i (integer c_i, common sp)."""
    if not terms:
        return T4(d, np.zeros((0, 4), dtype=np.int64), np.zeros(0, dtype=np.int64), sp)
    sps_ = {T.sp for c, T in terms}
    assert len(sps_) == 1
    d = terms[0][1].d
    idx = np.concatenate([T.idx for c, T in terms])
    obj = any(T.val.dtype == object for c, T in terms)
    vals = [(T.val.astype(object) if obj else T.val) * c for c, T in terms]
    return T4(d, idx, np.concatenate(vals), sps_.pop())


def scale(T, c):
    return T4(T.d, T.idx.copy(), T.val * c, T.sp)


class Engine:
    def __init__(self, L):
        self.L = L
        d = self.d = L.dim
        self.s = L.s
        G = L.G.astype(np.int64)
        self.G = G
        self.GG = sps.kron(G, G, format='csr').astype(np.int64)
        a, b, c, v = L.f_low
        self.f = (a, b, c, v)
        # Om: Fm[(a,c), x] = f(x,a,c);  Om_mat[(a,c),(b,e)] = Fm G Fm^T
        Fm = sps.csr_matrix((v, (b * d + c, a)), shape=(d * d, d))
        r1, c1, v1 = smm(Fm, G)
        FG = sps.csr_matrix((v1, (r1, c1)), shape=(d * d, d))
        r, cc, vv = smm(FG, Fm.T.tocsr())
        self.Om = from_mat(d, r, cc, vv, (0, 2, 1, 3), 1)
        # ff(a,b,c,e) = sum_x fu(a,b,x) f(x,c,e)
        FU = L.FU.astype(np.int64)
        FLr = sps.csr_matrix((v, (a, b * d + c)), shape=(d, d * d))
        r, cc, vv = smm(FU, FLr)
        self.ff = from_mat(d, r, cc, vv, (0, 1, 2, 3), 0)
        kc = L.kap.tocoo()
        ki, kj, kv = kc.row, kc.col, kc.data.astype(np.int64)
        n = len(kv)
        I, J = np.repeat(np.arange(n), n), np.tile(np.arange(n), n)
        idx = np.stack([ki[I], kj[I], ki[J], kj[J]], 1)
        self.kk = T4(d, idx, kv[I] * kv[J], 0)          # kappa_01 kappa_23

    def compose(self, A, B):
        r, c, v = smm(A.mat(), self.GG)
        AG = sps.csr_matrix((v, (r, c)), shape=(self.d ** 2, self.d ** 2))
        r, c, v = smm(AG, B.mat())
        return from_mat(self.d, r, c, v, (0, 1, 2, 3), A.sp + B.sp + 2)

    def raise_all(self, T):
        """T^{up}(a,b,c,e) * s^4 (integer), returned as T4 with sp = T.sp + 4 (value of the raised tensor)."""
        r, c, v = smm(self.GG, T.mat())
        X = sps.csr_matrix((v, (r, c)), shape=(self.d ** 2, self.d ** 2))
        r, c, v = smm(X, self.GG)
        return from_mat(self.d, r, c, v, (0, 1, 2, 3), T.sp + 4)

    def pair(self, A, Bup):
        """<A, B> with Bup = raise_all(B) (precomputed).  Exact Fraction."""
        ka, kb = A.keys(), Bup.keys()
        common, ia, ib = np.intersect1d(ka, kb, assume_unique=True, return_indices=True)
        tot = sum(int(x) * int(y) for x, y in zip(A.val[ia], Bup.val[ib]))
        return Fr(tot, self.s ** (A.sp + Bup.sp))

    # ------------------------------------------------------------------ colour weights
    def weights(self):
        Om = self.Om
        Om2 = self.compose(Om, Om)
        Om3 = self.compose(Om2, Om)
        W133 = scale(Om3, -1)
        W223 = self.compose(perm(Om2, (2, 1, 0, 3)), Om)
        self.Om2, self.Om3 = Om2, Om3
        return W133, W223


# ---------------------------------------------------------------------------- momentum polynomials
QB = ['uu', 'vv', 'ww', 'uv', 'uw', 'vw']


def br(a, b):
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
    out = np.zeros(6, dtype=np.int64)
    pos = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3, (0, 2): 4, (1, 2): 5}
    for i in range(3):
        for j in range(3):
            if l1[i] and l2[j]:
                out[pos[(min(i, j), max(i, j))]] += l1[i] * l2[j]
    return out


def alt_bb(W):
    """F^mu = sum_s sgn(s) quad([p_s0 p_s1],[p_s2 p_s3])_mu W o s   (mu over QB) -> list of 6 T4."""
    terms = [[] for _ in range(6)]
    for sig in itertools.permutations(range(4)):
        c = quad(br(sig[0], sig[1]), br(sig[2], sig[3]))
        if not c.any():
            continue
        P = perm(W, sig)
        sg = perm_sign(sig)
        for mu in range(6):
            if c[mu]:
                terms[mu].append((sg * int(c[mu]), P))
    return [lin(t, W.d, W.sp) for t in terms]


def beta_t(ff):
    """beta^mu = sum_{i<j} (-1)^{i+j} quad([p_k p_l],[p_k p_l])_mu ff o (i,j,k,l)."""
    terms = [[] for _ in range(6)]
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = quad(br(k, l), br(k, l))
        P = perm(ff, (i, j, k, l))
        for mu in range(6):
            if c[mu]:
                terms[mu].append(((-1) ** (i + j) * int(c[mu]), P))
    return [lin(t, ff.d, ff.sp) for t in terms]


def merge(Ts):
    """align several T4 (same sp) on the union of keys -> (keys, list of value arrays as Python ints)."""
    keys = np.unique(np.concatenate([T.keys() for T in Ts]))
    out = []
    for T in Ts:
        v = np.zeros(len(keys), dtype=object)
        pos = np.searchsorted(keys, T.keys())
        v[pos] = T.val.astype(object)
        out.append(v)
    return keys, out
