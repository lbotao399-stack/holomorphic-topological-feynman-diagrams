"""Part (c): H^3_(a,b) of G-invariant, translation-invariant cochains, 0 <= a, b <= 2, explicit structure constants.

Cochains of arity k: sum_t m_t(p_0..p_{k-2}) T_t(X_0..X_{k-1}) (p_{k-1} = -(p_0+..+p_{k-2}) eliminated), T_t sparse
integer k-tensors.  Bidegree (a,b): degree a in the components p_i^1, b in p_i^2.
  C^2_(a,b) = span alt( q0^m kappa ),  C^3_(a,b) = span alt( q0^m q1^m' T ),  T over a basis of (g^{*(x)3})^g
  (for T^3 = 1 algebras: T = f; for sl_N: f and d = tr(X{Y,Z}) computed from the matrices of lie.sl_tr).
  (c.s)(x_0..) = c(x_s0, x_s1, ..),  alt c = sum_s sgn(s) c.s
  (d c)(x_0..x_k) = sum_{i<j} (-1)^{i+j} c([x_i,x_j], x_0..^i..^j..x_k),  [X e^{pz}, Y e^{qz}] = [X,Y] e^{(p+q)z}
Ranks exact: rank mod p from elimination on the aligned coefficient rows, then the (rows - r) left-kernel vectors
of the pivot submatrix are verified exactly on all columns (Python integers), so rank over Q = r.
H^3_(a,b) = dim ker(d: C^3 -> C^4) - rank(d: C^2 -> C^3).
"""
import sys, time, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import scipy.sparse as sps
import sympy as sp
from fractions import Fraction as Fr
import lie, t4

PV = [[sp.Symbol('p%d_%d' % (i, c)) for c in (1, 2)] for i in range(6)]


class Tk:
    def __init__(self, d, idx, val):
        self.d = d
        self.idx = np.asarray(idx, dtype=np.int64)
        self.k = self.idx.shape[1]
        self.val = np.asarray(val, dtype=np.int64)
        k = self.keys()
        o = np.argsort(k, kind='stable')
        k = k[o]
        u, st = np.unique(k, return_index=True)
        v = np.add.reduceat(self.val[o], st) if len(k) else self.val
        idx = self.idx[o][st] if len(k) else self.idx
        m = v != 0
        self.idx, self.val = idx[m], v[m]

    def keys(self):
        key = np.zeros(len(self.idx), dtype=np.int64)
        for c in range(self.idx.shape[1]):
            key = key * self.d + self.idx[:, c]
        return key


def tperm(T, s):
    """T'(x_0..) = T(x_s0, x_s1, ..)."""
    idx = np.zeros_like(T.idx)
    for k in range(T.k):
        idx[:, s[k]] = T.idx[:, k]
    return Tk(T.d, idx, T.val.copy())


def bracket_first(L, T):
    """T'(a, b, y_1..) = sum_e fu(a,b,e) T(e, y_1..)."""
    d, k = L.dim, T.k
    rest = np.zeros(len(T.idx), dtype=np.int64)
    for c in range(1, k):
        rest = rest * d + T.idx[:, c]
    ncol = d ** (k - 1)
    M = sps.csr_matrix((T.val, (T.idx[:, 0], rest)), shape=(d, ncol))
    R = (L.FU.astype(np.int64) @ M).tocoo()
    idx = np.zeros((len(R.data), k + 1), dtype=np.int64)
    idx[:, 0], idx[:, 1] = R.row // d, R.row % d
    cc = R.col.copy()
    for c in range(k - 1, 0, -1):
        idx[:, c + 1] = cc % d; cc //= d
    return Tk(d, idx, R.data.astype(np.int64))


def poly_dict(expr, k):
    """expand with p_{k-1} eliminated -> dict exponent tuple (over p_0..p_{k-2}, 2 comps) -> int."""
    sub = {PV[k - 1][c]: -sum(PV[i][c] for i in range(k - 1)) for c in range(2)}
    gens = [PV[i][c] for i in range(k - 1) for c in range(2)]
    e = sp.expand(expr.subs(sub, simultaneous=True))
    if e == 0:
        return {}
    P = sp.Poly(e, *gens)
    return {m: int(c) for m, c in P.as_dict().items()}


class Coch:
    """arity k; terms: dict exponent tuple -> list of (int coefficient, Tk)."""
    def __init__(self, k):
        self.k = k
        self.terms = {}

    def add(self, expr, T, coef=1):
        for m, c in poly_dict(expr, self.k).items():
            self.terms.setdefault(m, []).append((coef * c, T))

    def flat(self):
        """dict (m, colour key) -> int."""
        out = {}
        for m, lst in self.terms.items():
            for c, T in lst:
                for key, v in zip(T.keys(), T.val):
                    kk = (m, int(key))
                    out[kk] = out.get(kk, 0) + c * int(v)
        return {k: v for k, v in out.items() if v}


def mono(q, a, b):
    """monomials of bidegree (a,b) in the vectors q (list of [q^1, q^2])."""
    n = len(q)
    out = []
    for ea in itertools.product(range(a + 1), repeat=n):
        if sum(ea) != a:
            continue
        for eb in itertools.product(range(b + 1), repeat=n):
            if sum(eb) != b:
                continue
            m = 1
            for i in range(n):
                m *= q[i][0] ** ea[i] * q[i][1] ** eb[i]
            out.append(m)
    return out


def alt_cochain(k, m_of, T):
    """alt( m(p_0..p_{k-1}) T ) with m_of(list of momentum vectors) -> expr."""
    c = Coch(k)
    for s in itertools.permutations(range(k)):
        sg = t4.perm_sign(s)
        expr = m_of([PV[s[i]] for i in range(k)])
        c.add(expr, tperm(T, s), sg)
    return c


def d_cochain(L, c):
    """(dc)(x_0..x_k) for c of arity k."""
    k = c.k
    out = Coch(k + 1)
    gens = [PV[i][cc] for i in range(k - 1) for cc in range(2)]
    for m, lst in c.terms.items():
        mono_expr = sp.Mul(*[g ** e for g, e in zip(gens, m)])
        for i, j in itertools.combinations(range(k + 1), 2):
            rest = [r for r in range(k + 1) if r not in (i, j)]
            q = [[PV[i][0] + PV[j][0], PV[i][1] + PV[j][1]]] + [PV[r] for r in rest]
            sub = {PV[a][cc]: q[a][cc] for a in range(k - 1) for cc in range(2)}
            expr = mono_expr.subs(sub, simultaneous=True)
            order = [i, j] + rest
            for coef, T in lst:
                TB = bracket_first(L, T)                    # legs (a, b, y_1..y_{k-1}) = (x_i, x_j, x_rest)
                # T'(x_0..x_k) = TB(x_i, x_j, x_rest): perm s with s[pos] = order[pos]
                out.add(expr, tperm(TB, order), (-1) ** (i + j) * coef)
    return out


def rank_exact(cochs, p=2147483629):
    flats = [c.flat() for c in cochs]
    keys = sorted(set().union(*[set(f) for f in flats])) if flats else []
    if not keys:
        return 0
    pos = {k: n for n, k in enumerate(keys)}
    n, ncol = len(flats), len(keys)
    M = np.zeros((n, ncol), dtype=object)
    for r, f in enumerate(flats):
        for k, v in f.items():
            M[r, pos[k]] = v
    # modular elimination (rows few)
    A = np.array([[int(x) % p for x in row] for row in M], dtype=np.int64)
    piv_cols, rows_used = [], []
    A = A.copy()
    r = 0
    for col in range(ncol):
        if r == n:
            break
        nzr = [i for i in range(r, n) if A[i, col] % p]
        if not nzr:
            continue
        i0 = nzr[0]
        A[[r, i0]] = A[[i0, r]]
        inv = pow(int(A[r, col]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for i in range(n):
            if i != r and A[i, col]:
                A[i] = (A[i] - A[i, col] * A[r]) % p
        piv_cols.append(col)
        r += 1
    # exact left kernel of M restricted to pivot columns, verified on all columns
    S = sp.Matrix([[M[i, c] for c in piv_cols] for i in range(n)])
    ker = S.T.nullspace()
    assert len(ker) == n - r
    for v in ker:
        den = sp.ilcm(*[sp.Rational(x).q for x in v])
        w = [int(x * den) for x in v]
        tot = sum(w[i] * M[i] for i in range(n) if w[i])
        assert not any(tot), "kernel vector fails on some column"
    return r


def f_tensor(L):
    a, b, c, v = L.f_low
    return Tk(L.dim, np.stack([a, b, c], 1), v)


def kappa_tensor(L):
    kc = L.kap.tocoo()
    return Tk(L.dim, np.stack([kc.row, kc.col], 1), kc.data)


def d_tensor_sl(L):
    """d(X,Y,Z) = tr(X(YZ+ZY)) from the matrices of lie.sl_tr."""
    B = L.basis
    idx, val = [], []
    d = L.dim
    for a in range(d):
        for b in range(d):
            YZ = {}
            for (i, j), v in B[a].items():
                for (j2, l), w in B[b].items():
                    if j2 == j:
                        YZ[(i, l)] = YZ.get((i, l), 0) + v * w
            for c in range(d):
                t = 0
                for (i, l), v in YZ.items():
                    w = B[c].get((l, i))
                    if w:
                        t += v * w
                if t:
                    idx.append((a, b, c)); val.append(t)
    T = Tk(d, np.array(idx), np.array(val))
    S = None
    terms = [tperm(T, s) for s in itertools.permutations(range(3))]
    idx = np.concatenate([x.idx for x in terms]); val = np.concatenate([x.val for x in terms])
    return Tk(d, idx, val)


def run(L, Ts, label):
    out = ["== %s (dim %d) ; colour basis of (g^{*(x)3})^g used: %s" % (L.name, L.dim, label)]
    kap = kappa_tensor(L)
    table = {}
    for a in range(3):
        for b in range(3):
            t0 = time.time()
            C2 = [alt_cochain(2, (lambda m: (lambda q: m.subs({PV[0][0]: q[0][0], PV[0][1]: q[0][1]}, simultaneous=True)))(m), kap)
                  for m in mono([PV[0]], a, b)]
            C3 = []
            for T in Ts:
                for m in mono([PV[0], PV[1]], a, b):
                    C3.append(alt_cochain(3, (lambda m: (lambda q: m.subs({PV[0][0]: q[0][0], PV[0][1]: q[0][1],
                                                                           PV[1][0]: q[1][0], PV[1][1]: q[1][1]}, simultaneous=True)))(m), T))
            dimC2 = rank_exact(C2)
            dimC3 = rank_exact(C3)
            dC2 = [d_cochain(L, c) for c in C2]
            rk2 = rank_exact(dC2)
            dC3 = [d_cochain(L, c) for c in C3]
            # rank of d on span(C3) = rank of images, given dim C3 = rank of C3 generators:
            # dim ker = dim C3 - rank d3; rank d3 = rank[images] computed on the generators (images of a spanning set).
            rk3 = rank_exact(dC3)
            H3 = (dimC3 - rk3) - rk2
            table[(a, b)] = (dimC2, dimC3, rk2, rk3, H3)
            out.append("   (a,b) = (%d,%d): dim C^2 = %d, dim C^3 = %d, rank d2 = %d, rank d3 = %d, dim Z^3 = %d, dim H^3 = %d (%.1fs)"
                       % (a, b, dimC2, dimC3, rk2, rk3, dimC3 - rk3, H3, time.time() - t0))
    return "\n".join(out), table


def dEf_check(L):
    """E_f = 2 (sum_i p_i^1 p_i^2) f = alt(q0^1 q0^2 f); d E_f != 0 ; d E_f against the closed formula
    dE_f = 2 sum_{i<j} (-1)^{i+j} (p_i^1 p_j^2 + p_i^2 p_j^1) kappa([X_i,X_j],[X_k,X_l])."""
    f = f_tensor(L)
    Ef = alt_cochain(3, lambda q: q[0][0] * q[0][1], f)
    Ef2 = Coch(3)
    for i in range(3):
        Ef2.add(2 * PV[i][0] * PV[i][1], f)
    same = Ef.flat() == Ef2.flat()
    dE = d_cochain(L, Ef)
    ff = t4.Engine(L).ff
    form = Coch(4)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        T = Tk(L.dim, ff.idx, ff.val)
        form.add(2 * (PV[i][0] * PV[j][1] + PV[i][1] * PV[j][0]) * (-1) ** (i + j), tperm(T, (i, j, k, l)))
    fl = dE.flat()
    return same, len(fl), fl == form.flat()


if __name__ == "__main__":
    for key in sys.argv[1:]:
        t0 = time.time()
        if key.startswith('sl') and key.endswith('tr'):
            L = lie.sl_tr(int(key[2:-2]))
            Ts = [f_tensor(L), d_tensor_sl(L)]
            label = "f, d = tr(X{Y,Z})"
        else:
            L = lie.build(key)
            Ts = [f_tensor(L)]
            label = "f (dim (g^{*(x)3})^g = 1, logs/a_classes.log)"
        txt, tab = run(L, Ts, label)
        print(txt, flush=True)
        same, n, okf = dEf_check(L)
        print("   E_f = alt(q0^1 q0^2 f) = 2(sum_i p_i^1 p_i^2) f: %s ; d E_f nonzero entries %d (d E_f != 0: %s) ; "
              "d E_f = 2 sum_{i<j} (-1)^{i+j} (p_i^1 p_j^2 + p_i^2 p_j^1) kappa([X_i,X_j],[X_k,X_l]): %s ; (%.1fs)"
              % (same, n, n > 0, okf, time.time() - t0), flush=True)
