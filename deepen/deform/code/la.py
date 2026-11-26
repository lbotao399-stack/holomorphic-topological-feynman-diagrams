"""Exact linear algebra over Q for cochain vectors.

rank_mod_p: Gaussian elimination modulo a prime (numpy int64, p < 2^31).
certified_rank: rank r from mod-p elimination, certified over Q by
   (a) an r x r minor with nonzero exact (Bareiss) determinant  -> rank >= r,
   (b) (#cols - r) exact kernel vectors, checked on all rows     -> rank <= r.
in_span: exact membership test v in col(M), with the exact solution.
"""
import numpy as np
from fractions import Fraction as Fr
from math import gcd

P1, P2 = 2147483629, 2147483587


def lcm(a, b):
    return a * b // gcd(a, b) if a and b else max(a, b)


def int_columns(cols):
    """list of columns (lists of Fractions) -> integer columns (each scaled by a positive rational)."""
    out = []
    for c in cols:
        L = 1
        for x in c:
            if x != 0:
                L = lcm(L, Fr(x).denominator)
        ci = [int(Fr(x) * L) for x in c]
        g = 0
        for x in ci:
            g = gcd(g, abs(x))
        if g > 1:
            ci = [x // g for x in ci]
        out.append(ci)
    return out


def matrix_from(cochains, keys=None):
    """cochains: list of dicts key -> Fraction (numeric N).  returns (keys, rows x cols list of Fractions)."""
    if keys is None:
        keys = sorted({k for c in cochains for k in c}, key=repr)
    idx = {k: i for i, k in enumerate(keys)}
    M = [[Fr(0)] * len(cochains) for _ in keys]
    for j, c in enumerate(cochains):
        for k, v in c.items():
            M[idx[k]][j] = Fr(v)
    return keys, M


def _modp(M, p):
    A = np.array([[int(x) % p for x in row] for row in M], dtype=np.int64)
    return A


def rref_mod_p(A, p):
    """returns pivot columns of A (numpy int64 mod p)."""
    A = A.copy() % p
    rows, cols = A.shape
    r = 0
    piv = []
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        i = r + nz[0]
        if i != r:
            A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr] = (A[nzr] - np.outer(col[nzr], A[r])) % p
        piv.append(c)
        r += 1
    return piv


def rank_mod_p(Mint, p):
    if len(Mint) == 0 or len(Mint[0]) == 0:
        return 0
    return len(rref_mod_p(_modp(Mint, p), p))


def bareiss_det(M):
    """exact determinant of a square integer matrix (list of lists)."""
    n = len(M)
    if n == 0:
        return 1
    A = [list(map(int, row)) for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def solve_square(A, b):
    """exact solution of A x = b, A square nonsingular (Fractions)."""
    n = len(A)
    M = [[Fr(x) for x in A[i]] + [Fr(b[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]


def _to_int_matrix(M):
    """rows x cols (Fractions/ints/numpy) -> (numpy object int matrix with integer columns, column scales)."""
    if isinstance(M, np.ndarray) and M.dtype != object:
        A = M.astype(object)
        return A, [Fr(1)] * M.shape[1]
    cols = list(zip(*M))
    Mi_cols = int_columns([list(c) for c in cols])
    scales = []
    for c, ci in zip(cols, Mi_cols):
        nz = next((t for t in range(len(c)) if Fr(c[t]) != 0), None)
        scales.append(Fr(ci[nz]) / Fr(c[nz]) if nz is not None else Fr(1))
    A = np.array(Mi_cols, dtype=object).T if Mi_cols else np.zeros((len(M), 0), dtype=object)
    return A, scales


def _as_int64(A):
    try:
        B = A.astype(np.int64)
        if A.dtype == object and not all(int(B.reshape(-1)[i]) == int(A.reshape(-1)[i]) for i in range(0, B.size, max(1, B.size // 1000))):
            return None
        return B
    except OverflowError:
        return None


def _modp_np(A, p):
    if A.dtype != object:
        return A % p
    B = None
    try:
        B = A.astype(np.int64)
        if np.abs(B).max() < 2 ** 62:
            return B % p
    except OverflowError:
        pass
    return np.array([[int(x) % p for x in row] for row in A], dtype=np.int64)


def _exact_matvec_zero(A, cols, x_int, extra_col, extra_coef):
    """checks A[:, cols] @ x_int + extra_coef * A[:, extra_col] == 0 exactly."""
    sub = A[:, list(cols) + [extra_col]]
    xv = [int(v) for v in x_int] + [int(extra_coef)]
    mx = max(abs(v) for v in xv)
    B = None
    try:
        B = sub.astype(np.int64)
        amax = int(np.abs(B).max()) if B.size else 0
    except OverflowError:
        B = None
    if B is not None and amax * mx * len(xv) < 2 ** 62:
        return not (B @ np.array(xv, dtype=np.int64)).any()
    s = sub.astype(object).dot(np.array(xv, dtype=object))
    return all(int(v) == 0 for v in s)


def certified_rank(M):
    """M: rows x cols (Fractions, ints or numpy int array).
    returns dict(rank, pivcols, pivrows, det, kernel, ranks_modp)."""
    if isinstance(M, np.ndarray):
        nrows, ncols = M.shape
    else:
        nrows = len(M)
        ncols = len(M[0]) if nrows else 0
    if nrows == 0 or ncols == 0:
        ker = [[Fr(int(i == j)) for i in range(ncols)] for j in range(ncols)]
        return dict(rank=0, pivcols=[], pivrows=[], det=1, kernel=ker, ranks_modp=[0, 0])
    A, scales = _to_int_matrix(M)
    res = []
    for p in (P1, P2):
        res.append(rref_mod_p(_modp_np(A, p), p))
    J = res[0] if len(res[0]) >= len(res[1]) else res[1]
    r = len(J)
    if r == 0:
        assert all(int(x) == 0 for x in A.reshape(-1))
        ker = [[Fr(int(i == j)) for i in range(ncols)] for j in range(ncols)]
        return dict(rank=0, pivcols=[], pivrows=[], det=1, kernel=ker, ranks_modp=[0, 0])
    AT = _modp_np(A[:, J].T.copy(), P1)
    I = rref_mod_p(AT, P1)
    assert len(I) == r
    sub = [[int(A[i, j]) for j in J] for i in I]
    det = bareiss_det(sub)
    assert det != 0, "minor vanishes over Q"
    ker = []
    for j in range(ncols):
        if j in J:
            continue
        x = solve_square(sub, [-int(A[i, j]) for i in I])
        D = 1
        for v in x:
            D = lcm(D, v.denominator)
        xi = [int(v * D) for v in x]
        assert _exact_matvec_zero(A, J, xi, j, D), "kernel check failed: rank over Q larger than mod-p rank"
        v = [Fr(0)] * ncols
        v[j] = Fr(1)
        for t in range(r):
            v[J[t]] = x[t]
        ker.append([v[c] * scales[c] for c in range(ncols)])
    return dict(rank=r, pivcols=J, pivrows=I, det=det, kernel=ker, ranks_modp=[len(x) for x in res])


def in_span(M, v, cert=None):
    """exact: is v in the column span of M?  returns (bool, x) with M x = v if True.
    M: rows x cols (lists or numpy int array); v: list of Fractions/ints."""
    if cert is None:
        cert = certified_rank(M)
    J, I = cert['pivcols'], cert['pivrows']
    if isinstance(M, np.ndarray):
        A = M
    else:
        A = None
    nrows = len(M)
    ncols = (M.shape[1] if A is not None else len(M[0]))
    if not J:
        return (all(Fr(x) == 0 for x in v), [Fr(0)] * ncols)
    get = (lambda i, j: int(A[i, j])) if A is not None else (lambda i, j: Fr(M[i][j]))
    sub = [[get(i, j) for j in J] for i in I]
    x = solve_square(sub, [Fr(v[i]) for i in I])
    # exact verification with integers: D (M_J x) == D v
    D = 1
    for t in list(x) + [Fr(y) for y in v]:
        D = lcm(D, Fr(t).denominator)
    xi = np.array([int(t * D) for t in x], dtype=object)
    vi = np.array([int(Fr(y) * D) for y in v], dtype=object)
    if A is not None:
        lhs = A[:, J].astype(object).dot(xi)
    else:
        Ms = []
        for i in range(nrows):
            row = [Fr(M[i][j]) for j in J]
            Ms.append(row)
        L = 1
        for row in Ms:
            for t in row:
                L = lcm(L, t.denominator)
        lhs_rows = np.array([[int(t * L) for t in row] for row in Ms], dtype=object)
        lhs = lhs_rows.dot(xi)
        vi = vi * L
    if not all(int(a_) == int(b_) for a_, b_ in zip(lhs, vi)):
        return (False, None)
    full = [Fr(0)] * ncols
    for t, j in enumerate(J):
        full[j] = x[t]
    return (True, full)


def column(c, keys):
    return [Fr(c.get(k, 0)) for k in keys]
