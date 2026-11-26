"""Exact linear algebra over Q (fractions.Fraction) and a small exterior algebra.

Matrices are lists of rows.  Differential forms on Schwinger space are dicts
{bitmask: coefficient}; bit e of the mask stands for dt_e and the monomial of a mask is the
wedge of its dt_e in increasing order of e.
"""
from fractions import Fraction as Fr

__all__ = [
    "Fr", "fmat", "fvec", "det", "solve", "inverse", "matmul", "transpose", "identity",
    "perm_sign", "shuffle_sign", "wedge", "wedge_all", "one_form", "form_scale", "form_add",
    "pfaffian_forms", "basic_coefficient",
]


def fvec(v):
    return [Fr(x) for x in v]


def fmat(rows):
    return [[Fr(x) for x in r] for r in rows]


def identity(n):
    return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(M):
    return [list(r) for r in zip(*M)] if M else []


def matmul(A, B):
    Bt = transpose(B)
    return [[sum((a * b for a, b in zip(row, col)), Fr(0)) for col in Bt] for row in A]


def det(M):
    """Determinant by Gaussian elimination over Q.  det of the 0 x 0 matrix is 1."""
    M = [[Fr(x) for x in r] for r in M]
    n = len(M)
    d = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        piv = M[i][i]
        d *= piv
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / piv
                Mi, Mr = M[i], M[r]
                for c in range(i, n):
                    Mr[c] -= f * Mi[c]
    return d


def solve(A, B):
    """Solve A X = B over Q.  A: n x n, B: n x k (list of rows).  Returns X as n x k."""
    n = len(A)
    k = len(B[0]) if n else 0
    aug = [[Fr(x) for x in A[i]] + [Fr(x) for x in B[i]] for i in range(n)]
    for i in range(n):
        p = next((r for r in range(i, n) if aug[r][i] != 0), None)
        if p is None:
            raise ZeroDivisionError("singular matrix")
        aug[i], aug[p] = aug[p], aug[i]
        piv = aug[i][i]
        aug[i] = [x / piv for x in aug[i]]
        for r in range(n):
            if r != i and aug[r][i] != 0:
                f = aug[r][i]
                ai, ar = aug[i], aug[r]
                for c in range(i, n + k):
                    ar[c] -= f * ai[c]
    return [row[n:] for row in aug]


def inverse(A):
    return solve(A, identity(len(A)))


def perm_sign(seq):
    """Sign of the permutation i -> seq[i] of 0..len-1."""
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        while seq[i] != i:
            j = seq[i]
            seq[i], seq[j] = seq[j], seq[i]
            s = -s
    return s


def shuffle_sign(first, n):
    """Sign of the permutation (first..., rest...) of 0..n-1, both parts in increasing order."""
    first = sorted(first)
    rest = [i for i in range(n) if i not in set(first)]
    order = first + rest
    inv = [0] * n
    for pos, i in enumerate(order):
        inv[i] = pos
    return perm_sign(inv)


# ---------------------------------------------------------------- exterior algebra

def one_form(coeffs):
    """sum_e coeffs[e] dt_e."""
    return {1 << e: Fr(c) for e, c in enumerate(coeffs) if c != 0}


def form_scale(a, c):
    c = Fr(c)
    return {m: c * v for m, v in a.items()} if c != 0 else {}


def form_add(a, b):
    out = dict(a)
    for m, v in b.items():
        w = out.get(m, 0) + v
        if w == 0:
            out.pop(m, None)
        else:
            out[m] = w
    return out


def _wedge_sign(a, b):
    """Sign of sorting (indices of a, indices of b); a and b disjoint masks."""
    s = 0
    bb = b
    while bb:
        low = bb & -bb
        j = low.bit_length() - 1
        s += bin(a >> (j + 1)).count("1")
        bb ^= low
    return -1 if s & 1 else 1


def wedge(a, b):
    out = {}
    for ma, va in a.items():
        for mb, vb in b.items():
            if ma & mb:
                continue
            m = ma | mb
            v = va * vb
            if _wedge_sign(ma, mb) < 0:
                v = -v
            w = out.get(m, 0) + v
            if w == 0:
                out.pop(m, None)
            else:
                out[m] = w
    return out


def wedge_all(forms):
    out = {0: Fr(1)}
    for f in forms:
        out = wedge(out, f)
        if not out:
            return {}
    return out


def pfaffian_forms(A, idx=None):
    """Pfaffian of an antisymmetric matrix A[r][s] of even forms (dicts), restricted to the index
    list idx:  Pf = sum_j (-1)^j A[i0][i_j] Pf(without i0, i_j),  j = 1..len-1."""
    if idx is None:
        idx = list(range(len(A)))
    if len(idx) == 0:
        return {0: Fr(1)}
    if len(idx) % 2:
        return {}
    i0 = idx[0]
    out = {}
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        term = wedge(A[i0][idx[j]], pfaffian_forms(A, rest))
        if (j - 1) % 2:
            term = form_scale(term, -1)
        out = form_add(out, term)
    return out


def basic_coefficient(form, t, e0=None):
    """For a basic form  w = D * iota_E(dt_0 ^ ... ^ dt_{N-1})  of degree N-1 (E = sum t_e d/dt_e),
    return D = (-1)^{e0} w[all but e0] / t[e0].  With e0 = None every e0 is used and compared."""
    N = len(t)
    full = (1 << N) - 1
    if e0 is not None:
        c = form.get(full ^ (1 << e0), Fr(0))
        return (-c if e0 & 1 else c) / t[e0]
    vals = []
    for e in range(N):
        c = form.get(full ^ (1 << e), Fr(0))
        vals.append((-c if e & 1 else c) / t[e])
    if any(v != vals[0] for v in vals):
        raise ArithmeticError("form is not basic: coefficients %s" % vals)
    return vals[0]
