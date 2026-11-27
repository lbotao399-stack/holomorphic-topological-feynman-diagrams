"""Shared helpers for explore/taxonomy/cs-class (logging, permutations, minors, mod-p linear algebra)."""
import os
import sys
import time
from fractions import Fraction as Fr
from itertools import permutations, combinations

P_ROOT = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + ''
HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (P_ROOT + '/code', P_ROOT + '/explore/twoloop-universal', P_ROOT + '/deepen/simple-all/code'):
    if _p not in sys.path:
        sys.path.insert(0, _p)


class Log:
    """print to stdout and to <script>.log next to the script."""

    def __init__(self, script):
        name = os.path.splitext(os.path.basename(script))[0] + os.environ.get('LOG_SUFFIX', '')
        self.path = os.path.join(HERE, name + '.log')
        self.f = open(self.path, 'w')
        self.t0 = time.time()
        self('# %s  %s' % (name, time.strftime('%Y-%m-%d %H:%M:%S')))

    def __call__(self, *a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(s + '\n')
        self.f.flush()

    def close(self):
        self('# done in %.1f s' % (time.time() - self.t0))
        self.f.close()


def parity(seq):
    """sign of the permutation given as a sequence of distinct comparable items (inversion count)."""
    s = list(seq)
    inv = 0
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            if s[i] > s[j]:
                inv += 1
    return -1 if inv % 2 else 1


def perm_from_lists(src, tgt):
    """sign of the permutation taking the list src to the list tgt (same elements)."""
    pos = {x: i for i, x in enumerate(tgt)}
    return parity([pos[x] for x in src])


def det_int(M):
    """exact determinant of a small integer/Fraction matrix (Bareiss for ints, else Fractions)."""
    n = len(M)
    if n == 0:
        return 1
    if all(isinstance(x, int) for r in M for x in r):
        A = [list(r) for r in M]
        sign, prev = 1, 1
        for k in range(n - 1):
            if A[k][k] == 0:
                p = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
                if p is None:
                    return 0
                A[k], A[p] = A[p], A[k]
                sign = -sign
            for i in range(k + 1, n):
                for j in range(k + 1, n):
                    A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            prev = A[k][k]
        return sign * A[n - 1][n - 1]
    A = [[Fr(x) for x in r] for r in M]
    d = Fr(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None:
            return Fr(0)
        if p != i:
            A[i], A[p] = A[p], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            if A[r][i] != 0:
                f = A[r][i] / A[i][i]
                for c in range(i, n):
                    A[r][c] -= f * A[i][c]
    return d


# ------------------------------------------------------------------ mod-p linear algebra
def rank_mod(rows, p):
    """rank over F_p of a list of integer rows (lists)."""
    M = [[x % p for x in r] for r in rows]
    if not M:
        return 0
    ncol = len(M[0])
    r = 0
    for c in range(ncol):
        piv = next((i for i in range(r, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                Mi, Mr = M[i], M[r]
                M[i] = [(a - f * b) % p for a, b in zip(Mi, Mr)]
        r += 1
        if r == len(M):
            break
    return r


def np_rank_mod(A, p):
    """rank over F_p of a numpy int64 matrix (p < 2^31), Gaussian elimination with numpy rows."""
    import numpy as np
    M = np.array(A, dtype=np.int64) % p
    nr, nc = M.shape
    r = 0
    for c in range(nc):
        if r == nr:
            break
        col = M[r:, c]
        nz = np.nonzero(col)[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        f = M[:, c].copy()
        f[r] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(f[nzr], M[r]) % p) % p
        r += 1
    return r


def np_solve_mod(A, b, p):
    """solve A x = b over F_p (A numpy int64, b vector); returns (consistent, x or None, rank)."""
    import numpy as np
    A = np.array(A, dtype=np.int64) % p
    b = np.array(b, dtype=np.int64).reshape(-1, 1) % p
    M = np.concatenate([A, b], axis=1)
    nr, nc = M.shape
    r = 0
    pivcols = []
    for c in range(nc - 1):
        if r == nr:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        piv = r + nz[0]
        if piv != r:
            M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        f = M[:, c].copy()
        f[r] = 0
        nzr = np.nonzero(f)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(f[nzr], M[r]) % p) % p
        pivcols.append(c)
        r += 1
    consistent = not np.any(M[r:, -1])
    if not consistent:
        return False, None, r
    x = np.zeros(nc - 1, dtype=np.int64)
    for i, c in enumerate(pivcols):
        x[c] = M[i, -1]
    return True, x, r


def ratrec(a, m):
    """rational reconstruction of a mod m."""
    a %= m
    if a == 0:
        return Fr(0)
    r0, r1, s0, s1 = m, a, 0, 1
    import math
    bound = math.isqrt(m // 2)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def crt_pair(a1, m1, a2, m2):
    import math
    g = math.gcd(m1, m2)
    assert g == 1
    t = ((a2 - a1) * pow(m1, -1, m2)) % m2
    return (a1 + m1 * t) % (m1 * m2), m1 * m2


PRIMES = [2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483543]
