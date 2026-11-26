"""numba kernel: P(t) = det[C_re Y^i_e]_{(i,r), e != 0} / (t_0 U^(n(H-1))) mod p at many points
(same algebra as numerator.BareNumerator.samples_mod; Y = U y, M ell = -C t y0)."""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def _inv(a, p):
    r = 1
    b = a % p
    e = p - 2
    while e:
        if e & 1:
            r = r * b % p
        b = b * b % p
        e >>= 1
    return r


@nb.njit(cache=True)
def _det_inplace(A, k, p):
    det = 1
    for c in range(k):
        piv = -1
        for r in range(c, k):
            if A[r, c] != 0:
                piv = r
                break
        if piv < 0:
            return 0
        if piv != c:
            for j in range(c, k):
                tmp = A[c, j]
                A[c, j] = A[piv, j]
                A[piv, j] = tmp
            det = p - det if det else 0
        d = A[c, c]
        det = det * d % p
        inv = _inv(d, p)
        for r in range(c + 1, k):
            f = A[r, c] * inv % p
            if f:
                for j in range(c + 1, k):
                    A[r, j] = (A[r, j] - f * A[c, j]) % p
    return det


@nb.njit(parallel=True, cache=True)
def samples(tpts, C, y0, n, H, p):
    B, N = tpts.shape
    out = np.zeros(B, dtype=np.int64)
    k = n * H
    for b in nb.prange(B):
        t = tpts[b]
        # M and rhs
        Mx = np.zeros((n, n + H), dtype=np.int64)
        for r in range(n):
            for s in range(n):
                acc = 0
                for e in range(N):
                    if C[r, e] != 0 and C[s, e] != 0:
                        acc += C[r, e] * C[s, e] * t[e]
                Mx[r, s] = acc % p
            for i in range(H):
                acc = 0
                for e in range(N):
                    if C[r, e] != 0 and y0[e, i] != 0:
                        acc = (acc - C[r, e] * (t[e] % p) % p * y0[e, i]) % p
                Mx[r, n + i] = acc % p
        # Gauss-Jordan mod p with determinant
        U = 1
        for c in range(n):
            piv = -1
            for r in range(c, n):
                if Mx[r, c] != 0:
                    piv = r
                    break
            if piv < 0:
                U = 0
                break
            if piv != c:
                for j in range(n + H):
                    tmp = Mx[c, j]
                    Mx[c, j] = Mx[piv, j]
                    Mx[piv, j] = tmp
                U = p - U if U else 0
            d = Mx[c, c]
            U = U * d % p
            inv = _inv(d, p)
            for j in range(n + H):
                Mx[c, j] = Mx[c, j] * inv % p
            for r in range(n):
                if r != c:
                    f = Mx[r, c]
                    if f:
                        for j in range(n + H):
                            Mx[r, j] = (Mx[r, j] - f * Mx[c, j]) % p
        if U == 0:
            out[b] = -1
            continue
        # Y = U y0 + C^T (U ell)
        Y = np.zeros((N, H), dtype=np.int64)
        for e in range(N):
            for i in range(H):
                acc = U * (y0[e, i] % p) % p
                for r in range(n):
                    if C[r, e] != 0:
                        acc = (acc + C[r, e] * (Mx[r, n + i] * U % p)) % p
                Y[e, i] = acc % p
        A = np.zeros((k, k), dtype=np.int64)
        for i in range(H):
            for r in range(n):
                for e in range(1, N):
                    A[i * n + r, e - 1] = (C[r, e] * Y[e, i]) % p
        d = _det_inplace(A, k, p)
        den = t[0] % p
        for _ in range(n * (H - 1)):
            den = den * U % p
        out[b] = d * _inv(den, p) % p
    return out


@nb.njit(parallel=True, cache=True)
def samples_ins(tpts, C, y0, n, H, p, ins_e, ins_i):
    """P(t) * prod_k (U y)^{ins_i[k]}_{ins_e[k]} mod p  (direct whole-graph route, for checks)"""
    B, N = tpts.shape
    out = np.zeros(B, dtype=np.int64)
    k = n * H
    m = len(ins_e)
    for b in nb.prange(B):
        t = tpts[b]
        Mx = np.zeros((n, n + H), dtype=np.int64)
        for r in range(n):
            for s in range(n):
                acc = 0
                for e in range(N):
                    if C[r, e] != 0 and C[s, e] != 0:
                        acc += C[r, e] * C[s, e] * t[e]
                Mx[r, s] = acc % p
            for i in range(H):
                acc = 0
                for e in range(N):
                    if C[r, e] != 0 and y0[e, i] != 0:
                        acc = (acc - C[r, e] * (t[e] % p) % p * y0[e, i]) % p
                Mx[r, n + i] = acc % p
        U = 1
        for c in range(n):
            piv = -1
            for r in range(c, n):
                if Mx[r, c] != 0:
                    piv = r
                    break
            if piv < 0:
                U = 0
                break
            if piv != c:
                for j in range(n + H):
                    tmp = Mx[c, j]
                    Mx[c, j] = Mx[piv, j]
                    Mx[piv, j] = tmp
                U = p - U if U else 0
            d = Mx[c, c]
            U = U * d % p
            inv = _inv(d, p)
            for j in range(n + H):
                Mx[c, j] = Mx[c, j] * inv % p
            for r in range(n):
                if r != c:
                    f = Mx[r, c]
                    if f:
                        for j in range(n + H):
                            Mx[r, j] = (Mx[r, j] - f * Mx[c, j]) % p
        if U == 0:
            out[b] = -1
            continue
        Y = np.zeros((N, H), dtype=np.int64)
        for e in range(N):
            for i in range(H):
                acc = U * (y0[e, i] % p) % p
                for r in range(n):
                    if C[r, e] != 0:
                        acc = (acc + C[r, e] * (Mx[r, n + i] * U % p)) % p
                Y[e, i] = acc % p
        A = np.zeros((k, k), dtype=np.int64)
        for i in range(H):
            for r in range(n):
                for e in range(1, N):
                    A[i * n + r, e - 1] = (C[r, e] * Y[e, i]) % p
        d = _det_inplace(A, k, p)
        den = t[0] % p
        for _ in range(n * (H - 1)):
            den = den * U % p
        val = d * _inv(den, p) % p
        for q in range(m):
            val = val * Y[ins_e[q], ins_i[q]] % p
        out[b] = val
    return out
