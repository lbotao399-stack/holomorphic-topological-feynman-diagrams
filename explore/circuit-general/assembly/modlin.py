"""Batched modular linear algebra (int64, p < 2^31), CRT and rational reconstruction.
det_mod, solve_mod, modinv_vec: copied from cores/fastker.py; ratrec: cores/modsolve.py."""
from fractions import Fraction as F

import numpy as np

PRIMES = [2147483629, 2147483587, 2147483579, 2147483563, 2147483549, 2147483543, 2147483497,
          2147483489, 2147483477, 2147483423, 2147483399, 2147483353, 2147483323, 2147483269,
          2147483249, 2147483237, 2147483179, 2147483171, 2147483137, 2147483123]


def _is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1 if i == 2 else 2
    return True


assert all(_is_prime(p) for p in PRIMES)


def modinv_vec(a, p):
    r = np.ones_like(a)
    b = a % p
    e = p - 2
    while e:
        if e & 1:
            r = r * b % p
        b = b * b % p
        e >>= 1
    return r


def det_mod(A, p):
    """batched determinant mod p of A (B, k, k) int64; destroys A"""
    B, k, _ = A.shape
    det = np.ones(B, dtype=np.int64)
    idx = np.arange(B)
    for c in range(k):
        col = A[:, c:, c]
        nz = col != 0
        has = nz.any(axis=1)
        piv = c + np.argmax(nz, axis=1)
        det[~has] = 0
        swap = piv != c
        if swap.any():
            rows_c = A[idx, c, :].copy()
            A[idx, c, :] = A[idx, piv, :]
            A[idx, piv, :] = rows_c
            det[swap] = (-det[swap]) % p
        d = A[:, c, c] % p
        det = det * d % p
        inv = modinv_vec(np.where(d == 0, 1, d), p)
        if c + 1 < k:
            f = A[:, c + 1:, c] * inv[:, None] % p
            A[:, c + 1:, :] = (A[:, c + 1:, :] - f[:, :, None] * A[:, c, None, :]) % p
    return det % p


def solve_mod(Mm, R, p):
    """batched solve Mm x = R mod p (Mm: (B,n,n), R: (B,n,h)); returns x and det Mm"""
    B, n, _ = Mm.shape
    A = np.concatenate([Mm, R], axis=2) % p
    det = np.ones(B, dtype=np.int64)
    idx = np.arange(B)
    for c in range(n):
        nz = A[:, c:, c] != 0
        piv = c + np.argmax(nz, axis=1)
        swap = piv != c
        if swap.any():
            rows_c = A[idx, c, :].copy()
            A[idx, c, :] = A[idx, piv, :]
            A[idx, piv, :] = rows_c
            det[swap] = (-det[swap]) % p
        d = A[:, c, c] % p
        det = det * d % p
        inv = modinv_vec(d, p)
        A[:, c, :] = A[:, c, :] * inv[:, None] % p
        for r in range(n):
            if r == c:
                continue
            f = A[:, r, c].copy()
            A[:, r, :] = (A[:, r, :] - f[:, None] * A[:, c, :]) % p
    return A[:, :, n:], det


def fmod(x, p):
    """Fraction -> residue mod p"""
    x = F(x)
    return x.numerator % p * pow(x.denominator % p, p - 2, p) % p


def ratrec(a, m):
    """rational reconstruction of a mod m (Wang)"""
    a %= m
    if a == 0:
        return F(0)
    bound = int((m // 2) ** 0.5)
    while bound * bound > m // 2:
        bound -= 1
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return F(r1, s1) if s1 > 0 else F(-r1, -s1)


class CRT:
    """accumulates residues of many keyed values over primes; reconstructs when stable"""

    def __init__(self):
        self.M = 1
        self.acc = {}
        self.last = None
        self.nprimes = 0

    def add(self, p, residues):
        """residues: dict key -> int mod p (missing keys = 0)"""
        keys = set(self.acc) | set(residues)
        inv = pow(self.M % p, p - 2, p)
        for k in keys:
            a = self.acc.get(k, 0)
            x = residues.get(k, 0)
            self.acc[k] = a + self.M * ((x - a) * inv % p)
        self.M *= p
        self.nprimes += 1

    def reconstruct(self):
        out = {}
        for k, a in self.acc.items():
            r = ratrec(a, self.M)
            if r is None:
                return None
            if r:
                out[k] = r
        return out

    def stable(self, min_primes=2):
        """returns reconstruction if it equals the previous one (and >= min_primes used)"""
        rec = self.reconstruct()
        ok = rec is not None and rec == self.last and self.nprimes >= min_primes
        self.last = rec
        return rec if ok else None


def dotmod(PF, V, p, check=True):
    """(PF @ V) mod p for int64 residues PF (G x M), V (M x L), M < 2^22: exact float64 BLAS on 16-bit limbs
    (products < 2^31, sums < 2^53); one entry re-checked in int64 arithmetic"""
    with np.errstate(all='ignore'):
        a1, a0 = (PF >> 16).astype(np.float64), (PF & 0xFFFF).astype(np.float64)
        b1, b0 = (V >> 16).astype(np.float64), (V & 0xFFFF).astype(np.float64)
        r11 = (a1 @ b1).astype(np.int64) % p
        r10 = ((a1 @ b0).astype(np.int64) + (a0 @ b1).astype(np.int64)) % p
        r00 = (a0 @ b0).astype(np.int64) % p
    s16, s32 = pow(2, 16, p), pow(2, 32, p)
    R = (r11 * s32 % p + r10 * s16 % p + r00) % p
    if check and PF.size and V.size:
        g = int(np.random.randint(PF.shape[0]))
        l = int(np.random.randint(V.shape[1]))
        if int((PF[g] * V[:, l] % p).sum() % p) != int(R[g, l]):
            raise ArithmeticError('float64 limb product not exact')
    return R
