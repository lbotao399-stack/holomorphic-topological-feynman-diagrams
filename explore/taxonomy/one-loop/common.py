"""Shared exact helpers for explore/taxonomy/one-loop (no project file is modified)."""
from fractions import Fraction as F
from itertools import permutations, combinations
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = HERE.parents[2]


class Log:
    def __init__(self, name):
        self.path = HERE / name
        self.lines = []

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.lines.append(s)

    def close(self):
        self.path.write_text("\n".join(self.lines) + "\n")


def det(M):
    """exact determinant (Fractions or ints), Gaussian elimination"""
    A = [[F(x) for x in row] for row in M]
    n = len(A)
    s = F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            A[c], A[p] = A[p], A[c]
            s = -s
        s *= A[c][c]
        inv = 1 / A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] * inv
            if f:
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return s


def inverse(M):
    n = len(M)
    A = [[F(x) for x in row] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def perm_sign(p):
    p = list(p)
    s = 1
    seen = [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def seq_sign(seq):
    """sign of the permutation sorting a sequence of distinct keys"""
    s = 1
    a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]:
                s = -s
    return s


# ------------------------------------------------------------------ exterior algebra on dt_0..dt_{N-1}
# a form is a dict {bitmask: coefficient}; degree = popcount

def _wsign(a, b):
    """sign of e_a ^ e_b -> e_{a|b} for disjoint masks (number of pairs i in a, j in b with i > j)"""
    s = 0
    x = a
    while x:
        lowb = x & -x
        i = lowb.bit_length() - 1
        s += bin(b & ((1 << i) - 1)).count("1")
        x ^= lowb
    return -1 if s & 1 else 1


def wedge(f, g, mod=None):
    out = {}
    for a, x in f.items():
        for b, y in g.items():
            if a & b:
                continue
            c = a | b
            v = _wsign(a, b) * x * y
            if mod:
                v %= mod
            out[c] = out.get(c, 0) + v
    if mod:
        return {k: v % mod for k, v in out.items() if v % mod}
    return {k: v for k, v in out.items() if v != 0}


def fadd(f, g, c=1, mod=None):
    out = dict(f)
    for k, v in g.items():
        out[k] = out.get(k, 0) + c * v
    if mod:
        return {k: v % mod for k, v in out.items() if v % mod}
    return {k: v for k, v in out.items() if v != 0}


def fscale(f, c, mod=None):
    if mod:
        return {k: (v * c) % mod for k, v in f.items() if (v * c) % mod}
    return {k: v * c for k, v in f.items() if v * c != 0}


def one_form(coeffs):
    return {1 << e: c for e, c in enumerate(coeffs) if c != 0}


def pfaffian(A, mod=None, idx=None, memo=None):
    """Pfaffian of an antisymmetric matrix whose entries are even forms (dicts); recursive expansion
    along the first row with memoization on index tuples."""
    if idx is None:
        idx = tuple(range(len(A)))
        memo = {}
    if len(idx) == 0:
        return {0: 1}
    if len(idx) % 2:
        return {}
    if idx in memo:
        return memo[idx]
    i = idx[0]
    out = {}
    for k, j in enumerate(idx[1:]):
        rest = tuple(x for x in idx[1:] if x != j)
        sub = pfaffian(A, mod, rest, memo)
        if not sub or not A[i][j]:
            continue
        term = wedge(A[i][j], sub, mod)
        out = fadd(out, term, -1 if k % 2 else 1, mod)
    memo[idx] = out
    return out
