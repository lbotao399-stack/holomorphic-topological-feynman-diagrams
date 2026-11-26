"""Closed form of the wheel radial masters (spokes q-1, rims 1):
I_q(W_n) = binom(2n-2,n-1) * 2/(2q-2)! * sum_{k=1}^{q-1} c_{q,k} zeta(2n-1-2k),
prod_{j=1}^{q-2} (x - j^2) = sum_k c_{q,k} x^(k-1); checked against wheel.integrate."""
from fractions import Fraction as F
from math import comb, factorial

import wheel


def closed(n, q):
    c = [F(1)]
    for j in range(1, q - 1):
        c = [F(0)] + c
        for k in range(len(c) - 1):
            c[k] -= j * j * c[k + 1]
    out = {}
    for k, ck in enumerate(c, start=1):
        if ck:
            w = 2 * n - 1 - 2 * k
            out[w] = out.get(w, 0) + comb(2 * n - 2, n - 1) * F(2, factorial(2 * q - 2)) * ck
    return {w: v for w, v in out.items() if v}


ok = 0
for n in range(3, 13):
    for q in range(2, n):
        assert closed(n, q) == wheel.integrate([q - 1] * n, [1] * n, q), (n, q)
        ok += 1
print('closed form = residue algorithm for %d (n,q), 3<=n<=12, 2<=q<n' % ok)
print('n=5:', {q: closed(5, q) for q in (2, 3, 4)})
