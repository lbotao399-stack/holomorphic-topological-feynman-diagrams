"""Integer-power wheel periods W_n (n>=3; W_3 = K_4) by the radial-angular method.

Pure-Fraction port of explore/nonSP/verify/full_ward/wheel_residues.py:
spoke powers p_i (hub - rim vertex i), rim powers r_i (rim vertex i - i+1), projective
period  int prod t^(a-1) / U^q,  sum a = n q.  Output {weight: Fraction}: weight 0 is the
rational part, weight w>=3 odd multiplies zeta(w).
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial, prod


# ---- univariate polynomials as tuples of Fractions (low degree first) ----
def pmul(a, b):
    if not a or not b:
        return ()
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    return tuple(c)


def ppow(a, n):
    r = (F(1),)
    for _ in range(n):
        r = pmul(r, a)
    return r


def pshift(a, h):
    """a(X + h)"""
    out = ()
    for k, x in reversed(list(enumerate(a))):
        out = pmul(out, (F(h), F(1))) if out else ()
        out = tuple(list(out) + [F(0)] * max(0, 1 - len(out)))
        out = (out[0] + x,) + out[1:] if out else (x,)
    return out


@lru_cache(maxsize=None)
def multiplicity(q):
    """M(s) = 2 s^2 prod_{j=1}^{lam-1} (s^2 - j^2) / (2 lam)!,  lam = q-1, as polynomial in s."""
    lam = q - 1
    p = (F(0), F(0), F(2))
    for j in range(1, lam):
        p = pmul(p, (F(-j * j), F(0), F(1)))
    return tuple(x / factorial(2 * lam) for x in p)


def indices(spoke, rim, q):
    n = len(rim)
    if len(spoke) != n or n < 3 or sum(spoke) + sum(rim) != n * q:
        raise ValueError('Expected n spoke and n rim powers of total n*q.')
    if any(not 1 <= v < q for v in (*spoke, *rim)) or sum(rim) <= q:
        raise ValueError('Divergent projective wheel powers.')
    c = 0
    nodes = []
    for i in range(n):
        c += 2 * spoke[i] + rim[i - 1] + rim[i] - 2 * q
        d = rim[i] - 1
        nodes.extend(range(c - d, c + d + 1, 2))
    assert c == 0
    minimum = min(nodes)
    assert all((v - minimum) % 2 == 0 for v in nodes)
    counts = Counter((v - minimum) // 2 for v in nodes)
    length = max(counts) + 1
    if length >= q:
        raise ValueError('Divergent wheel radial gap.')
    histogram = tuple(counts[j] for j in range(length))
    factor = prod(factorial(q - p - 1) * factorial(r - 1) for p, r in zip(spoke, rim))
    return min(histogram, histogram[::-1]), factor


@lru_cache(maxsize=None)
def summation(q, h, u):
    """sum_{s >= q-1} M(s)/(s-h)^u, polynomial parts dropped (they cancel in the total)."""
    coefficients = pshift(multiplicity(q), h)
    coefficients = list(coefficients) + [F(0)] * (u + 1)
    zeta = [F(0) for _ in range(u + 1)]
    for j in range(u):
        weight = u - j
        zeta[weight] = coefficients[j]
        zeta[0] -= coefficients[j] * sum((F(1, k ** weight) for k in range(1, q - 1 - h)), F(0))
    return tuple(zeta)


@lru_cache(maxsize=None)
def period(q, histogram):
    total = sum(histogram)
    parts = []
    for a, n in enumerate(histogram):
        if not n:
            continue
        den = (F(1),)
        for j, m in enumerate(histogram):
            if j != a and m:
                den = pmul(den, ppow((F(a - j), F(1)), m))
        coeff = []
        for i in range(n):
            value = (F(int(i == 0)) - sum(den[j] * coeff[i - j] for j in range(1, min(i, len(den) - 1) + 1))) / den[0]
            coeff.append(value)
            parts.append((a, n - i, value))
    pairs = Counter()
    for a, k, x in parts:
        for b, l, y in parts:
            pairs[a - b, k + l - 1] += (-1) ** (total + l) * comb(k + l - 2, k - 1) * x * y
    out = Counter()
    for (h, u), value in pairs.items():
        if value:
            for w, c in enumerate(summation(q, h, u)):
                if c:
                    out[w] += c * value
    assert all(not v or w == 0 or (w >= 3 and w % 2) for w, v in out.items()), out
    return tuple((w, v) for w, v in sorted(out.items()) if v)


def integrate(spoke, rim, q):
    histogram, factor = indices(tuple(spoke), tuple(rim), q)
    values = {w: factor * c for w, c in period(q, histogram)}
    assert all(w <= 2 * len(spoke) - 3 for w in values)
    return values


def radial_master(n, q):
    """binom(2n-2,n-1) sum_{s>=1} M(s)/s^(2n-1)  (spokes q-1, rims 1)"""
    M = multiplicity(q)
    out = Counter()
    for j, c in enumerate(M):
        if c:
            w = 2 * n - 1 - j
            assert w >= 2
            out[w] += c * comb(2 * n - 2, n - 1)
    return dict(out)


def as_label(values):
    lab = {}
    for w, v in values.items():
        lab['1' if w == 0 else 'z%d' % w] = v
    return lab


if __name__ == '__main__':
    import time
    t = time.time()
    for n in range(3, 11):
        for q in range(2, n):
            a = integrate([q - 1] * n, [1] * n, q)
            b = radial_master(n, q)
            assert a == b, (n, q, a, b)
    print('radial masters n=3..10 agree', round(time.time() - t, 2), 's')
    v = integrate([5, 5, 4, 5, 5], [3, 2, 2, 2, 2], 7)
    ref = {0: F(-45001, 43545600), 3: F(-116849, 29937600), 5: F(50837, 2993760), 7: F(-7, 594)}
    print('five-loop example', v, v == ref)
    print('WS_3 = K4 at q=2:', integrate([1, 1, 1], [1, 1, 1], 2))
    print('WS_4 at q=2:', integrate([1] * 4, [1] * 4, 2))
    print('W4 radial q=3:', integrate([2] * 4, [1] * 4, 3))
