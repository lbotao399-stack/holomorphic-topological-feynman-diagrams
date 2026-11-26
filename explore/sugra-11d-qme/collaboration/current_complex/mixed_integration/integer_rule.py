"""Integral-node Newton cubature for the third 7-adic pole, with no degree-eight hypothesis."""

from collections import Counter
from functools import cache
from fractions import Fraction as F
from math import factorial, comb, prod
from cubature import compositions
from pole_weights import value

PRIME = 7
MAX = 11


def vp(n):
    if not n:
        return 1000
    a = 0
    while n % 7 == 0:
        n //= 7
        a += 1
    return a


def val(v):
    return vp(v.numerator) - vp(v.denominator)


@cache
def falling_binomial(n):
    if n == 0:
        return (F(1),)
    p = falling_binomial(n - 1)
    return tuple(
        ((p[k - 1] if k else 0) - (n - 1) * (p[k] if k < n else 0)) / n
        for k in range(n + 1)
    )


@cache
def newton(alpha):
    p = [F(1)]
    for n in alpha:
        f = [factorial(i) * x for i, x in enumerate(falling_binomial(n))]
        out = [F(0)] * (len(p) + n)
        for i, x in enumerate(p):
            for j, y in enumerate(f):
                out[i + j] += x * y
        p = out
    return sum(x * factorial(5) / factorial(5 + i) for i, x in enumerate(p))


@cache
def layer(x, k):
    d = sum(x)
    if k < d:
        return F()
    s = F()
    for a in compositions(k - d, 5):
        u = tuple(v + j for v, j in zip(x, a))
        s += newton(tuple(sorted(u))) * prod(comb(i, j) for i, j in zip(u, x))
    return (-1) ** (k - d) * s


@cache
def weight(x, y):
    return sum(
        layer(x, r) * layer(y, s)
        for r in range(sum(x), MAX + 1)
        for s in range(sum(y), MAX + 1 - r)
    )


def angular_weight():
    # The degree bound for conductance alone is 10, since its divergent
    # current direction is common to all three paths.
    kmax = 10
    rows = Counter()
    for i in range(kmax + 1):
        for j in range(kmax + 1 - i):
            base = sum(
                u * v * value(a, b)
                for a, u in enumerate(falling_binomial(i))
                for b, v in enumerate(falling_binomial(j))
            )
            for x in range(i + 1):
                for y in range(j + 1):
                    rows[x, y] += (
                        base * (-1) ** (i + j - x - y) * comb(i, x) * comb(j, y)
                    )
    return rows


def residue(a, mod):
    return a.numerator * pow(a.denominator, -1, mod) % mod


def main():
    domains = [list(compositions(d, 5)) for d in range(MAX + 1)]
    reps = [[x for x in v if tuple(sorted(x)) == x] for v in domains]
    vol = Counter()
    for d in range(MAX + 1):
        for e in range(MAX + 1 - d):
            for x in reps[d]:
                for y in reps[e]:
                    v = weight(x, y)

                    def orbit(t):
                        return factorial(5) // prod(
                            factorial(n) for n in Counter(t).values()
                        )

                    vol[val(v)] += orbit(x) * orbit(y)
    ang = Counter(val(x) for x in angular_weight().values())
    print("resistance weight valuations", vol, "conductance weight valuations", ang)
    print(
        "nonzero tensor pairs at normalized 7^3",
        sum(n * m for a, n in vol.items() for b, m in ang.items() if a + b <= -3),
    )


if __name__ == "__main__":
    import sys

    part = sys.argv[1] if len(sys.argv) > 1 else "check"
    if part == "check":
        main()
    elif part == "conductance":
        for (x, y), v in angular_weight().items():
            w = residue(7 * v, 49)
            if w:
                print(x, y, w)
    elif part == "resistance":
        for d in range(MAX + 1):
            for x in compositions(d, 5):
                a = tuple(sorted(x))
                for s in range(MAX - d + 1):
                    for y in compositions(s, 5):
                        w = residue(7**3 * weight(a, tuple(sorted(y))), 49)
                        if w:
                            print(*x, 1 - d, *y, 1 - s, w)
    else:
        raise ValueError(part)
