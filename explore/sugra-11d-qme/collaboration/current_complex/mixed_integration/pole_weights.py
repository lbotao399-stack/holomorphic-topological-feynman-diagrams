"""Integer interpolation of the simple 7-adic pole of the 166 conductance integral."""

from fractions import Fraction
from math import comb, factorial, prod
from collections import defaultdict


def value(i, j):
    return Fraction(
        prod(range(1, 2 * i, 2)) * prod(range(1, 2 * j, 2)),
        prod(range(13, 13 + 2 * (i + j), 2)),
    )


def weights():
    basis = [[1]]
    for k in range(11):
        prev = basis[-1]
        basis.append(
            [
                (prev[i - 1] if i else 0) - k * (prev[i] if i <= k else 0)
                for i in range(k + 2)
            ]
        )
    points = defaultdict(Fraction)
    for r in range(12):
        for s in range(12 - r):
            newton = sum(
                a * b * value(i, j)
                for i, a in enumerate(basis[r])
                for j, b in enumerate(basis[s])
            ) / (factorial(r) * factorial(s))
            for u in range(r + 1):
                for v in range(s + 1):
                    points[u, v] += (
                        newton * (-1) ** (r - u + s - v) * comb(r, u) * comb(s, v)
                    )
    total = {}
    for k, w in points.items():
        w *= 49
        assert w.denominator % 7
        total[k] = (w.numerator * pow(w.denominator, -1, 49)) % 49
    for i in range(12):
        for j in range(12 - i):
            w = 49 * value(i, j)
            want = w.numerator * pow(w.denominator, -1, 49) % 49
            assert sum(a * x**i * y**j for (x, y), a in total.items()) % 49 == want
    return {key: w for key, w in total.items() if w}


if __name__ == "__main__":
    for (x, y), w in sorted(weights().items()):
        print(x, y, w)
