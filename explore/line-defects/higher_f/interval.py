"""Extract a single ordered coefficient by iterated polynomial antiderivatives."""

from collections import defaultdict
from fractions import Fraction as F
from math import comb
from coefficients import ONE, add, positive, times


def ordered(r, t, outputs):
    n = len(outputs) - 1
    if n < 1:
        raise ValueError("at least two outputs required")
    if (
        not positive(r, t, *outputs)
        or add(*outputs, times(n, ONE)) != add(r, t)
        or min(sum(r), sum(t)) < n
    ):
        return F()
    # Exponents are those of X1, X2 and the upper limit of integration.
    polynomial = {(0, 0, 0): F(1)}
    for k in reversed(outputs):
        result = defaultdict(F)
        weight = (sum(k) + 1) * comb(sum(k), k[0])
        for h in range(min(k[0], r[0]) + 1):
            for j in range(min(k[1], r[1]) + 1):
                remaining = sum(k) - h - j
                c = weight * comb(k[0], h) * comb(k[1], j)
                for i in range(remaining + 1):
                    v = (-1) ** i * comb(remaining, i) * c
                    for (a, b, degree), old in polynomial.items():
                        if a + h <= r[0] and b + j <= r[1]:
                            m = degree + h + j + i + 1
                            result[a + h, b + j, m] += v * old / m
        polynomial = {k: v for k, v in result.items() if v}
    value = F()
    for (a, b, _), v in polynomial.items():
        j = r[1] - b
        if 0 <= j <= n and r[0] - a == n - j:
            value -= (-1) ** j * comb(n, j) * v
    return value / (comb(sum(r), r[0]) * comb(sum(t), t[0]))
