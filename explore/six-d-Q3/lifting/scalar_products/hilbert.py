"""Graded restricted scalar dimensions, with an independent permutation average."""

from fractions import Fraction
from itertools import permutations
from math import factorial


def table(rank, degree):
    H = [[[0] * (degree + 1) for _ in range(degree + 1)] for _ in range(rank + 1)]
    H[0][0][0] = 1
    for a in range(degree + 1):
        for b in range(degree + 1):
            for n in range(1, rank + 1):
                for i in range(a, degree + 1):
                    for j in range(b, degree + 1):
                        H[n][i][j] += H[n - 1][i - a][j - b]
    raw = H[rank]
    return [
        [
            raw[i][j]
            - (raw[i - 1][j] if i else 0)
            - (raw[i][j - 1] if j else 0)
            + (raw[i - 1][j - 1] if i and j else 0)
            for j in range(degree + 1)
        ]
        for i in range(degree + 1)
    ]


def permutation_table(rank, degree):
    values = [[0] * (degree + 1) for _ in range(degree + 1)]
    for perm in permutations(range(rank)):
        remaining = set(range(rank))
        series = [1] + [0] * degree
        while remaining:
            i = min(remaining)
            cycle = 0
            while i in remaining:
                remaining.remove(i)
                cycle += 1
                i = perm[i]
            for j in range(cycle, degree + 1):
                series[j] += series[j - cycle]
        traceless = [series[k] - (series[k - 1] if k else 0) for k in range(degree + 1)]
        for i in range(degree + 1):
            for j in range(degree + 1):
                values[i][j] += traceless[i] * traceless[j]
    return [[Fraction(n, factorial(rank)) for n in row] for row in values]


if __name__ == "__main__":
    for n in range(2, 8):
        values = table(n, 7)
        assert values == permutation_table(n, 7)
        print(n, [values[i][i] for i in range(1, 7)])
