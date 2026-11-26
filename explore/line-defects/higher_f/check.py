"""Recurrence, Schwinger graph, ordered interval and special-family cross-checks."""

from fractions import Fraction as F
from itertools import product
from math import prod
from random import Random
from coefficients import coefficient, add, sub
from fan import graph, transferred
from interval import ordered


def integer_lists(total, length):
    if length == 1:
        yield (total,)
        return
    for i in range(total + 1):
        for rest in integer_lists(total - i, length - 1):
            yield (i, *rest)


def check(n, r, t, k, compact=False):
    expected = coefficient(r, t, k)
    assert (
        expected
        == coefficient(r, t, k, True)
        == ordered(r, t, k)
        == transferred(r, t, k)
    ), (n, r, t, k)
    if not compact:
        assert expected == graph(r, t, k), (n, r, t, k)
    return expected


def main():
    count = 0
    for n in (1, 2):
        for extra in range(4):
            for a in range(extra + 1):
                b = extra - a
                for kx in integer_lists(a, n + 1):
                    for ky in integer_lists(b, n + 1):
                        k = tuple(zip(kx, ky))
                        total = (a + n, b + n)
                        for i, j in product(range(total[0] + 1), range(total[1] + 1)):
                            r = (i, j)
                            t = sub(total, r)
                            check(n, r, t, k)
                            count += 1
    rng = Random(2602)
    for n in range(1, 10):
        for _ in range(80):
            k = [[0, 0] for j in range(n + 1)]
            for j in range(rng.randrange(8)):
                k[rng.randrange(n + 1)][rng.randrange(2)] += 1
            k = tuple(map(tuple, k))
            total = add(*k, (n, n))
            r = tuple(rng.randrange(j + 1) for j in total)
            check(n, r, sub(total, r), k, True)
        a = tuple(rng.randrange(4) for j in range(n + 1))
        expected = (-1) ** (n + 1) * prod(
            F(a[j] + 1, sum(a[j:]) + len(a) - j) for j in range(n)
        )
        assert (
            coefficient((0, n + sum(a)), (n, 0), tuple((0, i) for i in a)) == expected
        )
    print(
        count,
        "exhaustive labeled one/two-loop probes and 720 further probes through nine loops passed",
    )


if __name__ == "__main__":
    main()
