"""Genus-one bulk and marked contacts with one scalar descendant on the fivefold."""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
from random import Random

from one_descendant_hexagon import add, determinant, dot, transformed


def replaced(rows, at, wave):
    return determinant([wave if i == at else row for i, row in enumerate(rows)])


def contacts(scalars, vector, polarization):
    descendant = tuple(-x for x in add(*scalars, vector))
    div = dot(vector, polarization)
    marked, advected = transformed(scalars, vector, polarization, descendant)
    one_x = advected - div * determinant(scalars) ** 2
    two_x = -sum(replaced(scalars, i, vector) * replaced(scalars, i, descendant)
                 for i in range(5))
    bulk = 0
    for i, j in combinations(range(5), 2):
        pair = add(scalars[i], scalars[j])
        others = [p for l, p in enumerate(scalars) if l not in (i, j)]
        bulk += dot(pair, polarization) * determinant([pair, vector, *others]) ** 2
    terms = marked, Fraction(one_x, 4), -Fraction(bulk, 24), -Fraction(div * two_x, 12)
    assert sum(terms) == 0
    assert two_x == -sum(replaced(scalars, i, descendant) * replaced(scalars, i, vector)
                        for i in range(5))
    return terms


def polarized_primary_green():
    # Diagonal minus off-diagonal in the grounded unit-resistance-fraction simplex.
    g = lambda i, j: Fraction(min(i, j) * (6 - max(i, j)), 42)
    difference = sum(g(order[0], order[0]) - g(order[0], order[1])
                     for order in permutations(range(6)))
    return difference / (12 * factorial(5))


def higher_marks(scalars, v):
    k = tuple(-x for x in add(*scalars))
    d = dot(k, v)
    c = [dot(p, v) for p in scalars]
    subset4 = list(combinations(range(len(scalars)), 4))
    subset5 = list(combinations(range(len(scalars)), 5))
    h = {a: determinant([*(scalars[i] for i in a), k])**2 for a in subset4}
    bare = sum(h[a] * (2*d + sum(c[i] for i in a)) for a in subset4)
    scalar = sum(determinant([scalars[i] for i in b])**2 for b in subset5)
    primary = sum(c[i] * determinant([add(scalars[j], k) if i == j else scalars[j]
                                       for j in b])**2
                  for b in subset5 for i in range(len(scalars)))
    descendant = -d * (scalar + sum(h.values()))
    assert bare + primary + descendant == 0
    return (bare, primary, descendant)


def check():
    assert polarized_primary_green() == Fraction(1, 24)
    p = [tuple(int(i == j) for j in range(5)) for i in range(5)]
    example = contacts(p, p[0], p[1])
    assert example == (Fraction(1, 24), 0, -Fraction(1, 24), 0)
    print('Transverse probe (marked, one-x, bulk, two-x):', example)
    example = contacts(p, (1, 2, 3, 4, 5), p[1])
    assert all(example)
    print('Longitudinal probe:', example)
    rng = Random(1979)
    vec = lambda: tuple(rng.randint(-2, 2) for _ in range(5))
    for trial in range(40):
        rows = [vec() for i in range(5)]
        if trial % 3 == 0:
            rows[4] = add(rows[0], rows[1])
        contacts(rows, vec(), vec())
    for r in range(1, 5):
        values = higher_marks([vec() for _ in range(r + 4)], vec())
        assert values[0] and (r == 1 or values[1])
    print('Forty exact one-descendant tests and primary identities for F_1 through F_4 passed.')


if __name__ == '__main__':
    check()
