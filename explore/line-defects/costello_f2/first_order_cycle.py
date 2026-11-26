"""Four-flavor cycle for the principal symbol in Costello Proposition 16.5.2."""

from collections import Counter
from fractions import Fraction
from itertools import product

UNIT, Z1, Z2 = (0, 0), (1, 0), (0, 1)


# A product of two Weyl symbols in this calculation has at most one linear
# argument. Only its unit and z1*z2 coefficients can survive the tests.
def star(left, right, order):
    if order == 0:
        return {tuple(x + y for x, y in zip(left, right)): 1}
    if (left, right) == (Z1, Z2):
        return {UNIT: Fraction(1, 2)}
    if (left, right) == (Z2, Z1):
        return {UNIT: Fraction(-1, 2)}
    return {}


def bracket(a, b, order):
    i, j, u = a
    k, l, v = b
    result = Counter()
    if j == k:
        for w, c in star(u, v, order).items():
            result[i, l, w] += c
    if i == l:
        for w, c in star(v, u, order).items():
            result[k, j, w] -= c
    return result


def pair(field, matrix):
    return field[:2] == matrix[::-1] and field[2] == UNIT


def test_outputs(v, w, labels):
    x, y = labels
    return int(pair(v, x) and pair(w, y)) + int(pair(v, y) and pair(w, x))


def one_cochain(field, selected, f, g):
    i, j, u = field
    if u != selected:
        return []
    # All single-trace maps are linear combinations of this form.
    return [((i, k, f), (k, j, g)) for k in range(4)]


A = ((0, 1, Z1), (1, 2, Z2), ((2, 3), (3, 0)))
B = ((0, 1, Z2), (1, 2, Z1), ((2, 3), (3, 0)))
C = ((0, 1, Z1), (2, 3, Z2), ((1, 2), (3, 0)))


def differential(test, selected, f, g, order):
    x, y, labels = test

    def act(a, b):
        value = 0
        for v, w in one_cochain(b, selected, f, g):
            value += sum(
                c * test_outputs(u, w, labels) for u, c in bracket(a, v, order).items()
            )
            value += sum(
                c * test_outputs(v, u, labels) for u, c in bracket(a, w, order).items()
            )
        return value

    return (
        act(x, y)
        - act(y, x)
        - sum(
            c
            * sum(test_outputs(v, w, labels) for v, w in one_cochain(u, selected, f, g))
            for u, c in bracket(x, y, order).items()
        )
    )


def first_order(test):
    x, y, labels = test
    i, j, u = x
    k, l, v = y
    if (u, v) == (Z2, Z1):
        return -first_order((y, x, labels))
    assert (u, v) == (Z1, Z2)
    value = test_outputs((i, l, UNIT), (k, j, UNIT), labels)
    if j == k:
        value -= sum(test_outputs((i, t, UNIT), (t, l, UNIT), labels) for t in range(4))
    if i == l:
        value -= sum(test_outputs((k, t, UNIT), (t, j, UNIT), labels) for t in range(4))
    return value


def main():
    survivors = []
    for order, T, f, g in product(
        range(2), (Z1, Z2, (1, 1)), (UNIT, Z1, Z2), (UNIT, Z1, Z2)
    ):
        a, b, c = (differential(test, T, f, g, order) for test in (A, B, C))
        assert a - b == c
        if c:
            survivors.append((order, T, f, g, c))
    assert len(survivors) == 4
    assert all(order == 1 and c == Fraction(1, 2) for order, _, _, _, c in survivors)
    values = tuple(map(first_order, (A, B, C)))
    assert values == (-1, 1, 1)
    print(
        "All normalized one-cochain boundaries vanish. First-order bracket",
        values,
        "pairing",
        values[0] - values[1] - values[2],
    )


if __name__ == "__main__":
    main()
