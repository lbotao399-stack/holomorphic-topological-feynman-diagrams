"""Exact balanced planar symbol of the four-current in the mixed 4+5 class."""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import permutations
from math import factorial, prod
from mixed_current_lift import conditions, lift
from mixed_planar import normalize


@cache
def allocations(n):
    return [(i, j, n - i - j) for i in range(n + 1) for j in range(n - i + 1)]


def balance_rows(words, residual=(0, 0)):
    # Evaluate the marked momentum at minus the three ordinary momenta,
    # plus the residual. Take its indicated Taylor coefficient.
    rows = defaultdict(Counter)
    for kind, w in words:
        exponents = [(i, h + 2 - i) for _, h, i in w]
        mx, my = exponents.pop()
        rx, ry = residual
        if mx < rx or my < ry:
            continue
        degree = mx + my - rx - ry
        denom = (
            factorial(rx)
            * factorial(ry)
            * prod(factorial(i) * factorial(j) for i, j in exponents)
        )
        for dx in allocations(mx - rx):
            for dy in allocations(my - ry):
                dest = [(x + i, y + j) for (x, y), i, j in zip(exponents, dx, dy)]
                word, sign = normalize(dest)
                if sign:
                    denominator = denom * prod(map(factorial, dx + dy))
                    rows[kind, word][kind, w] += Fraction(
                        (-1) ** degree * sign, denominator
                    )
    return {k: dict((i, c) for i, c in row.items() if c) for k, row in rows.items()}


def symbol(current, waves, residual):
    # Independent direct alternating evaluation of the differentiated exponentials.
    values = Counter()
    for (kind, w), coefficient in current.items():
        for perm in permutations(range(3)):
            _, sign = normalize(perm)
            v = Fraction(sign * coefficient)
            for j, (_, h, i) in enumerate(w):
                a, b = i, h + 2 - i
                if j == 3:
                    a -= residual[0]
                    b -= residual[1]
                if min(a, b) < 0:
                    v = 0
                    break
                x, y = waves[perm[j] if j < 3 else 3]
                v *= Fraction(x**a * y**b, factorial(a) * factorial(b))
            values[kind] += v
    return tuple(values[k] for k in ("dx", "dy"))


def verify():
    def vanishing(words):
        return balance_rows(words).values()

    current, scalar, pairing, ncoords, rank, nrows = lift(vanishing)
    assert {w[-1][1] for _, w in current} == {0, 1, 2}
    assert {w[-1][1] for w in scalar} == {0, 1, 2}
    print(
        "balanced planar four-current",
        ncoords,
        rank,
        nrows,
        len(current),
        len(scalar),
        pairing,
    )
    for residual in ((0, 0), (0, 1), (1, 0)):
        rows = balance_rows(current, residual)
        coefficients = {
            k: sum(a * current[w] for w, a in row.items()) for k, row in rows.items()
        }
        coefficients = {k: v for k, v in coefficients.items() if v}
        if residual == (0, 0):
            assert not coefficients
        else:
            assert coefficients
        for pts in [((1, 2), (2, 1), (3, -2)), ((-2, 3), (1, -3), (4, 1))]:
            waves = (*pts, tuple(-sum(p[i] for p in pts) for i in range(2)))
            check = Counter()
            for (kind, w), v in coefficients.items():
                for perm in permutations(range(3)):
                    _, sign = normalize(perm)
                    check[kind] += (
                        v
                        * sign
                        * prod(
                            pts[j][0] ** a * pts[j][1] ** b
                            for (a, b), j in zip(w, perm)
                        )
                    )
            assert tuple(check[k] for k in ("dx", "dy")) == symbol(
                current, waves, residual
            )
        print(
            "residual", residual, "independent nonzero coefficients", len(coefficients)
        )

    fixture = ((1, 2), (2, 1), (3, -2), (-6, -1))
    assert symbol(current, fixture, (1, 0))[0] == Fraction(70312, 9)

    def first_jets(words):
        for residual in ((0, 0), (0, 1), (1, 0)):
            yield from balance_rows(words, residual).values()

    _, _, first, target, remainder = conditions(first_jets)
    assert len(first) == 222 and not remainder(target)
    assert (ncoords, rank, nrows, len(current), len(scalar), pairing) == (
        224,
        219,
        375,
        116,
        142,
        73920,
    )
    print(
        "balanced values and slopes have joint rank",
        len(first),
        "and force zero pairing",
    )
    return current, scalar, pairing


if __name__ == "__main__":
    verify()
