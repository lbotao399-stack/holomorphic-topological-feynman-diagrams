"""Bound the four-plus-five class with an integral five-vector transverse cycle."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod
from pathlib import Path
import sys

from complete_nodes import partitions
from ordered_shapes import determinant
from single_planar_nodes import fields
from theta_denominators import SHAPES
from theta_spectral_bound import wick_bound


ROOT = Path(__file__).resolve().parent
RATIO = Fraction(
    -1830654199083692458460335177513758417503026861512874941762533918423721792,
    3557,
)


def geometry(choice, mask, eta):
    momenta, polar = fields(0, 0)
    directions = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, -1, -1), (-1, -1, 1)]
    if choice == 2:
        directions[-1] = (1, 1, -1)
    function = ((1, 0, 0), (1, 0, -1), (0, 1, -1))[choice]
    momenta[3] = (*momenta[3][:2], *function)
    for j, p in enumerate(directions):
        pivot = next(i for i, x in enumerate(p) if x)
        axis = (pivot + 1 + ((mask >> j) & 1)) % 3
        u = [0] * 3
        u[axis], u[pivot] = p[pivot], -p[axis]
        momenta[4 + j] = (0, 0, *p)
        polar[4 + j] = (0, 0, *u)
    momenta[9] = tuple(-sum(col) for col in zip(*momenta[:9]))
    polar[9] = tuple(int(i == 2 + eta) for i in range(5))
    return momenta, polar


def physical_bound():
    rows = defaultdict(list)
    seen = set()
    for line in (ROOT / "four_small_weights.tsv").read_text().splitlines():
        w, choice, mask, eta = map(int, line.split())
        key = choice, mask, eta
        assert key not in seen and 0 <= mask < 32 and 0 <= eta < 3
        seen.add(key)
        _, polar = geometry(*key)
        rows[choice].append((w, polar))
    assert len(seen) == 201
    assert sum(abs(w) for group in rows.values() for w, _ in group) == 4460
    total = Fraction(0)
    for choice, samples in rows.items():
        momenta, _ = geometry(choice, 0, 0)
        kernel = wick_bound(momenta, samples)
        for shape, raw, _ in SHAPES:
            det = sum(abs(determinant(momenta, shape, order)) for order in partitions(shape))
            orders = prod(factorial(n - 1) for n in shape)
            total += kernel * Fraction(665280 * det * orders, 2 * raw)
    return abs(RATIO) * total


def verify(lines):
    from sympy import ZZ
    from sympy.polys.matrices import DomainMatrix

    samples = defaultdict(list)
    for line in lines:
        tag, *values = map(int, line.split())
        samples[tag].append(values)
    assert [len(samples[k]) for k in range(3)] == [100, 1, 1]
    old, small = samples[1][0], samples[2][0]
    assert old[:132] == small[:132] == [0] * 132
    closed = DomainMatrix.from_list(samples[0], ZZ).nullspace().to_Matrix()
    assert closed.shape == (89, 138)
    nonzero = False
    for v in closed.tolist():
        a = sum(x * int(y) for x, y in zip(old[132:], v))
        b = sum(x * int(y) for x, y in zip(small[132:], v))
        assert RATIO.denominator * a == RATIO.numerator * b
        nonzero |= bool(b)
    assert nonzero
    print("201 small polarizations, 132 boundary cancellations, 89 closed cochain ratios")


if __name__ == "__main__":
    verify(sys.stdin)
