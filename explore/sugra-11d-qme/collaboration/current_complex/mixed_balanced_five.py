"""Test all balanced transverse constraints on the mixed five-vector class."""

from collections import defaultdict
from functools import reduce
from itertools import permutations
from math import gcd
from sys import stdin
from mixed_divergence_test import relation
from mixed_probe_jets import matrix, pivots
from mixed_vertical import slots


def verify(lines):
    data = defaultdict(list)
    for line in lines:
        label, *row = map(int, line.split())
        data[label].append(row)
    closure, balanced, (target,), unmarked, first = [data[k] for k in range(5)]
    assert len(closure) == len(balanced) == len(unmarked) == len(first) == 169

    # Every marked column contains positive powers of the function momentum.
    # On setting f momentum to minus the sum of the five vector momenta,
    # the resulting function-independent column therefore vanishes on a
    # further restriction to five balanced vector momenta.
    rho = (2, 1, 0)
    multiplicity = 0
    for perm in permutations(rho):
        weight = tuple(x + y - z for x, y, z in zip((-2, -3, -3), rho, perm))
        sign = (-1) ** sum(i < j for n, i in enumerate(perm) for j in perm[n + 1 :])
        multiplicity += sign * sum(1 for _ in slots(8, 5, weight))
    assert multiplicity == 75
    assert len(pivots(unmarked)) == 16 and len(pivots(balanced)) == multiplicity - 16

    # The previous filtered CE calculation provides 88 boundaries and one
    # additional closed class in the 138-dimensional source.
    assert len(pivots(closure)) == 138 - 88 - 1
    equations = closure + balanced
    cols = pivots(equations + [target])
    rows = pivots(zip(*([r[j] for j in cols] for r in equations)))
    assert (len(cols), len(rows)) == (107, 106)
    square = matrix(
        [[equations[i][j] for j in cols] for i in rows] + [[target[j] for j in cols]]
    )
    values, scale = square.solve_den(matrix([[0]] * len(rows) + [[1]]))
    values = list(map(int, values.to_Matrix()))
    common = reduce(gcd, values, int(scale))
    values = [v // common for v in values]
    scale = int(scale) // common
    assert all(sum(r[j] * v for j, v in zip(cols, values)) == 0 for r in equations)
    assert sum(target[j] * v for j, v in zip(cols, values)) == scale and scale
    print(
        "unmarked multiplicity, rank on five balanced momenta, marked balanced rank",
        75,
        16,
        59,
    )
    print(
        "closed balanced witness support",
        sum(bool(v) for v in values),
        "detector",
        scale,
    )
    relation(equations + first, target, 132)


if __name__ == "__main__":
    verify(stdin)
