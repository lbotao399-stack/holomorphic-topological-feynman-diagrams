"""Prove the mixed detector annihilates every strictly conserved transverse cocycle."""

from collections import defaultdict
from functools import reduce
from math import gcd
from sys import stdin
from itertools import permutations
from mixed_vertical import slots

from mixed_probe_jets import matrix, pivots


def verify(lines):
    data = defaultdict(list)
    for line in lines:
        k, *v = map(int, line.split())
        data[k].append(v)
    (degrees,) = data[0]
    (detector,) = data[3]
    assert [degrees.count(j) for j in (1, 2, 3)] == [84, 45, 9]
    assert len(data[1]) == len(data[2]) == 172
    assert len(data[4]) == 172
    rho = (2, 1, 0)
    multiplicities = []
    for marked_degree in range(1, 6):
        hw = (marked_degree - 3, -3, -3)
        n = 0
        for perm in permutations(rho):
            weight = tuple(t + x - y for t, x, y in zip(hw, rho, perm))
            sign = (-1) ** sum(
                perm[i] < perm[j] for i in range(3) for j in range(i + 1, 3)
            )
            n += sign * sum(1 for _ in slots(9 - marked_degree, 4, weight))
        multiplicities.append(n)
    assert multiplicities == [31, 36, 22, 10, 1]
    dimensions = []
    for q in range(2, 9):
        dimension = 0
        for marked_degree in range(1, 10 - q):
            hw = (marked_degree - 3, -3, -3)
            for perm in permutations(rho):
                weight = tuple(t + x - y for t, x, y in zip(hw, rho, perm))
                sign = (-1) ** sum(
                    perm[i] < perm[j] for i in range(3) for j in range(i + 1, 3)
                )
                dimension += sign * sum(1 for _ in slots(9 - marked_degree, q, weight))
        dimensions.append(dimension)
    assert dimensions == [3, 30, 100, 150, 105, 31, 3]
    scalar_size = sum(multiplicities)
    print("scalar primitive multiplicities", multiplicities)
    relation(data[1] + data[2], detector, 110)
    equations = [r + [0] * scalar_size for r in data[1]]
    equations.extend(r + u for r, u in zip(data[2], data[4]))
    whole = detector + [0] * scalar_size
    relation(equations, whole, 123)
    assert len(pivots(data[4])) == 73 and len(pivots(data[5])) == 76
    print("scalar five-cohomology dimension", 150 - 73 - 76)
    if data[11]:
        closure = [r + [0] * scalar_size for r in data[1]]
        for second in (12, 13):
            relation(closure + data[11] + data[second], whole, 123)
        witness(closure + data[11] + data[14], whole, 120)


def witness(equations, detector, expected):
    cols = pivots(equations + [detector])
    rows = pivots(zip(*([r[j] for j in cols] for r in equations)))
    assert len(rows) == expected and len(cols) == expected + 1
    square = matrix(
        [[equations[i][j] for j in cols] for i in rows] + [[detector[j] for j in cols]]
    )
    rhs = matrix([[0] for _ in rows] + [[1]])
    num, den = square.solve_den(rhs)
    values = list(map(int, num.to_Matrix()))
    assert all(sum(r[j] * v for j, v in zip(cols, values)) == 0 for r in equations)
    assert sum(detector[j] * v for j, v in zip(cols, values)) == den and den
    print("linear and quartic scalar response can vanish on a detected closed current")


def relation(equations, detector, expected):
    # Modular pivots choose a rationally invertible square submatrix.
    # Integer equality on every column proves the constraint on every cochain.
    cols = pivots(equations)
    rows = pivots(zip(*([r[j] for j in cols] for r in equations)))
    assert len(cols) == len(rows) == expected
    square = matrix([[equations[i][j] for i in rows] for j in cols])
    target = matrix([[detector[j]] for j in cols])
    numer, den = square.solve_den(target)
    coefficients = list(map(int, numer.to_Matrix()))
    scale = reduce(gcd, coefficients, int(den))
    den = int(den) // scale
    coefficients = [v // scale for v in coefficients]
    for j, x in enumerate(detector):
        assert sum(c * equations[i][j] for c, i in zip(coefficients, rows)) == den * x
    print(
        "exact constraints",
        len(rows),
        "tested cochain coordinates",
        len(detector),
        "nonzero coefficients",
        sum(v != 0 for v in coefficients),
    )


if __name__ == "__main__":
    verify(stdin)
