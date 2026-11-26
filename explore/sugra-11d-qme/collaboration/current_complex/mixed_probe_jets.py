"""Check the two detector normalizations and which marked jets measure their class."""

from collections import defaultdict
from fractions import Fraction
from functools import reduce
from math import gcd
from sys import stdin

import numpy as np
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix


def pivots(rows):
    matrix = np.array([[int(x) % 251 for x in row] for row in rows], dtype=np.int64)
    result = []
    for j in range(matrix.shape[1]):
        r = len(result)
        nz = np.flatnonzero(matrix[r:, j])
        if len(nz):
            i = r + nz[0]
            matrix[[i, r]] = matrix[[r, i]]
            matrix[r] = matrix[r] * pow(int(matrix[r, j]), -1, 251) % 251
            matrix[r + 1 :] = (
                matrix[r + 1 :] - matrix[r + 1 :, j, None] * matrix[r]
            ) % 251
            result.append(j)
    return result


def matrix(rows):
    return DomainMatrix.from_list(rows, ZZ)


def verify(lines):
    data = defaultdict(list)
    for line in lines:
        kind, *row = map(int, line.split())
        data[kind].append(row)
    assert [len(data[i]) for i in range(6)] == [65, 110, 110, 110, 89, 94]
    outgoing = matrix(data[0])
    kernel = outgoing.nullspace()
    assert kernel.shape == (47, 84)
    closed = list(map(int, kernel.to_Matrix().row(0)))
    g = reduce(gcd, closed)
    closed = [x // g for x in closed]
    assert all(sum(a * b for a, b in zip(row, closed)) == 0 for row in data[0])
    blocks = []
    for s in (1, 2, 3):
        rows = data[s]
        independent = pivots(zip(*(r[216:] for r in rows)))
        assert len(independent) == (84, 45, 9)[s - 1]
        blocks.append(
            [
                r[:132] + [sum(a * b for a, b in zip(r[132:216], closed))]
                for r in (rows[i] for i in independent)
            ]
        )
    first_two = pivots(blocks[0] + blocks[1])
    assert len(first_two) == 89 and first_two[-1] == 132
    first_last = blocks[0] + blocks[2]
    columns = pivots(first_last)
    assert len(columns) == 87 and 132 not in columns
    independent = pivots(zip(*([r[j] for j in columns] for r in first_last)))
    square = matrix([[first_last[i][j] for j in columns] for i in independent])
    rhs = matrix([[first_last[i][-1]] for i in independent])
    primitive, den = square.solve_den(rhs)
    v = list(map(int, primitive.to_Matrix()))
    assert all(
        sum(r[j] * c for j, c in zip(columns, v)) == r[-1] * den for r in first_last
    )
    pairings = []
    for kind in (4, 5):
        totals = [sum(r[0] * r[j] for r in data[kind]) for j in range(1, 217)]
        assert totals[:132] == [0] * 132
        pairings.append(sum(i * j for i, j in zip(totals[132:], closed)))
    assert all(pairings)
    print(
        "Linear and quadratic ranks 88 and 89; linear and cubic primitive equality verified"
    )
    print("132 exact boundary pairings vanish for each detector")
    print("Old divided by single-momentum class normalization", Fraction(*pairings))


if __name__ == "__main__":
    verify(stdin)
