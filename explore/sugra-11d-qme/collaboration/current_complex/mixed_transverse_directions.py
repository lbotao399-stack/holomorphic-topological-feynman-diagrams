"""Five-vector currents and their coboundaries on three transverse directions."""

from collections import defaultdict
from functools import reduce
from math import gcd
from sys import stdin
import numpy as np
from mixed_probe_jets import matrix, pivots


def solve_modulo(square, rhs):
    n = len(square)
    system = np.column_stack((square, rhs))
    for col in range(n):
        pivot = col + np.flatnonzero(system[col:, col])[0]
        system[[pivot, col]] = system[[col, pivot]]
        system[col] = system[col] * pow(int(system[col, col]), -1, 251) % 251
        factors = system[:, col].copy()
        factors[col] = 0
        system = (system - factors[:, None] * system[col]) % 251
    return system[:, n:]


def four_lines(rows, closure):
    eq = np.array([[x % 251 for x in row] for row in rows], dtype=np.int64)
    incoming, currents = eq[:, :132], eq[:, 132:]
    cols = pivots(incoming)
    basis = pivots(incoming[:, cols].T)
    assert len(cols) == len(basis) == 86
    inv_times_rows = solve_modulo(incoming[basis][:, cols], currents[basis])
    residues = (currents - incoming[:, cols] @ inv_times_rows) % 251
    closed_rows = np.array([[x % 251 for x in row] for row in closure], dtype=np.int64)
    closure_cols = pivots(closed_rows)
    basis_closure = pivots(closed_rows[:, closure_cols].T)
    assert len(closure_cols) == 49
    relations = solve_modulo(
        closed_rows[basis_closure][:, closure_cols], closed_rows[basis_closure]
    )
    quotient = (residues - residues[:, closure_cols] @ relations) % 251
    witness = np.flatnonzero(np.any(quotient, axis=1))[0]

    square = matrix([[rows[i][j] for i in basis] for j in cols])
    numer, denom = square.solve_den(matrix([[rows[witness][j]] for j in cols]))
    coeff = list(map(int, numer.to_Matrix()))
    common = reduce(gcd, coeff, int(denom))
    coeff = [v // common for v in coeff]
    denom = int(denom) // common
    functional = [
        denom * x - sum(c * rows[i][j] for c, i in zip(coeff, basis))
        for j, x in enumerate(rows[witness])
    ]
    assert functional[:132] == [0] * 132
    assert len(pivots(closure + [functional[132:]])) == 50
    print(
        "four wave lines: exact boundary cancellations",
        132,
        "closure and augmented ranks",
        49,
        50,
    )
    print(
        "number of wave evaluations in the four-line detector",
        1 + sum(bool(x) for x in coeff),
    )


def verify(lines):
    data = defaultdict(list)
    for line in lines:
        label, *row = map(int, line.split())
        data[label].append(row)
    # All five independent wave scales and the three marked momentum
    # components have already been expanded as formal variables.
    assert not data[2]
    restriction = matrix(data[1])
    assert restriction.shape[1] == 138 and restriction.rank() == 18
    inputs = np.array([[x % 251 for x in r] for r in data[0]], dtype=np.int64)
    evaluations, boundaries = inputs[:, :138], inputs[:, 138:]
    basis_rows = pivots(evaluations.T)
    assert len(basis_rows) == 138
    # The selected integral evaluation minor is a unit at 251. This lets us
    # compute the rational differential in these coordinates modulo 251.
    differential = solve_modulo(evaluations[basis_rows], boundaries[basis_rows])
    assert np.array_equal(evaluations @ differential % 251, boundaries)
    restricted_boundaries = np.array(data[1], dtype=np.int64) @ differential % 251
    assert len(pivots(restricted_boundaries)) == 18
    print("three wave lines 2+2+1: current dimension and boundary rank", 18, 18)
    print("three wave lines 3+1+1: all 138 current symbols vanish")
    four_lines(data[3], data[4])


if __name__ == "__main__":
    verify(stdin)
