#!/usr/bin/env python3
"""Exact full-color certificate for the one-center two-loop cancellation.

The so(5) graph tensors are contracted from the six original cubic tensors,
without importing the repository's graph or cochain implementations.  All
24 input permutations and all six momentum coefficients are retained.
The mixed cubic bracket is checked in every central-input sector, including
the all-central sector.  Integer arithmetic is bounded before contraction.
"""

from itertools import combinations, permutations
from fractions import Fraction
import json

import numpy as np


def parity(order):
    return (-1) ** sum(order[i] > order[j]
                       for i in range(len(order))
                       for j in range(i + 1, len(order)))


PAIRS = tuple(combinations(range(5), 2))
DIM = len(PAIRS)
MONOMIALS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
# u=[p0,p1], v=[p0,p2], w=[p1,p2], p3=-p0-p1-p2.
AREAS = {
    (0, 1): (1, 0, 0), (0, 2): (0, 1, 0), (0, 3): (-1, -1, 0),
    (1, 2): (0, 0, 1), (1, 3): (1, 0, -1), (2, 3): (0, 1, 1),
}


def area(i, j):
    if i < j:
        return AREAS[i, j]
    return tuple(-x for x in AREAS[j, i])


def multiply_linear(a, b):
    return np.array([a[i] * b[j] + (a[j] * b[i] if i != j else 0)
                     for i, j in MONOMIALS], dtype=np.int64)


def permuted(tensor, order):
    # tensor(X_order[0],...,X_order[3]) as an array in the original slots.
    return tensor.transpose(tuple(order.index(i) for i in range(4)))


def alternate_area_pairs(tensor):
    answer = np.zeros((6,) + tensor.shape, dtype=np.int64)
    for order in permutations(range(4)):
        coefficient = multiply_linear(area(order[0], order[1]),
                                      area(order[2], order[3]))
        answer += parity(order) * coefficient[:, None, None, None, None] * permuted(tensor, order)
    return answer


def half_shuffle_square(seed):
    answer = np.zeros((6,) + seed.shape, dtype=np.int64)
    for first in combinations(range(4), 2):
        second = tuple(i for i in range(4) if i not in first)
        order = first + second
        coefficient = multiply_linear(area(first[0], first[1]),
                                      area(second[0], second[1]))
        answer += parity(order) * coefficient[:, None, None, None, None] * permuted(seed, order)
    assert not np.any(answer % 2)
    return answer // 2


def build_so5():
    basis = np.zeros((DIM, 5, 5), dtype=np.int64)
    for a, (i, j) in enumerate(PAIRS):
        basis[a, i, j], basis[a, j, i] = 1, -1
    structure = np.zeros((DIM, DIM, DIM), dtype=np.int64)
    for a in range(DIM):
        for b in range(DIM):
            bracket = basis[a] @ basis[b] - basis[b] @ basis[a]
            for c, (i, j) in enumerate(PAIRS):
                structure[a, b, c] = bracket[i, j]
            assert np.array_equal(np.einsum('c,cij->ij', structure[a, b], basis), bracket)
    assert np.max(np.abs(structure)) == 1
    assert np.array_equal(structure, -structure.swapaxes(0, 1))
    assert np.array_equal(structure, -structure.swapaxes(1, 2))
    killing = np.einsum('aij,bji->ab', structure, structure)
    assert np.array_equal(killing, -6 * np.eye(DIM, dtype=np.int64))
    # Jacobi in the complete adjoint tensor, not sampled external inputs.
    ff = np.einsum('abe,ecd->abcd', structure, structure)
    assert not np.any(ff + ff.transpose(1, 2, 0, 3) + ff.transpose(2, 0, 1, 3))
    return structure


def graph_certificate(structure):
    # At metric delta, every original graph tensor entry is a sum of at
    # most DIM**7 signed unit contributions.  This and subsequent sums fit.
    integer_bound = 24 * 2 * DIM**7
    assert integer_bound < 2**63
    graph_133 = -np.einsum('xyz,ayu,buw,czs,dsr,xwr->abcd',
                          *([structure] * 6), optimize=True)
    graph_223 = -np.einsum('xyz,axu,byw,czs,dsr,uwr->abcd',
                          *([structure] * 6), optimize=True)
    F133 = alternate_area_pairs(graph_133)
    F223 = alternate_area_pairs(graph_223)
    identity = np.eye(DIM, dtype=np.int64)
    Gkk = alternate_area_pairs(np.einsum('ab,cd->abcd', identity, identity))
    ff = np.einsum('abe,cde->abcd', structure, structure)
    beta = np.zeros_like(F223)
    for i, j in combinations(range(4), 2):
        k, l = tuple(v for v in range(4) if v not in (i, j))
        coefficient = multiply_linear(area(k, l), area(k, l))
        beta += ((-1) ** (i+j)) * coefficient[:, None, None, None, None] * permuted(ff, (i, j, k, l))
    # Killing=-6 delta: F_K=-F_delta/6, beta_K=-6 beta_delta,
    # G_K=36 G_delta.  These are the denominator-cleared identities.
    assert not np.any(F133 + F223 - 18 * beta)
    assert not np.any(F223 + 36 * beta - 12 * Gkk)
    # The two invariant cochains are independent over Q.
    flat_b, flat_g = beta.ravel(), Gkk.ravel()
    first = next(i for i, v in enumerate(flat_b) if v)
    second = next(i for i in range(len(flat_b))
                  if int(flat_b[first]) * int(flat_g[i]) != int(flat_b[i]) * int(flat_g[first]))
    minor = int(flat_b[first]) * int(flat_g[second]) - int(flat_b[second]) * int(flat_g[first])
    return {"graph_color_momentum_entries": int(F223.size),
            "full_graph_identity": "F223 = -beta - Gkk/18 (Killing)",
            "sum_identity": "F133 + F223 = beta/2 (Killing)",
            "independence_minor_in_delta_units": minor,
            "integer_bound": integer_bound}, beta, F133, F223


def central_certificate(structure):
    count, central = DIM + 1, DIM
    mixed = np.zeros((count, count, count), dtype=np.int64)
    pure = np.zeros_like(mixed)
    for a in range(DIM):
        for order in set(permutations((a, a, central))):
            mixed[order] = -6  # Killing on the simple factor.
    pure[central, central, central] = 1
    extended_structure = np.zeros_like(mixed)
    extended_structure[:DIM, :DIM, :DIM] = structure
    # Verify the complete classical differential, before forming the square.
    for cubic in (mixed, pure):
        differential = np.zeros((3, count, count, count, count), dtype=np.int64)
        insertion = np.einsum('abe,ecd->abcd', extended_structure, cubic)
        for i, j in combinations(range(4), 2):
            k, l = tuple(v for v in range(4) if v not in (i, j))
            differential += ((-1) ** (i+j)) * np.array(area(k, l))[:, None, None, None, None] * permuted(insertion, (i, j, k, l))
        assert not np.any(differential)
    # Four independently stored coefficients: chi²/ell, chi²,
    # chi*gamma/ell, gamma²/ell.  Internal simple inverse metric=-delta/6.
    seed0 = np.einsum('ab,cd->abcd', mixed[:, :, central], mixed[central])
    seed1_numerator = -np.einsum('abe,ecd->abcd', mixed[:, :, :DIM], mixed[:DIM])
    assert not np.any(seed1_numerator % 6)
    seed1 = seed1_numerator // 6
    seed2 = (np.einsum('ab,cd->abcd', mixed[:, :, central], pure[central])
             + np.einsum('ab,cd->abcd', pure[:, :, central], mixed[central]))
    seed3 = np.einsum('ab,cd->abcd', pure[:, :, central], pure[central])
    squares = tuple(half_shuffle_square(seed) for seed in (seed0, seed1, seed2, seed3))
    # Independent direct alternation formula gives precisely the same 1/8.
    for seed, square in zip((seed0, seed1, seed2, seed3), squares):
        assert np.array_equal(alternate_area_pairs(seed), 8 * square)
    sectors = {}
    for c in range(5):
        colors = [index for index in np.ndindex((count,) * 4)
                  if index.count(central) == c]
        sectors[c] = len(colors)
        if c in (1, 3, 4):
            for index in colors:
                assert all(not np.any(square[(slice(None),) + index]) for square in squares)
    simple_g = alternate_area_pairs(np.einsum('ab,cd->abcd',
                                             -6*np.eye(DIM, dtype=np.int64),
                                             -6*np.eye(DIM, dtype=np.int64)))
    assert np.array_equal(8*squares[0][:, :DIM, :DIM, :DIM, :DIM], simple_g)
    for square in squares[1:]:
        assert not np.any(square[:, :DIM, :DIM, :DIM, :DIM])
    area_product = multiply_linear(area(0, 1), area(2, 3))
    for a in range(DIM):
        for b in range(DIM):
            expected = (-6 if a == b else 0) * area_product
            assert not np.any(squares[0][:, a, b, central, central])
            assert np.array_equal(squares[1][:, a, b, central, central], -expected)
            assert np.array_equal(squares[2][:, a, b, central, central], expected)
            assert not np.any(squares[3][:, a, b, central, central])
    # At gamma=ell*chi, all mixed sectors cancel, including all permutations.
    assert not np.any(squares[1] + squares[2])
    return {"central_sector_color_counts": sectors,
            "full_bracket_color_momentum_entries": int(sum(x.size for x in squares)),
            "central_condition": "gamma = ell*chi", "all_sectors_verified": True}, squares


def main():
    structure = build_so5()
    graphs, beta, F133, F223 = graph_certificate(structure)
    central, squares = central_certificate(structure)
    # Full order-hbar² QME at ell=epsilon*U=1 and gamma=chi.
    # Multiply q+epsilon*o2-d*i2 by 276480*pi^4.  Extend the
    # original graphs and primitive by zero on central external fields.
    residual = -40 * sum(squares)
    residual[:, :DIM, :DIM, :DIM, :DIM] += 3*(-3*F133 + 2*F223) + 702*beta
    assert not np.any(residual)
    central["full_qme_entries_verified"] = int(residual.size)
    assert Fraction(-1, 18) / 384 == Fraction(-1, 6912)
    assert -(Fraction(-1, 48) - Fraction(1, 160)) / 64 == Fraction(13, 30720)
    print(json.dumps({"arithmetic": "exact signed integers and fractions", "so5": graphs,
                      "central_extension": central,
                      "chi_squared_pi4_at_ell_epsilon_U_one": "-1/6912",
                      "i2_coefficient_pi4_at_epsilon_U_one": "13/30720"}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
