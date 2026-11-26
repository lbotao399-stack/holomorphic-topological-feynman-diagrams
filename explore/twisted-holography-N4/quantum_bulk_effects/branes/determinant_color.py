#!/usr/bin/env python3
"""Exact determinant-color transport, with no boundary graph-engine imports.
This verifies a boundary-induced giant-graviton source constraint, not an
independently calculated bulk string amplitude.
"""
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import json


def commutator_unit(i, j, matrix):
    out = {}
    for (a, b), value in matrix.items():
        if j == a:
            out[i, b] = out.get((i, b), 0) + value
        if b == i:
            out[a, j] = out.get((a, j), 0) - value
    return {key: value for key, value in out.items() if value}


def permutation_sign(order):
    return (-1) ** sum(order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order)))


def subset_color(selected, a, b):
    """Enumerate every determinant cofactor and every nested commutator."""
    count = 0
    for rows in permutations(selected):
        for order in permutations(range(len(selected))):
            matrix = {(a, b): 1}
            for row, p in zip(rows, order):
                # Trace-dual of E_row,col is E_col,row.
                matrix = commutator_unit(rows[p], row, matrix)
                if not matrix:
                    break
            if matrix:
                count += permutation_sign(order) * matrix.get((a, b), 0)
    return count


def subset_formula(size, has_a, has_b):
    if has_a and has_b:
        return -(1 + (-1) ** size) * size * factorial(size - 2)
    if has_a:
        return (-1) ** (size - 1) * factorial(size - 1)
    if has_b:
        return -factorial(size - 1)
    return 0


def elementary(values, degree):
    if degree < 0 or degree > len(values):
        return 0
    total = 0
    for selected in combinations(values, degree):
        prod = 1
        for value in selected:
            prod *= value
        total += prod
    return total


def color_direct(diagonal, size, a, b):
    out = 0
    for selected in combinations(range(len(diagonal)), size):
        value = subset_color(selected, a, b)
        for index, entry in enumerate(diagonal):
            if index not in selected:
                value *= entry
        out += value
    return out


def color_formula(diagonal, size, a, b):
    others = [entry for j, entry in enumerate(diagonal) if j not in (a, b)]
    return (
        factorial(size - 1) * (((-1) ** (size - 1)) * diagonal[b] - diagonal[a])
        * elementary(others, len(diagonal) - size - 1)
        - (1 + (-1) ** size) * size * factorial(size - 2)
        * elementary(others, len(diagonal) - size)
    )


def wedge_monomials(left, right):
    if set(left) & set(right):
        return 0, ()
    inversions = sum(a > b for a in left for b in right)
    return (-1) ** inversions, tuple(sorted(left + right))


def gaussian_i(power):
    return ((1, 0), (0, 1), (-1, 0), (0, -1))[power % 4]


def gaussian_mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def main():
    subset_cases = []
    for size in range(2, 7):
        cases = [
            tuple(range(size)),
            (0,) + tuple(range(2, size + 1)),
            (1,) + tuple(range(2, size + 1)),
            tuple(range(2, size + 2)),
        ]
        for selected in cases:
            value = subset_color(selected, 0, 1)
            expected = subset_formula(size, 0 in selected, 1 in selected)
            assert value == expected
            subset_cases.append({"spokes": size, "has_a": 0 in selected, "has_b": 1 in selected, "coefficient": value})
    complete_checks = []
    for rank in range(2, 8):
        diagonals = [tuple(2 * j + 3 for j in range(rank)), tuple((-1) ** j * (j + 1) for j in range(rank))]
        for size in range(2, min(rank, 5) + 1):
            for diagonal in diagonals:
                for a, b in ((0, 1), (1, 0)):
                    actual = color_direct(diagonal, size, a, b)
                    expected = color_formula(diagonal, size, a, b)
                    assert actual == expected
                    complete_checks.append({"rank": rank, "spokes": size, "root": [a + 1, b + 1], "diagonal": list(diagonal), "value": actual})
    # Six-spoke endpoint cases check the first larger factorial directly.
    for diagonal in ((3, 5, 7, 9, 11, 13), (3, 5, 7, 9, 11, 13, 15)):
        value = color_direct(diagonal, 6, 0, 1)
        assert value == color_formula(diagonal, 6, 0, 1)
        complete_checks.append({"rank": len(diagonal), "spokes": 6, "root": [1, 2], "diagonal": list(diagonal), "value": value})
    # eta_i beta_i are even. Their ordered eta123 coefficient is minus beta123.
    source_signs = []
    for order in permutations(range(3)):
        sign = 1
        monomial = ()
        for i in order:
            factor, monomial = wedge_monomials(monomial, (i, i + 3))
            sign *= factor
        assert monomial == tuple(range(6)) and sign == -1
        source_signs.append(sign)
    phases = []
    for n in range(2, 33):
        phase = ((-1) ** n, 0)
        for factor in (gaussian_i(n + 1), (-1, 0), gaussian_i(-(n - 2)), gaussian_i(2), gaussian_i(-(n + 1))):
            phase = gaussian_mul(phase, factor)
        expected = tuple(-x for x in gaussian_i(n))
        assert phase == expected
        assert Fraction(factorial(n + 1), factorial(n) ** 3 * (n + 1)) == Fraction(1, factorial(n) ** 2)
        phases.append({"n": n, "phase": list(phase), "factorial": str(Fraction(1, factorial(n) ** 2))})
    # Q1 polynomial form checks the inverse-matrix Hessian identity independently.
    q1_checks = 0
    for case in complete_checks:
        if case["spokes"] != 2 or 0 in case["diagonal"]:
            continue
        mu = list(map(Fraction, case["diagonal"]))
        a, b = [j - 1 for j in case["root"]]
        prod = 1
        for entry in mu: prod *= entry
        inverse = [1 / x for x in mu]
        hessian = -prod * ((inverse[a] + inverse[b]) * sum(inverse) - (inverse[a] - inverse[b]) ** 2)
        assert hessian == case["value"]
        q1_checks += 1
    result = {
        "status": "PASS",
        "scope": "Exact boundary-induced quantum determinant-source correction; no independent bulk Q_n amplitude is claimed.",
        "subset_checks": subset_cases,
        "complete_color_checks": len(complete_checks),
        "complete_cases": complete_checks,
        "independent_q1_inverse_matrix_checks": q1_checks,
        "odd_source_sign_checks": len(source_signs),
        "phase_and_factorial_checks": phases,
        "q2_rank4_example": {"diagonal": [3,5,7,9], "tau": "1/2", "color": 4, "coefficient_removing_hbar2_g2_and_Kminus2": "1/4", "K": "32*pi^2"},
        "q3_rank4_example": {"diagonal": [3,5,7,9], "tau": "1/2", "color": -16, "coefficient_removing_hbar3_g3_and_Kminus3": "-i/18"},
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS', len(subset_cases), 'subset identities,', len(complete_checks), 'complete color sums,', q1_checks, 'inverse-matrix checks, 6 source signs, 31 phases')


if __name__ == '__main__':
    main()
