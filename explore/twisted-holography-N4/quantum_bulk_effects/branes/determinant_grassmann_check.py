#!/usr/bin/env python3
"""Independent Grassmann-matrix verification of all-order determinant colors.

No inherited fan/color engine is imported. This checks the determinant color
identity, not an independently computed quantum bulk amplitude.
"""
from itertools import combinations
from math import factorial
from pathlib import Path
import json


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for mask, coefficient in polynomial.items():
            result[mask] = result.get(mask, 0) + coefficient
    return {mask: coefficient for mask, coefficient in result.items() if coefficient}


def scale(polynomial, value):
    return {mask: coefficient * value for mask, coefficient in polynomial.items()
            if coefficient * value}


def multiply(left, right):
    result = {}
    for a, ca in left.items():
        for b, cb in right.items():
            if a & b:
                continue
            inversions = sum((b & ((1 << j) - 1)).bit_count()
                             for j in range(a.bit_length()) if a & (1 << j))
            mask = a | b
            result[mask] = result.get(mask, 0) + (-1) ** inversions * ca * cb
    return {mask: coefficient for mask, coefficient in result.items() if coefficient}


def matrix_multiply(left, right):
    rank = len(left)
    return [[add(*(multiply(left[i][k], right[k][j]) for k in range(rank)))
             for j in range(rank)] for i in range(rank)]


def matrix_add(left, right, sign=1):
    return [[add(a, scale(b, sign)) for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(left, right)]


def matrix_scale(polynomial, value):
    return [[multiply(polynomial, entry) for entry in row] for row in value]


def scalar_product(values):
    value = 1
    for entry in values:
        value *= entry
    return value


def elementary(values, degree):
    if degree < 0 or degree > len(values):
        return 0
    return sum(scalar_product(part) for part in combinations(values, degree))


def gaussian_integral(polynomial, eigenvalues):
    """Top coefficient after multiplication by product(1+mu_i barpsi_i psi_i)."""
    result = 0
    for mask, coefficient in polynomial.items():
        remaining = []
        for i, eigenvalue in enumerate(eigenvalues):
            pair = (mask >> (2 * i)) & 3
            if pair in (1, 2):
                break
            if pair == 0:
                remaining.append(eigenvalue)
        else:
            # Each missing pair has even degree, hence no reordering sign.
            result += coefficient * scalar_product(remaining)
    return result


def closed_core(eigenvalues, a, b, spokes):
    loop = spokes - 1
    others = [x for i, x in enumerate(eigenvalues) if i not in (a, b)]
    if spokes % 2:
        return (factorial(loop) * (eigenvalues[b] - eigenvalues[a])
                * elementary(others, len(eigenvalues) - loop - 2))
    return (-factorial(loop) * (eigenvalues[a] + eigenvalues[b])
            * elementary(others, len(eigenvalues) - loop - 2)
            - 2 * (loop + 1) * factorial(loop - 1)
            * elementary(others, len(eigenvalues) - loop - 1))


def main():
    matrix_checks = 0
    integration_checks = 0
    samples = []
    for rank in range(2, 8):
        # Generator order is barpsi_1,psi_1,...,barpsi_N,psi_N.
        bilinear = [[multiply({1 << (2 * i + 1): 1}, {1 << (2 * j): 1})
                     for j in range(rank)] for i in range(rank)]
        scalar = add(*({3 << (2 * i): 1} for i in range(rank)))
        assert matrix_multiply(bilinear, bilinear) == matrix_scale(scalar, bilinear)
        matrix_checks += 1
        for a, b in ((0, 1), (rank - 1, 0)):
            unit = [[{} for _ in range(rank)] for _ in range(rank)]
            unit[a][b] = {0: 1}
            left = matrix_multiply(bilinear, unit)
            right = matrix_multiply(unit, bilinear)
            middle = matrix_multiply(left, bilinear)
            power = [{0: 1}]
            for _ in range(rank + 1):
                power.append(multiply(power[-1], scalar))
            iterated = unit
            for spokes in range(1, rank + 2):
                iterated = matrix_add(matrix_multiply(bilinear, iterated),
                                      matrix_multiply(iterated, bilinear), -1)
                if spokes == 1:
                    continue
                expected = matrix_add(
                    matrix_scale(power[spokes - 1], matrix_add(left, right, (-1) ** spokes)),
                    matrix_scale(power[spokes - 2], middle), -(1 + (-1) ** spokes))
                assert iterated == expected
                matrix_checks += 1
                assignments = [list(range(1, rank + 1)), [2] * rank,
                               [0 if i == 2 else (-1) ** i * (i + 1) for i in range(rank)]]
                for eigenvalues in assignments:
                    # Completeness gives sum t_A(barpsi t_A psi)=-tau B.
                    actual = (-1) ** spokes * gaussian_integral(iterated[a][b], eigenvalues)
                    wanted = closed_core(eigenvalues, a, b, spokes)
                    assert actual == wanted, (rank, a, b, spokes, eigenvalues, actual, wanted)
                    integration_checks += 1
                    if rank == 4 and (a, b) == (0, 1) and eigenvalues == [1, 2, 3, 4]:
                        samples.append({"rank": rank, "spokes": spokes,
                                        "eigenvalues": eigenvalues, "color_core": actual})
    # Rank thresholds and every loop order in a longer independent scalar check.
    threshold_checks = 0
    for loop in range(2, 41):
        rank = loop + (2 if loop % 2 == 0 else 1)
        values = list(range(1, rank + 1))
        assert closed_core(values, 0, 1, loop + 1) != 0
        assert closed_core(values[:-1], 0, 1, loop + 1) == 0
        threshold_checks += 2
    report = {
        "status": "PASS",
        "scope": "Independent exact Grassmann-matrix proof check of determinant fan colors; not a quantum bulk amplitude.",
        "ranks_directly_checked": list(range(2, 8)),
        "matrix_identity_checks": matrix_checks,
        "gaussian_integral_checks": integration_checks,
        "rank_threshold_checks": threshold_checks,
        "samples": samples,
        "convention": "barpsi_i before psi_i; normalized Gaussian integral equals det(M); color core divides out tau^spokes",
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
