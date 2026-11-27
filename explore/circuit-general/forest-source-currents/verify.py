"""Independent rational finite-difference and complete-spinor certificate checks."""
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def clean(values):
    return {key: value for key, value in values.items() if value}


def add_labels(result, values, weight):
    for label, value in values.items():
        result[label] = result.get(label, F()) + weight * F(value)


def falling(n):
    result = [1]
    for k in range(n):
        following = [0] * (len(result) + 1)
        for i, value in enumerate(result):
            following[i] -= k * value
            following[i + 1] += value
        result = following
    return result


def newton_coefficients(samples, points, basis):
    by_point = dict(zip(map(tuple, points), samples))
    polynomial = {tuple(power): {} for power in basis}
    for i, j in basis:
        difference = {}
        for a in range(i + 1):
            for b in range(j + 1):
                add_labels(difference, by_point[a, b],
                           F((-1) ** (i + j - a - b) * comb(i, a) * comb(j, b),
                             factorial(i) * factorial(j)))
        for a, first in enumerate(falling(i)):
            for b, second in enumerate(falling(j)):
                if first and second:
                    add_labels(polynomial[a, b], difference, first * second)
    return [clean(polynomial[tuple(power)]) for power in basis]


def rational_rows(rows):
    return [{label: F(value) for label, value in row.items()} for row in rows]


def polynomial_add(target, values, basis, coefficient=1, shift=(0, 0)):
    for (i, j), value in zip(basis, values):
        power = (i + shift[0], j + shift[1])
        for label, number in value.items():
            key = power, label
            target[key] = target.get(key, F()) + coefficient * F(number)


def verify():
    document = json.loads((HERE / "certificate.json").read_text())
    if not document["complete"] or document["residual_periods"]:
        raise ArithmeticError("the source-current family is incomplete")
    edges = [[0, 1], [1, 3], [0, 2], [2, 4], [0, 5], [0, 6],
             [3, 5], [3, 6], [4, 5], [4, 6], [5, 6]]
    expected_basis = [[i, degree - i] for degree in range(7) for i in range(degree + 1)]
    if document["edges"] != edges or document["source_support"] != [1, 2, 5, 6]:
        raise ArithmeticError("the original graph or supported-source scope changed")
    if document["vertices"] != 7 or document["marked"] != [5, 6]:
        raise ArithmeticError("the original marked graph changed")
    if document["structurally_zero_source_vertices"] != [0, 3, 4]:
        raise ArithmeticError("a structural source constraint changed")
    from certify import dependency_hashes, verification_hashes
    if document["dependency_sha256"] != dependency_hashes():
        raise ArithmeticError("an exact dependency or the mandatory dependency list changed")
    if document["verification_sha256"] != verification_hashes():
        raise ArithmeticError("a complete replay implementation changed")
    points, basis = document["interpolation_points"], document["interpolation_monomials"]
    samples, reconstruction = document["samples"], document["reconstruction"]
    if points != expected_basis or basis != expected_basis or len(samples) != 28:
        raise ArithmeticError("a complete interpolation node is missing")
    if document["direct_points"] != [point for point in expected_basis if point[0] <= point[1]]:
        raise ArithmeticError("the independently integrated source-node list changed")
    if reconstruction["frame_monomials"] != expected_basis:
        raise ArithmeticError("the stored frame polynomial basis changed")
    for (a, b), sample in zip(points, samples):
        expected_sources = [["0", "0"], ["1", "0"], ["0", "1"], ["0", "0"],
                            ["0", "0"], [str(a), str(b)], [str(-1 - a), str(-1 - b)]]
        if sample["full_balanced_sources"] != expected_sources:
            raise ArithmeticError("a full interpolation source assignment changed")
        if not sample["all_values_finite"] or sample["residual_periods"]:
            raise ArithmeticError("an interpolation sample retained a scalar residual")
    bare = newton_coefficients([row["bare"] for row in samples], points, basis)
    if bare != rational_rows(reconstruction["bare_frame_coefficients"]):
        raise ArithmeticError("independent finite differences changed a bare coefficient")
    if any(value for power, value in zip(basis, bare) if sum(power) == 6):
        raise ArithmeticError("the complete bare polynomial exceeded its degree bound")
    current = []
    for edge in range(11):
        spinor = []
        for entry in range(2):
            found = newton_coefficients([row["all_original_first_current_spinors"][edge][entry]
                                         for row in samples], points, basis)
            if found != rational_rows(reconstruction["full_current_frame_coefficients"][edge][entry]):
                raise ArithmeticError("independent finite differences changed a complete current coefficient")
            spinor.append(found)
        current.append(spinor)
    independently_lifted = []
    for first, second in current:
        first_by_power = dict(zip(map(tuple, basis), first))
        second_by_power = dict(zip(map(tuple, basis), second))
        if first_by_power[0, 6] or second_by_power[6, 0]:
            raise ArithmeticError("the highest current entries lack the required coordinate factors")
        high = {}
        for i in range(6):
            j = 5 - i
            coefficient = first_by_power[i + 1, j]
            if coefficient != second_by_power[i, j + 1]:
                raise ArithmeticError("the complete degree-six identity b*f6=a*g6 failed")
            if coefficient:
                high[i, j] = coefficient
        independently_lifted.append([
            [[list(power), values] for power, values in zip(basis, entry)
             if sum(power) <= 5 and values] for entry in (first, second)] +
            [[[list(power), values] for power, values in sorted(high.items())]])
    stored_lift = [[[ [power, {label: F(value) for label, value in values.items()}]
                      for power, values in direction] for direction in edge]
                   for edge in reconstruction["denominator_free_original_source_covariants"]]
    if independently_lifted != stored_lift:
        raise ArithmeticError("a denominator-free full source-spinor coefficient differs from the independent lift")
    reconstruction_hash = hashlib.sha256(json.dumps(reconstruction, sort_keys=True,
                                                    separators=(",", ":")).encode()).hexdigest()
    if reconstruction_hash != document["complete_reconstruction_sha256"]:
        raise ArithmeticError("the complete reconstructed coefficient archive changed")
    for vertex in range(7):
        for entry in range(2):
            difference = {}
            for edge, (left, right) in enumerate(document["edges"]):
                polynomial_add(difference, current[edge][entry], basis,
                               int(left == vertex) - int(right == vertex))
            if vertex == 1 and entry == 0 or vertex == 2 and entry == 1:
                polynomial_add(difference, bare, basis, -1)
            if vertex == 5:
                polynomial_add(difference, bare, basis, -1, (1, 0) if entry == 0 else (0, 1))
            if vertex == 6:
                polynomial_add(difference, bare, basis, 1)
                polynomial_add(difference, bare, basis, 1, (1, 0) if entry == 0 else (0, 1))
            if clean(difference):
                raise ArithmeticError("a complete polynomial source identity failed")
    for entry in [bare] + [row for spinor in current for row in spinor]:
        if any(set(value) - {"1", "z3", "z5"} for value in entry):
            raise ArithmeticError("a coefficient left the recorded finite period channels")
    from independent_controls import check
    check()
    from source_currents import evaluate, bracket, add, record
    physical = ((0, 0), (1, 2), (3, -1), (0, 0), (0, 0), (5, 3), (-9, -4))
    expected = {"1": "3327759589/46656000", "z3": "50942017/324000", "z5": "-45619/180"}
    if evaluate(reconstruction, physical)["bare"] != expected:
        raise ArithmeticError("the independent earlier physical-quartic value changed")
    original = ((0, 0), (1, 0), (0, 1), (0, 0), (0, 0), (F(2, 3), F(-5, 2)),
                (F(-5, 3), F(3, 2)))
    value = evaluate(reconstruction, original)
    for matrix in (((2, 1), (1, 2)), ((0, 1), (1, 0)), ((3, 0), (0, -2))):
        determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
        transformed = [tuple(sum(F(matrix[i][j]) * source[j] for j in range(2))
                             for i in range(2)) for source in original]
        actual = evaluate(reconstruction, transformed)
        if actual["bare"] != record(add([(determinant ** 5, value["bare"])])):
            raise ArithmeticError("the scalar determinant covariance failed")
        for edge in range(11):
            for entry in range(2):
                expected = record(add((determinant ** 5 * matrix[entry][j],
                                       value["all_original_first_current_spinors"][edge][j]) for j in range(2)))
                if actual["all_original_first_current_spinors"][edge][entry] != expected:
                    raise ArithmeticError("the full free-spinor covariance failed")
    print("Independent finite differences: all 644 frame coefficients; 14 polynomial source identities; "
          "all 11 complete spinors, fresh and degenerate controls, and nonunit source covariance pass.")


if __name__ == "__main__":
    verify()
