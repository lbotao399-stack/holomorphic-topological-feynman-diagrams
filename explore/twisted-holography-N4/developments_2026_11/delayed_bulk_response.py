#!/usr/bin/env python3
"""Exact four-loop onset of a finite descendant polarization in flat BCOV.

Standard library only. The existing complete planar genus-four periods are
read as data, without importing their graph engine. New donor cancellation,
field-level admissibility, Fourier factorials, and signs are checked here.
The low-genus test independently evaluates the genus-one stable partitions.
"""

from collections import Counter
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import factorial, gcd, prod
from pathlib import Path
import ast
import json


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = HERE.parent / "code/planar_four_loop.py"


def weak(total, length):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in weak(total - first, length - 1):
            yield (first, *rest)


def rank(rows, width):
    rows = [[F(x) for x in row] for row in rows]
    pivot = 0
    for column in range(width):
        chosen = next((i for i in range(pivot, len(rows)) if rows[i][column]), None)
        if chosen is None:
            continue
        rows[pivot], rows[chosen] = rows[chosen], rows[pivot]
        leading = rows[pivot][column]
        rows[pivot] = [x / leading for x in rows[pivot]]
        for i, row in enumerate(rows):
            if i != pivot and row[column]:
                factor = row[column]
                rows[i] = [x - factor * y for x, y in zip(row, rows[pivot])]
        pivot += 1
    return pivot


def multiply(left, right):
    result = Counter()
    for a, x in left.items():
        for b, y in right.items():
            result[tuple(i + j for i, j in zip(a, b))] += x * y
    return {a: x for a, x in result.items() if x}


def polynomial_power(poly, exponent):
    result = {(0, 0, 0): 1}
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def derivative(poly):
    result = Counter()
    for a, coefficient in poly.items():
        for i, value in enumerate(a):
            if value:
                b = tuple(x - (j == i) for j, x in enumerate(a))
                result[b] += value * coefficient
    return {a: x for a, x in result.items() if x}


def normalized_coefficients(poly):
    coefficients = {a: int(c) * prod(map(factorial, a)) for a, c in poly.items()}
    divisor = 0
    for value in coefficients.values():
        divisor = gcd(divisor, abs(value))
    return {a: c // divisor for a, c in coefficients.items()}, divisor


def classification():
    quadratic = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1,
                 (1, 1, 0): -1, (1, 0, 1): -1, (0, 1, 1): -1}
    cubic = multiply(multiply({(1, 0, 0): 2, (0, 1, 0): -1, (0, 0, 1): -1},
                              {(1, 0, 0): -1, (0, 1, 0): 2, (0, 0, 1): -1}),
                     {(1, 0, 0): -1, (0, 1, 0): -1, (0, 0, 1): 2})
    vandermonde = multiply(multiply({(1, 0, 0): 1, (0, 1, 0): -1},
                                    {(0, 1, 0): 1, (0, 0, 1): -1}),
                           {(0, 0, 1): 1, (1, 0, 0): -1})
    polynomials = {0: {(0, 0, 0): 1}, 2: quadratic, 3: cubic,
                   4: polynomial_power(quadratic, 2),
                   6: polynomial_power(vandermonde, 2)}
    rows = []
    normalized = {}
    for total in range(13):
        orbits = sorted({tuple(sorted(a, reverse=True)) for a in product(range(5), repeat=3)
                         if sum(a) == total}, reverse=True)
        equations = []
        for remaining in weak(total - 1, 3) if total else ():
            if max(remaining) > 4:
                continue
            row = []
            for orbit in orbits:
                row.append(sum(tuple(sorted(tuple(x + (i == j) for j, x in enumerate(remaining)),
                                            reverse=True)) == orbit for i in range(3)))
            equations.append(row)
        dimension = len(orbits) - rank(equations, len(orbits))
        assert dimension == int(total in polynomials), (total, dimension)
        if dimension:
            poly = polynomials[total]
            assert not derivative(poly)
            assert all(max(a) <= 4 for a in poly)
            coeff, divisor = normalized_coefficients(poly)
            for a, c in coeff.items():
                assert all(coeff[b] == c for b in set(permutations(a)))
            normalized[total] = coeff
        rows.append({"total_level": total, "symmetric_space_dimension": len(orbits),
                     "donor_cancellation_rank": len(orbits) - dimension,
                     "kernel_dimension": dimension})
    return rows, normalized


def donor_check(coefficients, total):
    tests = 0
    assignments = 0
    for count in range(1, total + 1):
        for remaining in weak(total - count, 3):
            value = 0
            for places in product(range(3), repeat=count):
                after = tuple(k + places.count(i) for i, k in enumerate(remaining))
                value += coefficients.get(after, 0)
                assignments += 1
            assert value == 0, (count, remaining, value)
            tests += 1
    return {"nonempty_donor_counts": list(range(1, total + 1)),
            "residual_level_tests": tests, "individual_labeled_assignments": assignments}


def read_planar_data():
    """Read constants from an AST; never execute or import the supplied engine."""
    def evaluate(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return node.value
        if isinstance(node, (ast.List, ast.Tuple)):
            return tuple(evaluate(x) for x in node.elts)
        if isinstance(node, ast.Dict):
            return {evaluate(k): evaluate(v) for k, v in zip(node.keys, node.values)}
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "Fraction" and not node.keywords):
            return F(*(evaluate(x) for x in node.args))
        raise ValueError(ast.dump(node))
    assignments = {}
    for node in ast.parse(SOURCE.read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in ("LEVELS", "EXPECTED"):
                    assignments[target.id] = evaluate(node.value)
    return dict(zip(assignments["LEVELS"], assignments["EXPECTED"][4]))


def determinant(matrix):
    """Bareiss determinant, with exact integer divisions."""
    if not matrix:
        return 1
    matrix = [list(map(int, row)) for row in matrix]
    size = len(matrix)
    previous, sign = 1, 1
    for k in range(size - 1):
        pivot = next((i for i in range(k, size) if matrix[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            sign = -sign
        divisor = matrix[k][k]
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                value = matrix[i][j] * divisor - matrix[i][k] * matrix[k][j]
                assert value % previous == 0
                matrix[i][j] = value // previous
            matrix[i][k] = 0
        previous = divisor
    return sign * matrix[-1][-1]


def fiber_degree(vectors):
    genus, remainder = divmod(len(vectors), 3)
    assert not remainder and all(len(p) == 3 for p in vectors)
    triples = {sum(1 << i for i in group): determinant([vectors[i] for i in group]) ** 2
               for group in combinations(range(len(vectors)), 3)}

    @cache
    def unordered(mask):
        if not mask:
            return 1
        first = mask & -mask
        return sum(value * unordered(mask ^ triple) for triple, value in triples.items()
                   if value and triple & first and triple & mask == triple)
    return factorial(genus) * unordered((1 << len(vectors)) - 1)


def gram_finite_difference(vectors):
    genus = len(vectors) // 3
    answer = 0
    for mask in range(1 << len(vectors)):
        rows = [p for i, p in enumerate(vectors) if mask >> i & 1]
        gram = [[sum(p[i] * p[j] for p in rows) for j in range(3)] for i in range(3)]
        answer += (-1) ** (len(vectors) - mask.bit_count()) * determinant(gram) ** genus
    return answer


def four_loop_fourier():
    axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    surviving = []
    total = F()
    for a, b, c in weak(10, 3):
        root = (4 - a, 3 - b, 4 - c)
        vectors = (*((axes[0],) * a), *((axes[1],) * (b + 1)),
                   *((axes[2],) * c), root)
        degree = fiber_degree(vectors)
        if degree:
            surviving.append({"primary_scalar_counts": [a, b, c], "root": list(root),
                              "fiber_degree": degree, "root_vector_pairing": root[0]})
        total += F(root[0] * degree, factorial(a) * factorial(b) * factorial(c))
    assert {tuple(row["primary_scalar_counts"]) for row in surviving} == {(3, 3, 4), (4, 2, 4), (4, 3, 3)}
    assert all(row["fiber_degree"] == factorial(4) ** 3 for row in surviving)
    assert total == 16
    generic = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0),
               (0, 1, 1), (1, 0, 1), (2, 1, 1), (-1, 2, 1),
               (1, 1, -2), (2, -1, 1), (1, 2, 3), (3, 1, -1))
    degree = fiber_degree(generic)
    assert degree == gram_finite_difference(generic) and degree > 0
    return {"all_primary_count_triples_checked": 66, "surviving": surviving,
            "root_weighted_factorial_sum": str(total),
            "independent_twelve_vector_gram_check": degree}


@cache
def partitions(count, groups):
    if count == 0:
        return ((),) if groups == 0 else ()
    if groups == 0 or groups > count:
        return ()
    answer = []
    for partition in partitions(count - 1, groups):
        for i in range(groups):
            answer.append((*partition[:i], (*partition[i], count - 1), *partition[i + 1:]))
    for partition in partitions(count - 1, groups - 1):
        answer.append((*partition, (count - 1,)))
    return tuple(answer)


@cache
def genus_one_psi(levels):
    levels = tuple(sorted(levels))
    if not levels or sum(levels) != len(levels):
        return F()
    if levels[0]:
        return F(factorial(len(levels) - 1), 24)
    rest = levels[1:]
    return sum((genus_one_psi(tuple(sorted((*rest[:i], k - 1, *rest[i + 1:]))))
                for i, k in enumerate(rest) if k), F())


def genus_one_vertex(momenta, levels):
    assert len(momenta) == len(levels) and not any(map(sum, zip(*momenta)))
    dimension = len(momenta[0])
    assert sum(levels) == len(momenta) - dimension
    answer = F()
    contributing = 0
    for blocks in partitions(len(momenta), dimension + 1):
        ordinary, elliptic = [], []
        for block in blocks:
            degree = sum(levels[i] for i in block)
            if degree == len(block) - 1:
                ordinary.append(block)
            elif degree == len(block):
                elliptic.append(block)
            else:
                break
        else:
            if len(elliptic) != 1 or len(ordinary) != dimension:
                continue
            weight = genus_one_psi(tuple(levels[i] for i in elliptic[0]))
            for block in ordinary:
                weight *= F(factorial(len(block) - 1),
                            prod(factorial(levels[i]) for i in block))
            rows = [tuple(sum(momenta[i][j] for i in block) for j in range(dimension))
                    for block in ordinary]
            value = weight * determinant(rows) ** 2
            answer += value
            contributing += int(value != 0)
    return answer, contributing


def one_loop_control():
    # A full-rank fivefold input, with (a,b,c)=(1,1,2) and balancing root (3,2,2).
    momenta = ((1, 0, -4, -4, -4), (0, 1, 0, 0, 0), (-1, -1, 0, 0, 0),
               (0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 0, 1),
               (0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 3, 2, 2))
    orbits = (((4, 0, 0), 2), ((3, 1, 0), -1), ((2, 2, 0), 1))
    entries, total = [], F()
    expected = {(4, 0, 0): F(7, 3), (3, 1, 0): F(28, 3), (2, 2, 0): F(14)}
    for orbit, coefficient in orbits:
        for row in sorted(set(permutations(orbit))):
            value, count = genus_one_vertex(momenta, (*row, *(0 for _ in momenta[3:])))
            assert value == expected[orbit]
            entries.append({"levels": list(row), "weight": coefficient,
                            "ordinary_positive_vertex": str(value),
                            "nonzero_stable_partitions": count})
            total += coefficient * value
    assert len(entries) == 12 and total == 0
    return {"all_distinct_weighted_tuples_checked": len(entries), "entries": entries,
            "sum": str(total),
            "method": "all twelve tuples evaluated independently by stable partitions and integer determinants"}


def main():
    dimensions, coefficients = classification()
    donor = {str(k): donor_check(coefficients[k], k) for k in (2, 3, 4, 6)}
    data = read_planar_data()
    combined = sum((value * data[tuple(sorted(k, reverse=True))]
                    for k, value in coefficients[4].items()), F())
    assert combined == -F(31, 1024192512000)
    classical = sum((F(value * factorial(4), prod(map(factorial, k)))
                     for k, value in coefficients[4].items()), F())
    assert classical == 0
    fourier = four_loop_fourier()
    heat_stripped = -combined * F(fourier["root_weighted_factorial_sum"])
    assert heat_stripped == F(31, 64012032000)
    result = {
        "status": "PASS",
        "arithmetic": "exact integers and rational numbers; standard library only",
        "scope": "selected polarized local elementary-field response in fixed flat BCOV coordinates",
        "source_revision": "aa0d2f1725dff135494eef132f5e92a8c0401a42",
        "minimal_scalar_levels": [0, 1, 2, 3, 4],
        "classification": dimensions,
        "admissible_symmetric_total_levels": [0, 2, 3, 4, 6],
        "degree_four_weights_per_distinct_permutation": {
            "4,0,0": 2, "3,1,0": -1, "2,2,0": 1, "2,1,1": 0},
        "labeled_donor_checks": donor,
        "classical_multinomial_sum": str(classical),
        "lower_genus_vertical_rank": [
            {"genus": g, "forgotten_marks": 2*g+4, "complex_rank_bound": 3*g,
             "strict_deficit": 4-g} for g in (1, 2, 3)],
        "one_loop_direct_control": one_loop_control(),
        "inherited_complete_planar_periods": {",".join(map(str, k)): str(v) for k, v in data.items()},
        "new_planar_polarization": str(combined),
        "four_loop_fourier": fourier,
        "labeled_scalar_vertex": str(combined * factorial(4)**3),
        "labeled_elementary_response": str(-combined * factorial(4)**3),
        "fixed_fourier_heat_stripped_response": str(heat_stripped),
        "heat_normalization": "(2*pi*i)^20 Q4 phi; elementary sign -(-1)^4; root pairing +1",
        "degree_six_status": "all lower genera vanish; first possible genus-six coefficient not evaluated",
        "inherited_input_sha256": {str(SOURCE.relative_to(ROOT)): sha256(SOURCE.read_bytes()).hexdigest()},
        "not_recomputed": "the full genus-four master-period integrations",
    }
    destination = HERE / "delayed_bulk_response.json"
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
