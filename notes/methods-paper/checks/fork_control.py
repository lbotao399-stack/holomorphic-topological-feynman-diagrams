#!/usr/bin/env python3
"""Independently replay the complete six-variable three-loop fork control.

Constructs all seven original current spinors, their complete chord-current
Jacobian, every polynomial moment and every Kirchhoff identity.  An explicit
positive resistance map supplies independent original-network and orientation
controls.  The default recomputes and compares the stored exact certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import sys

import sympy as sp


PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT / "code"))
from htk import brute  # noqa: E402
from htk.graph import Graph  # noqa: E402


VARIABLES = sp.symbols("r1:7")
r1, r2, r3, r4, r5, r6 = VARIABLES
EDGES = ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
EDGE_NAMES = ("12", "13", "2A", "2B", "3A", "3B", "AB")
SOURCES = ((1, 0), (0, 0), (0, 0), (0, 1), (-1, -1))
TREE = (0, 1, 2, 3)
COTREE = (4, 5, 6)


def rational(value):
    numerator, denominator = sp.cancel(value).as_numer_denom()
    return Fraction(int(numerator), int(denominator))


def polynomial_data(polynomial):
    return [{"powers": list(powers), "coefficient": str(rational(coefficient))}
            for powers, coefficient in sp.Poly(sp.expand(polynomial), *VARIABLES).terms()
            if coefficient]


def integrate(polynomial):
    """Two intervals and two ordinary two-simplices; every monomial kept."""
    total = Fraction(0)
    for powers, coefficient in sp.Poly(sp.expand(polynomial), *VARIABLES).terms():
        term = rational(coefficient)
        term *= Fraction(1, (powers[0] + 1) * (powers[1] + 1))
        for first, second in ((2, 3), (4, 5)):
            term *= Fraction(factorial(powers[first]) * factorial(powers[second]),
                             factorial(powers[first] + powers[second] + 2))
        total += term
    return total


def currents():
    first = sp.Matrix((r5, -r6))
    second = sp.Matrix((1 - r5, r6))
    chord = sp.Matrix((r5, 1 - r6))
    return (r1 * first + r2 * second,
            (1 - r1) * first + (1 - r2) * second,
            r1 * first - r3 * chord,
            r2 * second + r3 * chord,
            (1 - r1) * first - r4 * chord,
            (1 - r2) * second + r4 * chord,
            (1 - r3 - r4) * chord)


def resistances():
    """Invert both star--triangle maps and all three parallel splittings."""
    reduced_1A, reduced_1B, reduced_AB = 1 - r5 - r6, r5, r6
    triangle_first = (reduced_1A / r1, reduced_1B / r2, reduced_AB / r3)
    triangle_second = (reduced_1A / (1 - r1), reduced_1B / (1 - r2), reduced_AB / r4)

    def star(triangle):
        a, b, c = triangle
        total = a + b + c
        return tuple(sp.factor(value) for value in (a * b / total,
                                                    a * c / total,
                                                    b * c / total))

    first, second = star(triangle_first), star(triangle_second)
    return (first[0], second[0], first[1], first[2], second[1], second[2],
            reduced_AB / (1 - r3 - r4))


def point_control(point, source_1, source_A, all_currents, density, times, time_derivatives):
    substitution = dict(zip(VARIABLES, map(sp.Rational, point)))
    exact_times = [rational(value.subs(substitution)) for value in times]
    assert all(value > 0 for value in exact_times)
    full_sources = (tuple(source_1), (0, 0), (0, 0), tuple(source_A),
                    tuple(-a - b for a, b in zip(source_1, source_A)))
    graph = Graph(5, EDGES)
    original = brute.currents(graph, exact_times, full_sources)
    expected = []
    for current in all_currents:
        coefficients = [rational(value.subs(substitution)) for value in current]
        expected.append([coefficients[0] * source_1[a] + coefficients[1] * source_A[a]
                         for a in range(2)])
    assert original == expected
    time_matrix = sp.Matrix([[sp.Rational(exact_times[edge].numerator,
                                          exact_times[edge].denominator)]
                             + [value.subs(substitution)
                                for value in time_derivatives[edge, :]]
                            for edge in range(7)])
    projective_jacobian = rational(time_matrix.det())
    assert projective_jacobian > 0
    bracket = source_1[0] * source_A[1] - source_1[1] * source_A[0]
    chart_density = rational(density.subs(substitution)) * bracket ** 3
    original_density = brute.density(graph, exact_times, full_sources).rational()
    assert original_density * projective_jacobian == chart_density
    return {"ratio_point": list(map(str, point)),
            "full_balanced_sources": [list(row) for row in full_sources],
            "source_bracket": bracket, "original_resistances": list(map(str, exact_times)),
            "all_original_currents": [[str(value) for value in row] for row in original],
            "positive_projective_jacobian": str(projective_jacobian),
            "original_density": str(original_density),
            "complete_pulled_back_density": str(chart_density)}


def build():
    all_currents = currents()
    density = (r1 - r2) * (1 - r3 - r4) * r5 * (1 - r5 - r6)
    current_matrix = sp.Matrix([list(all_currents[edge][a].diff(variable)
                                     for variable in VARIABLES)
                                for edge in COTREE for a in range(2)])
    current_determinant = sp.factor(current_matrix.det())
    assert sp.expand(current_determinant + density) == 0
    assert integrate(density) == 0
    expected = [(-1, -2), (1, 2), (2, -1), (-3, -1), (-2, 1), (3, 1), (0, 0)]
    moments = [[integrate(density * value) for value in current] for current in all_currents]
    assert moments == [[Fraction(value, 8640) for value in row] for row in expected]
    polynomial_conservation = []
    integrated_conservation = []
    for vertex in range(5):
        for a in range(2):
            current_sum = sum((int(tail == vertex) - int(head == vertex)) * all_currents[edge][a]
                              for edge, (tail, head) in enumerate(EDGES))
            residual = sp.expand(current_sum - SOURCES[vertex][a])
            assert residual == 0
            integrated = sum((int(tail == vertex) - int(head == vertex)) * moments[edge][a]
                             for edge, (tail, head) in enumerate(EDGES))
            assert integrated == SOURCES[vertex][a] * integrate(density) == 0
            polynomial_conservation.append({"vertex": vertex, "source_direction": a,
                                             "complete_residual": "0"})
            integrated_conservation.append({"vertex": vertex, "source_direction": a,
                                              "complete_residual": str(integrated)})
    times = resistances()
    time_derivatives = sp.Matrix(times).jacobian(VARIABLES)
    controls = [point_control((Fraction(1, 3), Fraction(2, 3), Fraction(1, 4),
                               Fraction(1, 5), Fraction(1, 5), Fraction(1, 4)),
                              (1, 0), (0, 1), all_currents, density, times, time_derivatives),
                point_control((Fraction(2, 5), Fraction(3, 4), Fraction(1, 6),
                               Fraction(2, 7), Fraction(2, 7), Fraction(1, 3)),
                              (2, -1), (1, 3), all_currents, density, times, time_derivatives)]
    dependencies = [Path(__file__), PAPER / "sections/networks.tex",
                    ROOT / "code/htk/brute.py", ROOT / "code/htk/graph.py",
                    ROOT / "code/htk/CONVENTIONS.md"]
    return {"schema": "methods_paper_fork_control_v1",
            "original_vertex_labels": ["1", "2", "3", "A", "B"],
            "ordered_oriented_edges": [list(edge) for edge in EDGES],
            "edge_names": list(EDGE_NAMES), "marked": [3, 4],
            "source_scope": "lambda_1 and lambda_A arbitrary complete spinors; lambda_2=lambda_3=0; lambda_B=-lambda_1-lambda_A",
            "symbolic_basis_source_rows": [list(row) for row in SOURCES],
            "source_factor_in_density": "<lambda_1,lambda_A>^3",
            "variables": list(map(str, VARIABLES)),
            "domain": "(0,1)^2 times Delta^2 times Delta^2",
            "full_current_polynomials_in_original_source_directions": [
                {"edge": name, "lambda_1_coefficient": polynomial_data(current[0]),
                 "lambda_A_coefficient": polynomial_data(current[1])}
                for name, current in zip(EDGE_NAMES, all_currents)],
            "complete_density_polynomial": polynomial_data(density),
            "full_chord_current_jacobian": [[str(sp.expand(entry)) for entry in row]
                                            for row in current_matrix.tolist()],
            "complete_chord_current_jacobian_determinant": polynomial_data(current_determinant),
            "orientation": {"tree_edges": list(TREE), "cotree_edges": list(COTREE),
                            "current_area_sign": -1, "projective_chart_sign": 1,
                            "determinant_plus_density_polynomial": "0"},
            "complete_resistance_map": list(map(str, times)),
            "bare_kernel": "0",
            "all_first_current_spinors_after_bracket_cubed": [
                {"edge": name, "lambda_1": str(row[0]), "lambda_A": str(row[1])}
                for name, row in zip(EDGE_NAMES, moments)],
            "all_polynomial_kirchhoff_identities": polynomial_conservation,
            "all_integrated_kirchhoff_identities": integrated_conservation,
            "original_network_controls": controls,
            "dependencies_sha256": {str(path.relative_to(ROOT)):
                                    hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in dependencies}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = build()
    path = PAPER / "verification/fork_control.json"
    if args.write:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print("Wrote", path)
    else:
        assert json.loads(path.read_text()) == result, "fork certificate differs from exact replay"
    print("PASS: full six-variable current-area determinant, all seven free spinors,")
    print("      ten polynomial and ten integrated Kirchhoff identities, and two original-network orientation controls.")


if __name__ == "__main__":
    main()
