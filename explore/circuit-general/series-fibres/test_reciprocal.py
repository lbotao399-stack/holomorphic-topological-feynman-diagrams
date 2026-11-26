"""Independent exact controls for reciprocal support and CRT certificates.

Run: python3 explore/circuit-general/series-fibres/test_reciprocal.py
Small complete numerator dictionaries are obtained from the original
Fraction interpolator, then integrated monomial by monomial. They do not
use the reciprocal sampler or its interpolation grid.
"""

from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "cores"))

from htk import brute
from htk.graph import Graph, fan, polygon
from htk.series import SeriesReduction
from htk.fibres import integrate_fibre
from reciprocal import reciprocal_plan, reciprocal_gauge_counts, integer_certificate
import numer


def fixtures():
    first = fan(2)
    first_sources = [[i + 1, (-1) ** i * (i + 2)] for i in range(first.nv - 1)]
    yield SeriesReduction(first, first_sources)
    graph = Graph(5, [(0, 4), (4, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
    yield SeriesReduction(graph, [[1, 2], [3, -2], [-1, 1], [2, 4]])


def original_reciprocal(reduction):
    """Reference push-forward of the old full Fraction coefficient dictionary."""
    polynomial = numer.bare_numerator(reduction.G, reduction.full_sources)
    result = {}
    for exponent, coefficient in polynomial.items():
        reciprocal = []
        weight = coefficient
        for chain in reduction.chains:
            degree = sum(exponent[edge] for edge in chain.edges)
            reciprocal.append(reduction.n - 1 - degree)
            numerator = 1
            for edge in chain.edges:
                numerator *= factorial(exponent[edge])
            weight *= F(numerator, factorial(degree + chain.length - 1))
        key = tuple(reciprocal)
        result[key] = result.get(key, F(0)) + weight
    return {key: value for key, value in result.items() if value}


def value_at(polynomial, point):
    total = F(0)
    for exponent, coefficient in polynomial.items():
        value = coefficient
        for coordinate, power in zip(point, exponent):
            value *= F(coordinate) ** power
        total += value
    return total


class ReciprocalChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reductions = list(fixtures())
        cls.references = [original_reciprocal(reduction) for reduction in cls.reductions]

    def test_exact_support_and_period_powers(self):
        for reduction, reference in zip(self.reductions, self.references):
            self.assertTrue(reference)
            plan = reciprocal_plan(reduction)
            exponents = np.array(list(reference), dtype=np.int64)
            self.assertTrue(plan.support_mask(exponents).all())
            powers = plan.period_powers(exponents)
            self.assertTrue((powers > 0).all())
            self.assertTrue((powers.sum(axis=1)
                             == reduction.n * (reduction.n + reduction.H)).all())
            self.assertEqual(plan.total_degree,
                             (reduction.n - 1) * (plan.core_vertices - 2))

    def test_all_grid_coefficients_against_original(self):
        prime = 1000003
        for reduction, reference in zip(self.reductions, self.references):
            for gauge in range(len(reduction.chains)):
                plan = reciprocal_plan(reduction, gauge=gauge)
                samples = []
                for point in plan.sample_points():
                    value = value_at(reference, point)
                    samples.append(value.numerator * pow(value.denominator, -1, prime) % prime)
                coefficients = plan.lower_set.newton_to_monomial(np.array(samples, dtype=np.int64), prime)
                actual = {tuple(int(x) for x in exponent): int(coefficient)
                          for exponent, coefficient in zip(plan.exponents(), coefficients)
                          if coefficient}
                expected = {exponent: coefficient.numerator * pow(coefficient.denominator, -1, prime) % prime
                            for exponent, coefficient in reference.items()}
                self.assertEqual(actual, {key: value for key, value in expected.items() if value})
                self.assertTrue((coefficients[~plan.support_mask()] == 0).all())

    def test_complete_fibre_values_and_jacobian(self):
        for reduction, reference in zip(self.reductions, self.references):
            conductances = tuple(F(i + 2) for i in range(len(reduction.chains)))
            resistances = tuple(1 / value for value in conductances)
            fibre = integrate_fibre(reduction, resistances)
            centre = tuple(tuple(F(1, chain.length) for _ in chain.edges) for chain in reduction.chains)
            U = brute.symanzik_U(reduction.G, reduction.times(resistances, centre))
            expected = fibre.value * U ** (reduction.n + reduction.H)
            for conductance, resistance, chain in zip(conductances, resistances, reduction.chains):
                expected *= conductance ** (reduction.n - 1) / resistance ** (chain.length - 1)
            self.assertEqual(value_at(reference, conductances), expected)

    def test_integrality_height_and_rational_sources(self):
        for reduction, reference in zip(self.reductions, self.references):
            certificate = integer_certificate(reduction)
            json.dumps(certificate)  # All evidence is persistable without custom conversion.
            for coefficient in reference.values():
                scaled = coefficient * certificate["denominator"]
                self.assertEqual(scaled.denominator, 1)
                self.assertLessEqual(abs(scaled.numerator), certificate["coefficient_bound"])
            rational = SeriesReduction(reduction.G,
                                       [[value / 3 for value in row] for row in reduction.full_sources])
            rational_certificate = integer_certificate(rational)
            self.assertEqual(rational_certificate["source_denominator"], 3)
            self.assertEqual(rational_certificate["denominator"],
                             certificate["denominator"] * 3 ** (reduction.n * reduction.H))
            self.assertEqual(rational_certificate["coefficient_bound"], certificate["coefficient_bound"])

    def test_grid_nonsingularity_bound(self):
        for reduction in self.reductions:
            plan = reciprocal_plan(reduction)
            core = Graph(plan.core_vertices, plan.core_edges)
            for conductances in plan.sample_points():
                laplacian_determinant = det_direct(core, conductances)
                self.assertGreater(laplacian_determinant, 0)
                self.assertLessEqual(laplacian_determinant, plan.grid_laplacian_bound)

    def test_k5_counts_and_certificate(self):
        record = json.loads((HERE / "examples.json").read_text())["examples"][1]
        reduction = SeriesReduction(Graph(record["vertices"], record["original_oriented_edges"]),
                                    record["full_balanced_sources"])
        plan = reciprocal_plan(reduction)
        self.assertEqual(plan.lower_set.M, 683834)
        self.assertEqual(int(plan.support_mask().sum()), 518124)
        self.assertEqual(plan.flat_count, 52)
        self.assertEqual(plan.grid_laplacian_bound, 162000)
        self.assertEqual(reciprocal_gauge_counts(reduction), (683834,) * 10)
        certificate = integer_certificate(reduction)
        self.assertEqual(certificate["moment_denominator"], 216000)
        self.assertEqual(certificate["coefficient_bound"],
                         25384925169104535878250163645356935162028057805002375168000)
        self.assertEqual(certificate["coefficient_bound_bits"], 195)
        self.assertEqual(certificate["crt_modulus_must_exceed"], 2 * certificate["coefficient_bound"])

    def test_validation_and_budgets(self):
        reduction = self.reductions[-1]
        with self.assertRaises(TypeError):
            reciprocal_plan(None)
        with self.assertRaises(ValueError):
            reciprocal_plan(SeriesReduction(polygon(3), [[1, 2], [2, 1]]))
        with self.assertRaises(ValueError):
            reciprocal_plan(reduction, max_subsets=1)
        with self.assertRaises(MemoryError):
            reciprocal_plan(reduction, max_points=1)
        with self.assertRaises(ValueError):
            reciprocal_plan(reduction, gauge=100)
        with self.assertRaises(TypeError):
            reciprocal_plan(reduction, gauge=True)
        with self.assertRaises(ValueError):
            integer_certificate(reduction, max_tree_candidates=1)
        with self.assertRaises(ValueError):
            integer_certificate(reduction, max_permanent_size=1)
        with self.assertRaises(ValueError):
            integer_certificate(reduction, max_moment_terms=1)
        plan = reciprocal_plan(reduction)
        with self.assertRaises(ValueError):
            plan.sample_points(2, 1)
        with self.assertRaises(ValueError):
            plan.support_mask(np.zeros((1, len(reduction.chains)), dtype=float))


def det_direct(core, conductances):
    # Different circuit formulation from reciprocal support enumeration.
    from htk.linalg import det
    return det(brute.laplacian(core, [1 / F(int(g)) for g in conductances]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
