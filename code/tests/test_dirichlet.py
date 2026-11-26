"""Independent rational moment checks for the series-fibre cubature.

Run: python3 code/tests/test_dirichlet.py
The references are elementary Dirichlet factorial moments, not another
cubature implementation.  The ten-interval test covers every independent
coordinate monomial of total degree <= 5, including all mixed moments.
"""

import os
import sys
import unittest
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from htk.dirichlet import (CubatureLimitError, monomial_integral,
                          product_simplex_rule, simplex_rule)


def powers_through(variables, degree):
    for total in range(degree + 1):
        for positions in combinations_with_replacement(range(variables), total):
            powers = [0] * variables
            for position in positions:
                powers[position] += 1
            yield tuple(powers)


def reference_moment(blocks):
    """Independent standard-simplex Lebesgue moment formula."""
    result = F(1)
    for powers in blocks:
        numerator = 1
        for power in powers:
            numerator *= factorial(power)
        result *= F(numerator, factorial(sum(powers) + len(powers) - 1))
    return result


def split_powers(dimensions, powers, all_barycentric):
    offset = 0
    blocks = []
    for dimension in dimensions:
        size = dimension + 1 if all_barycentric else dimension
        block = tuple(powers[offset:offset + size])
        blocks.append(block if all_barycentric else (0,) + block)
        offset += size
    return tuple(blocks)


class DirichletChecks(unittest.TestCase):
    def check_moments(self, dimensions, degree, *, all_barycentric=False):
        rule = product_simplex_rule(dimensions, degree)
        variables = sum(dimensions) + (len(dimensions) if all_barycentric else 0)
        # Reuse exact powers and evaluate only nonconstant factors.
        table = []
        for point, weight in rule.nodes:
            coordinates = tuple(x for block in point
                                for x in (block if all_barycentric else block[1:]))
            table.append((weight, tuple(tuple(x ** power for power in range(degree + 1))
                                        for x in coordinates)))
        count = 0
        for powers in powers_through(variables, degree):
            factors = tuple((axis, power) for axis, power in enumerate(powers) if power)
            actual = F(0)
            for weight, node_powers in table:
                value = weight
                for axis, power in factors:
                    value *= node_powers[axis][power]
                actual += value
            blocks = split_powers(dimensions, powers, all_barycentric)
            expected = reference_moment(blocks)
            self.assertEqual(actual, expected, (dimensions, powers))
            count += 1
        return count

    def test_simplex_moments_at_higher_degree(self):
        self.assertEqual(self.check_moments((1,), 9, all_barycentric=True), 55)
        self.assertEqual(self.check_moments((2,), 7, all_barycentric=True), 120)
        self.check_moments((4,), 5, all_barycentric=True)

    def test_product_moments(self):
        for dimensions in ((1, 1), (1, 2), (2, 2), (1, 1, 1), (0, 1, 2)):
            self.check_moments(dimensions, 5, all_barycentric=True)

    def test_sparse_ten_intervals(self):
        dimensions = (1,) * 10
        rule = product_simplex_rule(dimensions, 5)
        self.assertEqual(len(rule), 221)
        self.assertEqual(self.check_moments(dimensions, 5), 3003)
        self.assertEqual(rule.integrate(lambda point: 1), 1)
        # A high-dimensional single simplex needs fewer points still.
        self.assertEqual(len(simplex_rule(10, 5)), 78)

    def test_interior_points_weights_and_volume(self):
        for dimensions, degree in (((0,), 0), ((1,), 9), ((2,), 7),
                                   ((0, 2, 1, 0), 5), ((1,) * 10, 5)):
            rule = product_simplex_rule(dimensions, degree)
            self.assertEqual(len(rule.nodes), len({point for point, _ in rule.nodes}))
            for point, weight in rule.nodes:
                self.assertIsInstance(weight, F)
                self.assertTrue(weight)
                self.assertEqual(len(point), len(dimensions))
                for block, dimension in zip(point, dimensions):
                    self.assertEqual(len(block), dimension + 1)
                    self.assertTrue(all(isinstance(x, F) and x > 0 for x in block))
                    self.assertEqual(sum(block), 1)
            self.assertEqual(sum(weight for _, weight in rule.nodes), rule.volume)
            self.assertEqual(rule.exact_degree, 2 * (degree // 2) + 1)
        self.assertTrue(any(weight < 0 for _, weight in simplex_rule(2, 3).nodes))

    def test_even_degree_guarantee_and_product_order(self):
        rule = simplex_rule(1, 4)
        self.assertEqual(rule.degree, 4)
        self.assertEqual(rule.exact_degree, 5)
        self.assertEqual(rule.integrate(lambda point: point[0][0] ** 5), F(1, 6))
        forward = product_simplex_rule((1, 2), 5)
        backward = product_simplex_rule((2, 1), 5)
        self.assertEqual({tuple(reversed(point)): weight for point, weight in forward.nodes},
                         dict(backward.nodes))

    def test_empty_and_zero_dimensional_products(self):
        for dimensions in ((), (0,), (0, 0, 0)):
            rule = product_simplex_rule(dimensions, 1000000)
            self.assertEqual(len(rule), 1)
            self.assertEqual(rule.volume, 1)
            self.assertEqual(rule.integrate(lambda point: F(7, 11)), F(7, 11))
        self.assertEqual(monomial_integral(()), 1)
        self.assertEqual(monomial_integral(((19,), (3,))), 1)

    def test_monomial_integral(self):
        for dimensions in ((1,), (2,), (1, 2)):
            variables = sum(dimensions) + len(dimensions)
            for powers in powers_through(variables, 5):
                blocks = split_powers(dimensions, powers, True)
                self.assertEqual(monomial_integral(blocks), reference_moment(blocks))
        self.assertEqual(monomial_integral(((2, 1, 0), (0, 3))), F(1, 240))

    def test_invalid_input_and_budgets(self):
        for dimension in (-1,):
            with self.assertRaises(ValueError):
                simplex_rule(dimension, 3)
        for dimension in (1.5, True, "2"):
            with self.assertRaises(TypeError):
                simplex_rule(dimension, 3)
        with self.assertRaises(ValueError):
            simplex_rule(2, -1)
        for degree in (3.0, False):
            with self.assertRaises(TypeError):
                simplex_rule(2, degree)
        for keyword in ("max_nodes", "max_terms"):
            with self.assertRaises(ValueError):
                simplex_rule(1, 3, **{keyword: 0})
        with self.assertRaises(CubatureLimitError):
            product_simplex_rule((1,) * 10, 5, max_nodes=220)
        with self.assertRaises(CubatureLimitError):
            simplex_rule(2, 7, max_terms=5)
        with self.assertRaises(CubatureLimitError):
            simplex_rule(100, 1, max_terms=10)
        with self.assertRaises(CubatureLimitError):
            product_simplex_rule((0,) * 11, 0, max_terms=10)
        with self.assertRaises(CubatureLimitError):
            simplex_rule(1, 10 ** 12, max_terms=10)
        with self.assertRaises(ValueError):
            monomial_integral(((),))
        with self.assertRaises(ValueError):
            monomial_integral(((1, -1),))
        with self.assertRaises(TypeError):
            monomial_integral(((1, 2.0),))
        with self.assertRaises(TypeError):
            simplex_rule(1, 1).integrate(lambda point: 0.5)

    def test_immutable_cached_rules(self):
        rule = product_simplex_rule((1, 2), 5)
        self.assertIs(rule, product_simplex_rule(iter((1, 2)), 5))
        with self.assertRaises(FrozenInstanceError):
            rule.degree = 4
        with self.assertRaises(TypeError):
            rule.nodes[0] = ()
        # A cache hit at a larger budget must not bypass a smaller budget.
        with self.assertRaises(CubatureLimitError):
            product_simplex_rule((1, 2), 5, max_nodes=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
