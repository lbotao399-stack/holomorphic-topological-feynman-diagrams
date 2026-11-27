"""Independent exact controls of the all-forest geometric expansion."""

from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from htk.forest_cube_expansion import (ForestCubeExpansion,
                                     ForestCubeExpansionBudgetExceeded,
                                     _fraction, _q)
from htk.graph import Graph


TRIANGLE = ((0, 1), (0, 2), (1, 2))
THREE = ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4))
SOURCES = ((1, 0), (0, 1), (-1, -1))


class GeometricCubeChecks(unittest.TestCase):
    def test_finite_tail_identity_and_monotonicity(self):
        # Independent finite subtraction from the complete binomial sum.
        for trees in range(1, 10):
            centre = F(trees+1, 2)
            ratio = F(trees-1, trees+1)
            for power in range(1, 9):
                previous = F(1)
                for order in range(12):
                    subtraction = 1-sum((F(comb(power+k-1, k))*ratio**k/centre**power
                                          for k in range(order+1)), F(0))
                    positive = ratio**(order+1)*sum((F(comb(order+j, j))/centre**j
                                                     for j in range(power)), F(0))
                    self.assertEqual(positive, subtraction)
                    self.assertGreaterEqual(positive, 0)
                    self.assertLessEqual(positive, previous)
                    previous = positive

    def test_triangle_every_term_and_exact_infinite_sum(self):
        series = ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES, marked=(1, 2)).series()
        rows = list(series.terms(12))
        for row in rows:
            expected_term = F(0) if row.order % 2 else -F(3, 8*4**(row.order//2))
            self.assertEqual(row.term, expected_term)
            self.assertEqual(abs(row.partial_sum+F(1, 2)), F(1, 2*4**(row.order//2+1)))
            self.assertLessEqual(row.lower, -F(1, 2))
            self.assertGreaterEqual(row.upper, -F(1, 2))
        # -3/8 sum(j>=0) 4**(-j) = -1/2, independently of period tables.
        self.assertEqual(-F(3, 8)/(1-F(1, 4)), -F(1, 2))
        self.assertEqual(series.required_order(F(1, 1000)), 18)
        result = series.enclose(F(1, 1000))
        self.assertEqual(result.partial_sum, -F(1048575, 2097152))
        self.assertEqual(result.error_bound, F(87, 262144))
        self.assertLessEqual(result.upper-result.lower, F(1, 1000))

    def test_mixed_current_polynomial_and_complete_spinor(self):
        engine = ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES, marked=(1, 2))
        polynomial = {(0,)*6: F(1), (1, 0, 0, 0, 0, 0): F(1),
                      (1, 0, 0, 1, 0, 0): -F(2)}
        # Original simplex moments: -1/2-1/6+2/24=-7/12.
        exact = -F(7, 12)
        with patch("htk.forest_periods.CoreCurrentPeriods._source_numerator",
                   side_effect=AssertionError("global original numerator used")), \
             patch("htk.forest_periods._research_module",
                   side_effect=AssertionError("period backend used")):
            result = engine.series(polynomial=polynomial).enclose(F(1, 1000))
            self.assertLessEqual(result.lower, exact)
            self.assertGreaterEqual(result.upper, exact)
            self.assertLessEqual(abs(result.partial_sum-exact), result.error_bound)
            # y_01=(t_02,-t_12) and simplex density=-1.
            for entry, expected in zip(engine.spinor(0), (-F(1, 6), F(1, 6))):
                value = entry.enclose(F(1, 1000))
                self.assertLessEqual(value.lower, expected)
                self.assertGreaterEqual(value.upper, expected)

    def test_augmented_determinant_retains_all_three_loop_sources(self):
        source = ((F(2, 3), -F(1, 7)), (F(1, 2), F(3, 4)),
                  (-F(1, 5), F(2, 9)), (F(7, 11), F(4, 13)))
        engine = ForestCubeExpansion(Graph(5, THREE), source, marked=(3, 4), check=0)
        for index in (0, 1, 36, 72, 142, 143):
            cube = engine.atlas.chart(index)
            denominator, currents, numerator, power = engine._chart(cube)
            for point in (tuple(F(j+2, 2*j+5) for j in range(6)), (F(1),)*6, (F(0),)*6):
                values = tuple(map(_q, point))
                divisor = _fraction(denominator(*values))
                self.assertGreaterEqual(divisor, 1)
                self.assertLessEqual(divisor, engine.atlas.tree_count)
                expected = tuple(value for row in cube.currents(point) for value in row)
                self.assertEqual(tuple(_fraction(current(*values))/divisor for current in currents), expected)
                self.assertEqual(_fraction(numerator(*values))/divisor**power, cube.density(point))
            self.assertLessEqual(power, 7)

    def test_original_edge_reversal_and_permutation_at_every_order(self):
        graph = Graph(3, TRIANGLE)
        source = ((1, 2), (3, -1), (-4, -1))
        original = list(ForestCubeExpansion(graph, source).series({(0, 0): 1}).terms(6))
        changed = graph.reverse_edges((0,)).permute_edges((2, 0, 1))
        transformed = list(ForestCubeExpansion(changed, source).series({(1, 0): 1}).terms(6))
        for a, b in zip(original, transformed):
            self.assertEqual(a.term, -b.term)
            self.assertEqual(a.partial_sum, -b.partial_sum)
            # The complete moments have the exact sign.  Coefficientwise
            # majorants can differ after chart-variable relabeling and
            # cancellation, but both must enclose the independent value.
            self.assertLessEqual(a.lower, -F(7, 3))
            self.assertGreaterEqual(a.upper, -F(7, 3))
            self.assertLessEqual(b.lower, F(7, 3))
            self.assertGreaterEqual(b.upper, F(7, 3))

    def test_zero_loop_zero_input_and_replayable_complete_specification(self):
        engine = ForestCubeExpansion(Graph(2, ((1, 0),)), ((2, 3), (-2, -3)), marked=(0, 1))
        series = engine.series({(0, 0): 2, (0, 1): 1})
        rows = list(series.terms(4))
        self.assertEqual([row.term for row in rows], [-F(12), F(0), F(0), F(0), F(0)])
        self.assertTrue(all(row.lower == row.upper == -12 for row in rows))
        self.assertEqual(engine.series(polynomial={}).enclose(F(1, 1000)).error_bound, 0)
        specification = ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES).series().specification()
        self.assertEqual(specification["absolute_ratio"], "1/2")
        self.assertEqual(specification["complete_cubes"], 3)
        self.assertEqual(specification["components"], [{
            "power": 3, "denominator": [[[0, 0], "1"], [[0, 1], "1"], [[1, 0], "1"]],
            "numerator": [[[0, 0], "-3"]], "absolute_moment": "3"}])

    def test_explicit_budgets_and_atomic_retry(self):
        small = ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES, max_terms=1)
        with self.assertRaises(ForestCubeExpansionBudgetExceeded):
            small.series().specification()
        engine = ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES, max_order=2)
        with self.assertRaises(ForestCubeExpansionBudgetExceeded):
            engine.series().required_order(F(1, 1000))
        with self.assertRaises(ForestCubeExpansionBudgetExceeded):
            list(engine.series().terms(3))
        series = engine.series()
        original = engine._chart
        calls = 0
        def fail_after_polynomial_component(cube):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise ForestCubeExpansionBudgetExceeded("deliberate complete-atlas interruption")
            denominator, currents, unused, unused_power = original(cube)
            return denominator, currents, engine.context.constant(7), 0
        with patch.object(engine, "_chart", side_effect=fail_after_polynomial_component):
            with self.assertRaises(ForestCubeExpansionBudgetExceeded):
                series.specification()
        self.assertIsNone(series.components)
        self.assertEqual(series.constant, 0)
        self.assertEqual(series.specification()["constant"], "0")
        self.assertEqual(next(series.terms(0)).term, -F(3, 8))
        with self.assertRaises(TypeError):
            ForestCubeExpansion(Graph(3, TRIANGLE), SOURCES, max_order=True)
        with self.assertRaises(TypeError):
            engine.series().enclose(0.1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
