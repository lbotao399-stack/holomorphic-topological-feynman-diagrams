"""Independent exact checks of sourced-core finite periods and rational series."""
from fractions import Fraction as F
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk.graph import Graph
from htk.forest_periods import (CoreCurrentPeriods, ExactPeriod, ExactForestSeries,
                                ForestPeriodBudgetExceeded, MissingForestPeriod)


FORK = Graph(5, [(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)])
SOURCES = ((1, 2), (3, -1), (0, 0), (5, 3), (-9, -4))


def word(graph, *locations):
    exponent = [0] * (2 * graph.ne)
    for edge, component in locations:
        exponent[2 * edge + component] += 1
    return tuple(exponent)


class ForestPeriodTests(unittest.TestCase):
    def setUp(self):
        self.reducer = CoreCurrentPeriods(FORK, SOURCES, check=2)

    def tearDown(self):
        self.reducer.close()

    def test_independent_fork_values_and_aggregation(self):
        # Independently obtained by exact Newton interpolation from the
        # original grounded Kirchhoff density, followed by wheel residues.
        self.assertEqual(self.reducer.integrate({word(FORK): F(1)}), F(1715, 192))
        first = self.reducer.integrate({word(FORK, (0, 0)): F(1)})
        self.assertEqual(first, {"1": -F(2772371, 43200), "z3": F(2597, 50)})
        mixed = self.reducer.integrate({word(FORK, (0, 0), (6, 1)): F(1)})
        self.assertEqual(mixed, {"1": -F(9359, 120), "z3": F(637, 10)})
        other = self.reducer.integrate({word(FORK, (2, 0), (3, 1)): F(1)})
        self.assertEqual(other, {"1": -F(526309, 21600), "z3": F(11123, 600)})
        numerator = self.reducer._bare
        total = self.reducer.integrate({word(FORK, (0, 0), (6, 1)): F(3, 7),
                                        word(FORK, (2, 0), (3, 1)): -F(2, 9)})
        self.assertEqual(total, F(3, 7) * mixed - F(2, 9) * other)
        self.assertIs(self.reducer._bare, numerator)
        self.assertLess(self.reducer.last_stats["core_period_orbits"],
                        self.reducer.last_stats["inserted_numerator_terms"])

    def test_rational_sources_and_original_orientation(self):
        scaled = [[F(row[0], 2), F(row[1], 3)] for row in SOURCES]
        # H=2 bare homogeneity is n in each component; y^0 adds one
        # first-component degree. Reversing its edge changes that insertion.
        graph = FORK.reverse_edges([0, 3])
        reducer = CoreCurrentPeriods(graph, scaled, check=2)
        try:
            value = reducer.integrate({word(graph, (0, 0)): F(1)})
            expected = ExactPeriod({"1": -F(2772371, 43200), "z3": F(2597, 50)})
            self.assertEqual(value, -F(1, 2 ** 4 * 3 ** 3) * expected)
        finally:
            reducer.close()

    def test_current_bracket_cancels_U_before_periods(self):
        polynomial = {word(FORK, (2, 0), (3, 1)): F(1),
                      word(FORK, (2, 1), (3, 0)): -F(1)}
        value = self.reducer.integrate(polynomial)
        self.assertEqual(value, {"1": -F(84721, 7200), "z3": F(686, 25)})
        self.assertEqual(self.reducer.last_stats["current_U_cancellations"], {2: 1})
        self.assertEqual(self.reducer.last_stats["denominator_powers"], [6])
        # Verify the reduced denominator pointwise against the original
        # Kirchhoff currents, independently of period identities.
        from htk import brute
        from flint import fmpq
        numerator = self.reducer._schwinger(polynomial)[6]
        point = [F(i + 2, i + 1) for i in range(FORK.ne)]
        current = brute.currents(FORK, point, SOURCES)
        expected = brute.density(FORK, point, SOURCES).rational() * (
            current[2][0] * current[3][1] - current[2][1] * current[3][0])
        actual = F(str(numerator(*(fmpq(value.numerator, value.denominator) for value in point))))
        actual /= brute.symanzik_U(FORK, point) ** 6
        self.assertEqual(actual, expected)

    def test_mixed_degrees_and_kirchhoff_zero(self):
        value = self.reducer.integrate({word(FORK): F(3, 2),
                                        word(FORK, (0, 0)): F(-2, 3)})
        expected = F(3, 2) * F(1715, 192) + F(-2, 3) * ExactPeriod(
            {"1": -F(2772371, 43200), "z3": F(2597, 50)})
        self.assertEqual(value, expected)
        # The two original edges at vertex zero point outwards, and its
        # first source component is one. Cancelling their exact U factors
        # must combine with the constant part before period evaluation.
        zero = self.reducer.integrate({word(FORK, (0, 0)): F(1),
                                       word(FORK, (1, 0)): F(1), word(FORK): F(-1)})
        self.assertEqual(zero, F(0))
        self.assertEqual(self.reducer.last_stats["backend"], "exact_zero")

    def test_series_parallel_polynomial_sum(self):
        graph = Graph(3, [(0, 1), (0, 2), (1, 2)])
        reducer = CoreCurrentPeriods(graph, [(1, 2), (3, -1), (-4, -1)])
        first = reducer.integrate({word(graph, (0, 0)): F(1)})
        second = reducer.integrate({word(graph, (1, 1), (2, 0)): F(1)})
        value = reducer.integrate({word(graph, (0, 0)): F(2),
                                  word(graph, (1, 1), (2, 0)): -F(3, 5)})
        self.assertEqual(value, 2 * first - F(3, 5) * second)
        self.assertEqual(reducer.last_stats["backend"], "series_parallel_generic")
        series = reducer.series({word(graph): F(1)})
        self.assertEqual(series.interval(1), (F(7, 2), F(7, 2)))

    def test_private_imports_and_finite_arithmetic(self):
        sentinel = SimpleNamespace(sentinel=True)
        path = list(sys.path)
        with patch.dict(sys.modules, {"core": sentinel, "kernel": sentinel}):
            value = self.reducer.integrate({word(FORK, (0, 0)): F(1)})
            self.assertIs(sys.modules["core"], sentinel)
            self.assertIs(sys.modules["kernel"], sentinel)
        self.assertEqual(sys.path, path)
        self.assertEqual(value + (-value), F(0))
        self.assertEqual((value + F(2))["1"], -F(2772371, 43200) + 2)
        self.assertEqual(value.as_sympy().free_symbols, set())

    def test_exact_series_and_independent_rational_period(self):
        # On a triangle U=t0+t1+t2; fixing U=1 leaves a simplex of area 1/2.
        series = ExactForestSeries([(0, 1), (0, 2), (1, 2)], [(3, {(0, 0, 0): F(7)})])
        rows = list(series.terms(20))
        self.assertEqual([row.order for row in rows], list(range(21)))
        for row in rows:
            self.assertLessEqual(row.lower, F(7, 2))
            self.assertLessEqual(F(7, 2), row.upper)
            self.assertIsInstance(row.term, F)
            self.assertIsInstance(row.error_bound, F)
        self.assertLess(rows[-1].error_bound, rows[0].error_bound)
        self.assertEqual(sum(row.term for row in rows), rows[-1].partial_sum)
        # The alternative exact binomial expansion has an independent
        # recurrence for its rational coefficients and the same integral.
        lower, upper = series.interval(20, method="binomial")
        self.assertLessEqual(lower, F(7, 2))
        self.assertLessEqual(F(7, 2), upper)

    def test_missing_data_is_an_exact_series_or_explicit_error(self):
        with patch("htk.forest_periods._PeriodWorker.integrate", side_effect=MissingForestPeriod("test missing table")):
            result = self.reducer.integrate({word(FORK): F(1)})
            self.assertIsInstance(result, ExactForestSeries)
            row = next(result.terms(0))
            self.assertLessEqual(row.lower, F(1715, 192))
            self.assertLessEqual(F(1715, 192), row.upper)
            self.assertEqual(self.reducer.last_stats["backend"], "exact_convergent_series")
            with self.assertRaises(MissingForestPeriod):
                self.reducer.integrate({word(FORK): F(1)}, allow_series=False)

    def test_worker_timeout_allows_a_fresh_retry(self):
        import selectors
        with patch.object(selectors.DefaultSelector, "select", return_value=[]):
            with self.assertRaises(ForestPeriodBudgetExceeded):
                self.reducer.integrate({word(FORK): F(1)})
        self.assertIsNone(self.reducer._worker)
        self.assertEqual(self.reducer.integrate({word(FORK): F(1)}), F(1715, 192))
        self.reducer._worker.process.kill()
        self.reducer._worker.process.wait()
        with self.assertRaises(ArithmeticError):
            self.reducer.integrate({word(FORK): F(1)})
        self.assertIsNone(self.reducer._worker)
        self.assertEqual(self.reducer.integrate({word(FORK): F(1)}), F(1715, 192))

    def test_budgets_and_invalid_inputs_are_not_zero(self):
        reducer = CoreCurrentPeriods(FORK, SOURCES, max_matrix_entries=1)
        with self.assertRaises(ForestPeriodBudgetExceeded):
            reducer.integrate({word(FORK): F(1)})
        with self.assertRaises(TypeError):
            self.reducer.integrate({word(FORK): 0.5})
        with self.assertRaises(ValueError):
            self.reducer.integrate({(0, 0): F(1)})
        series = ExactForestSeries([(0, 1), (0, 2), (1, 2)], [(3, {(0, 0, 0): F(1)})], max_edges=2)
        with self.assertRaises(ForestPeriodBudgetExceeded):
            next(series.terms(0))
        series = ExactForestSeries([(0, 1), (0, 2), (1, 2)], [(3, {(0, 0, 0): F(1)})], max_terms=1)
        with self.assertRaises(ForestPeriodBudgetExceeded):
            list(series.terms(1))


if __name__ == "__main__":
    unittest.main()
