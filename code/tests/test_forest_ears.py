"""Independent sourced-ear checks: exact currents, periods, moments and fallback.

Run: python3 code/tests/test_forest_ears.py
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from flint import fmpq_mpoly_ctx
from htk import Graph
from htk.forest_ears import (EarBundle, ForestEarMoments, EarBudgetExceeded,
                             _source_weights, _fraction)
from htk.forest_periods import (CoreCurrentPeriods, ExactForestSeries, MissingForestPeriod,
                                ForestPeriodBudgetExceeded)
from htk.linalg import perm_sign
from htk.spred import sp_reduce, plan_generic


CORE = ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8))


def fixture(count):
    edges = list(CORE)
    rows = [list(map(F, row)) for row in SOURCES]
    for i in range(count):
        vertex = len(rows)
        source = (F(i + 2), F(3 * i - 3))
        edges.extend(((1, vertex), (vertex, 3)))
        rows.append(list(source))
        for a in range(2):
            rows[4][a] -= source[a]
    return Graph(len(rows), edges), tuple(map(tuple, rows)), (3, 4)


def monomial(graph, word):
    exponents = [0] * (2 * graph.ne)
    for (edge, component), power in word.items():
        exponents[2 * edge + component] = power
    return {tuple(exponents): F(1)}


class ForestEarChecks(unittest.TestCase):
    def test_source_moment_functional_against_independent_box_integrals(self):
        coordinates = ((F(1), F(0)), (F(0), F(1)), (F(2), F(-1)))
        beta = ((1, 2), (0, 1), (2, 0))
        weights = _source_weights(coordinates, beta, 4)
        context = fmpq_mpoly_ctx.get(("r0", "r1", "r2"))
        r0, r1, r2 = context.gens()
        for a in range(5):
            for b in range(5 - a):
                polynomial = (r0 + 2 * r2) ** a * (r1 - r2) ** b
                expected = F(0)
                for exponent, coefficient in polynomial.terms():
                    value = _fraction(coefficient)
                    for power, (alpha, second) in zip(exponent, beta):
                        value *= F(factorial(int(power) + alpha) * factorial(second),
                                   factorial(int(power) + alpha + second + 1))
                    expected += value
                self.assertEqual(sum(weight * point[0] ** a * point[1] ** b for point, weight in weights), expected)

    def test_complete_bundle_maps_and_all_convention_changes(self):
        rng = random.Random(10429)
        graph, sources, marked = fixture(2)
        for changed in (False, True):
            if changed:
                edge_order, labels = list(range(graph.ne)), list(range(graph.nv))
                rng.shuffle(edge_order)
                rng.shuffle(labels)
                changed_graph = graph.reverse_edges(range(0, graph.ne, 2)).permute_edges(edge_order).relabel(labels)
                rows = [None] * graph.nv
                for vertex, row in enumerate(sources):
                    rows[labels[vertex]] = row
                ends, protected = (labels[1], labels[3]), tuple(labels[v] for v in marked)
            else:
                changed_graph, rows, ends, protected = graph, sources, (1, 3), marked
            for endpoints in (ends, ends[::-1]):
                EarBundle(changed_graph, endpoints, marked=protected).check(rows, 2)

    def test_removed_edge_words_against_complete_original_graph_periods(self):
        graph, sources, marked = fixture(1)
        words = ({}, {(7, 0): 1}, {(8, 1): 1}, {(7, 0): 1, (8, 1): 1})
        direct = CoreCurrentPeriods(graph, sources, marked, check=1)
        try:
            with ForestEarMoments(graph, sources, marked, check=1) as ears:
                for word in words:
                    expected = direct.integrate(monomial(graph, word), allow_series=False)
                    self.assertEqual(ears.moment(word, allow_series=False), expected)
                    self.assertEqual(ears.last_stats["reduced_loops"], 3)
                    self.assertTrue(ears.last_stats["all_original_sources_retained"])
                mixed = monomial(graph, {})
                mixed.update({exponent: F(-2, 7) for exponent in monomial(graph, {(7, 0): 1})})
                self.assertEqual(ears.integrate(mixed, allow_series=False),
                                 direct.integrate(mixed, allow_series=False))
        finally:
            direct.close()

    def test_nonparallel_two_source_bundle(self):
        graph, sources, marked = fixture(2)
        # Independently completed on the original five-loop graph, before
        # ear reduction: its numerator has 99,901 monomials.
        expected = {"1": F(-38335280693, 6048000), "z3": F(355647577, 84000)}
        with ForestEarMoments(graph, sources, marked, check=1) as ears:
            self.assertEqual(ears.moment(allow_series=False), expected)
            self.assertEqual(ears.last_stats["source_variables_by_stage"], [2])
            self.assertEqual(ears.last_stats["terminal_source_evaluations"], 21)
            self.assertEqual(ears.last_stats["aggregated_period_terms"], 114)

    def test_high_degree_original_words_and_opposite_parallel_sources(self):
        graph = Graph(4, [(0, 1), (0, 2), (2, 1), (0, 3), (3, 1)])
        word = {(0, 0): 1, (1, 0): 2, (2, 1): 1, (3, 1): 2, (4, 0): 1}
        cases = ((((3, 1), (-4, -2), (2, -3), (-1, 4)), F(-44219303, 1512000)),
                 (((3, 1), (-2, F(-5, 2)), (2, -3), (-3, F(9, 2))), F(39100061, 1152000)))
        order, reversals = (4, 1, 3, 0, 2), {0, 2, 3}
        for sources, expected in cases:
            with ForestEarMoments(graph, sources, (0, 1), check=1) as ears:
                self.assertEqual(ears.moment(word), expected)
            changed = graph.reverse_edges(reversals).permute_edges(order)
            changed_word = {(order.index(edge), component): power for (edge, component), power in word.items()}
            sign = perm_sign(order) * (-1) ** sum(power for (edge, _), power in word.items() if edge in reversals)
            self.assertEqual(sign, -1)
            with ForestEarMoments(changed, sources, (0, 1), check=1) as ears:
                self.assertEqual(ears.moment(changed_word), sign * expected)

    def test_multiple_bundles_against_source_complete_sp_oracle(self):
        graph = Graph(5, [(0, 1), (0, 2), (1, 2), (0, 3), (3, 1), (0, 4), (4, 2)])
        sources = ((1, 2), (3, -1), (-8, -4), (2, -3), (2, 6))
        word = {(3, 0): 1, (6, 1): 1}
        reduction = sp_reduce(graph, D=2, plan=plan_generic(graph))
        expected = reduction.integrate(sources, T=0, powers=word, method="generic", check=1)
        with ForestEarMoments(graph, sources, (1, 2), check=1) as ears:
            self.assertEqual(ears.moment(word), expected)
            self.assertGreater(len(ears.steps), 1)
            self.assertEqual(ears.reduced_graph.loops, 0)
            self.assertEqual(ears.series(word).interval(1), (expected, expected))

    def test_exact_fallback_zero_and_explicit_budgets(self):
        graph, sources, marked = fixture(1)
        with ForestEarMoments(graph, sources, marked, check=0) as ears:
            explicit = ears.series({(8, 0): 1})
            with patch("htk.forest_ears._PeriodWorker.integrate", side_effect=MissingForestPeriod("test unavailable")):
                fallback = ears.moment({(8, 0): 1})
                self.assertIsInstance(fallback, ExactForestSeries)
                self.assertEqual(fallback.edges, explicit.edges)
                self.assertEqual(fallback.components, explicit.components)
                with self.assertRaises(MissingForestPeriod):
                    ears.moment({(8, 0): 1}, allow_series=False)
            self.assertEqual(ears.moment(), {"1": F(-591590731, 483840), "z3": F(661159, 840)})
        with self.assertRaises(EarBudgetExceeded):
            EarBundle(graph, (1, 3), marked=marked, max_source_samples=2).pushforward(sources)
        with self.assertRaises(EarBudgetExceeded):
            EarBundle(graph, (1, 3), marked=marked, max_terms=1).pushforward(sources)
        for option in ("max_terms", "max_source_samples", "max_tree_candidates", "max_period_terms",
                       "max_series_edges", "period_timeout", "max_matrix_entries"):
            for invalid in (0, -1, True, 1.5):
                with self.assertRaises(ValueError):
                    ForestEarMoments(graph, sources, marked, **{option: invalid})
        rows = [list(row) for row in sources]
        for a in range(2):
            rows[4][a] += rows[5][a]
        rows[5] = [F(0), F(0)]
        with ForestEarMoments(graph, rows, marked, check=1, max_source_samples=1) as ears:
            self.assertEqual(ears.moment({(7, 0): 2, (8, 1): 1}), F(0))
            self.assertEqual(ears.last_stats["terminal_source_evaluations"], 0)


if __name__ == "__main__":
    unittest.main()
