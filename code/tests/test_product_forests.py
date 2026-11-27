"""Exact broader-product dispatcher controls; no unsupported source is erased.

Run: python3 code/tests/test_product_forests.py
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk import Graph, ExactPeriod, ExactForestSeries
from htk.forest_model import NotProductForest, forest_graph
from htk.product_forests import (ProductForestMoments, UnsupportedProductForestReduction,
                                 product_forest_moment)
from htk.spred import SPReduction


HIGHER_GRAPH = Graph(5, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (2, 4), (3, 4)])
HIGHER_SOURCES = [[1, 2], [3, -1], [-2, 4], [5, 3], [-7, -8]]


class ProductForestChecks(unittest.TestCase):
    def test_higher_superpotential_full_sources_and_mixed_words(self):
        engine = ProductForestMoments(HIGHER_GRAPH, HIGHER_SOURCES, marked=(3, 4))
        self.assertFalse(engine.structure.cubic)
        self.assertEqual(engine.method, "series_parallel_generic")
        self.assertEqual(engine.sources, tuple(tuple(map(F, row)) for row in HIGHER_SOURCES))
        words = (None, {(0, 0): 1}, {(0, 0): 1, (6, 1): 2})
        density_calls = []
        original = SPReduction.density

        def checked_density(reduction, *args, **kwargs):
            density_calls.append(kwargs)
            return original(reduction, *args, **kwargs)

        with patch.object(SPReduction, "density", checked_density):
            actual = engine.moments(words)
            self.assertEqual(engine.moment(words[-1]), actual[-1])
        self.assertEqual(actual, (F(55625, 432), F(-31083, 80), F(4369801, 8400)))
        self.assertEqual(len(density_calls), 1)
        self.assertEqual(density_calls[0]["method"], "generic")
        self.assertEqual(density_calls[0]["check"], 2)

    def test_cubic_topology_with_extra_source_uses_general_series_reduction(self):
        G, roots = forest_graph((("A", ("A", "B")),))
        # Vertex 1 has full degree three and is sourced, so no source-free star may remove it.
        sources = [[1, 2], [3, -1], [5, 3], [-9, -4]]
        engine = ProductForestMoments(G, sources, marked=(2, 3))
        self.assertTrue(engine.structure.cubic)
        self.assertEqual(engine.method, "series_parallel_generic")
        reference = SPReduction(G).integrate(sources, T=0, powers={(0, 1): 2}, method="generic", check=2)
        self.assertEqual(engine.moment({(0, 1): 2}), reference)
        self.assertEqual(engine.sources[1], (F(3), F(-1)))

    def test_cubic_root_sources_use_local_boundary_messages(self):
        G, roots = forest_graph((("A", "B"),))
        sources = [[1, 2], [5, 3], [-6, -5]]
        engine = ProductForestMoments(G, sources, marked=(1, 2))
        self.assertEqual(engine.method, "cubic_boundary_moments")
        expected = SPReduction(G).integrate(sources, T=0, method="generic", check=2)
        self.assertEqual(engine.moment(), expected)
        self.assertEqual(product_forest_moment(G, sources, marked=(1, 2)), expected)

    def test_previous_non_series_parallel_rejections_have_complete_values(self):
        fork, _ = forest_graph(((("A", "B"), ("A", "B")),))
        sourced_fork = [[1, 2], [3, -1], [0, 0], [5, 3], [-9, -4]]
        with ProductForestMoments(fork, sourced_fork, marked=(3, 4)) as engine:
            self.assertEqual(engine.method, "sourced_core_periods")
            self.assertEqual(engine.moment(), F(1715, 192))
            self.assertEqual(engine.moment({(0, 0): 1}),
                             ExactPeriod({"1": F(-2772371, 43200), "z3": F(2597, 50)}))
            self.assertTrue(engine.last_stats["all_original_sources_retained"])
        noncubic = Graph(7, [(0, 1), (0, 2), (0, 3), (0, 4),
                            (1, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
        sources = [[1, 2], [3, -1], [-2, 4], [0, 0], [0, 0], [5, 3], [-7, -8]]
        with ProductForestMoments(noncubic, sources, marked=(5, 6)) as engine:
            self.assertEqual(engine.method, "sourcefree_contractions")
            self.assertEqual(engine.moment(), 0)
            self.assertEqual(engine.moment({(6, 0): 1}), F(82364363, 2721600))
            self.assertEqual(engine.last_stats["contraction"]["paired_steps"], 1)

    def test_explicit_series_has_exact_terms_and_rigorous_finite_bounds(self):
        fork, _ = forest_graph(((("A", "B"), ("A", "B")),))
        sources = [[1, 2], [3, -1], [0, 0], [5, 3], [-9, -4]]
        with ProductForestMoments(fork, sources, (3, 4)) as engine:
            series = engine.series({(0, 0): 1})
            self.assertIsInstance(series, ExactForestSeries)
            terms = list(series.terms(2))
        self.assertEqual([row.order for row in terms], [0, 1, 2])
        self.assertEqual(terms[-1].partial_sum, sum(row.term for row in terms))
        self.assertTrue(all(isinstance(row.error_bound, F) and row.error_bound >= 0 for row in terms))
        # Integral-test bounds on zeta(3), all rational. They enclose the
        # independently known complete period, not merely a decimal sample.
        zeta_sum = sum((F(1, k**3) for k in range(1, 33)), F(0))
        lower = F(-2772371, 43200) + F(2597, 50) * (zeta_sum + F(1, 2 * 33**2))
        upper = F(-2772371, 43200) + F(2597, 50) * (zeta_sum + F(1, 2 * 32**2))
        for row in terms:
            self.assertLessEqual(row.lower, lower)
            self.assertGreaterEqual(row.upper, upper)

    def test_finite_fast_paths_terminate_without_global_current_expansion(self):
        # Four independent forks have a four-term boundary-message route,
        # whereas generic sequential contraction would exceed this budget.
        count = 4
        A, B = 3 * count, 3 * count + 1
        sources = [[0, 0] for _ in range(3 * count)] + [[1, 0], [-1, -count]]
        edges = []
        for block in range(count):
            root, left, right = 3 * block, 3 * block + 1, 3 * block + 2
            sources[root] = [0, 1]
            edges.extend(((root, left), (root, right), (left, A), (left, B),
                          (right, A), (right, B)))
        edges = sorted(edges + [(A, B)])
        graph = Graph(3 * count + 2, edges)
        powers = {(edges.index((3 * block + 1, A)), 0): 1 for block in range(count)}
        with ProductForestMoments(graph, sources, (A, B), max_terms=2000, check=0) as engine:
            expected = engine.moment(powers)
            self.assertEqual(expected, F(1, 9558692723662258176000))
            self.assertEqual(engine.last_stats["largest_polynomial_terms"], 4)
            with patch("htk.product_forests.ForestContractions", side_effect=AssertionError("fast path lost")):
                series = engine.series(powers)
            self.assertEqual(series.interval(2), (expected, expected))
            self.assertEqual(engine.last_stats["backend"], "exact_terminating_series")
            self.assertEqual(engine.last_stats["finite_backend"], "cubic_boundary_moments")
            self.assertEqual(engine.moment(powers), expected)
            self.assertEqual(engine.last_stats["backend"], "cubic_boundary_moments")
        with ProductForestMoments(HIGHER_GRAPH, HIGHER_SOURCES, (3, 4)) as engine:
            with patch("htk.product_forests.ForestContractions", side_effect=AssertionError("SP fast path lost")):
                series = engine.series()
            self.assertEqual(series.interval(1), (F(55625, 432), F(55625, 432)))
            self.assertEqual(engine.last_stats["finite_backend"], "series_parallel_generic")

    def test_source_formats_and_invalid_data(self):
        expected = F(55625, 432)
        for sources in (HIGHER_SOURCES, HIGHER_SOURCES[:-1], dict(enumerate(HIGHER_SOURCES))):
            self.assertEqual(product_forest_moment(HIGHER_GRAPH, sources, (3, 4)), expected)
        for sources in ([[0, 0]] * 4 + [[1, 0]], [[0]] * 5, [[0.0, 0.0]] * 5):
            with self.assertRaises((TypeError, ValueError)):
                ProductForestMoments(HIGHER_GRAPH, sources, (3, 4))
        with self.assertRaises(NotProductForest):
            ProductForestMoments(Graph(3, [(0, 1), (1, 2)]), [[0, 0]] * 3, (1, 2))
        engine = ProductForestMoments(HIGHER_GRAPH, HIGHER_SOURCES, (3, 4))
        for powers in ({(-1, 0): 1}, {(0, 2): 1}, {(0, 0): -1}, {(0, 0): 1.5}):
            with self.assertRaises((TypeError, ValueError)):
                engine.moment(powers)
        triangle, _ = forest_graph((("A", "B"),))
        triangle_sources = [[1, 2], [5, 3], [-6, -5]]
        for option in ("max_terms", "max_states", "max_tree_candidates", "max_period_terms",
                       "max_series_edges", "period_timeout", "max_matrix_entries"):
            for invalid in (0, -1, True, F(3, 2)):
                for graph, sources, marked in ((triangle, triangle_sources, (1, 2)),
                                                (HIGHER_GRAPH, HIGHER_SOURCES, (3, 4))):
                    with self.assertRaises((TypeError, ValueError)):
                        ProductForestMoments(graph, sources, marked, **{option: invalid})


if __name__ == "__main__":
    unittest.main()
