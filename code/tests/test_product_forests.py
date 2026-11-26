"""Exact broader-product dispatcher controls; no unsupported source is erased.

Run: python3 code/tests/test_product_forests.py
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk import Graph
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

    def test_sourced_non_series_parallel_cases_raise_without_returning_zero(self):
        fork, _ = forest_graph(((("A", "B"), ("A", "B")),))
        sourced_fork = [[1, 2], [3, -1], [0, 0], [5, 3], [-9, -4]]
        with self.assertRaises(UnsupportedProductForestReduction):
            ProductForestMoments(fork, sourced_fork, marked=(3, 4))
        noncubic = Graph(7, [(0, 1), (0, 2), (0, 3), (0, 4),
                            (1, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
        sources = [[1, 2], [3, -1], [-2, 4], [0, 0], [0, 0], [5, 3], [-7, -8]]
        with self.assertRaises(UnsupportedProductForestReduction):
            ProductForestMoments(noncubic, sources, marked=(5, 6))

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


if __name__ == "__main__":
    unittest.main()
