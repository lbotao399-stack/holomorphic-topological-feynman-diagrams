"""Complete current-cube controls, independent of global numerator expansion."""

from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from htk import brute
from htk.forest_cubes import (ForestCubes, ForestCubeBudgetExceeded,
                              _Interval, _det_interval)
from htk.graph import Graph
from htk.linalg import det


THREE = ((0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4))
FIVE = ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
        (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))


def full_sources(vertices):
    rows = tuple((F(vertex+1, vertex+2), F((vertex+1)**2+1, vertex+3))
                 for vertex in range(vertices-1))
    return rows + (tuple(-sum(row[a] for row in rows) for a in range(2)),)


def independent_jacobian(cube, point):
    """Differentiate the ORIGINAL grounded Kirchhoff equations, entry by entry."""
    graph = cube.atlas.G
    times = cube.time_map(point)
    currents = brute.currents(graph, times, cube.atlas.sources)
    incidence = graph.incidence()
    result = [[[F(0)]*cube.dimension for _ in range(2)] for _ in graph.edges]
    for direction in range(cube.dimension):
        logarithmic = [F(exponent[direction], point[direction]) for exponent in cube.time_exponents]
        source_derivative = [[sum(F(incidence[edge][vertex])*logarithmic[edge]*currents[edge][a]
                                  for edge in range(graph.ne)) for a in range(2)]
                             for vertex in range(graph.nv-1)]
        response = brute.currents(graph, times, source_derivative)
        for edge in range(graph.ne):
            for a in range(2):
                result[edge][a][direction] = response[edge][a]-logarithmic[edge]*currents[edge][a]
    return tuple(tuple(tuple(row) for row in edge) for edge in result)


class ForestCubeChecks(unittest.TestCase):
    def test_complete_three_loop_atlas_and_original_orientations(self):
        atlas = ForestCubes(Graph(5, THREE), full_sources(5), marked=(3, 4), check=0)
        self.assertEqual((atlas.tree_count, atlas.chart_count), (24, 144))
        self.assertEqual(atlas.stats["ordinary_Hepp_orders"], 5040)
        point = tuple(F(j+2, 2*j+5) for j in range(atlas.dimension))
        for cube in atlas.charts():
            self.assertTrue(all(power >= 0 for row in cube.weight_exponents for power in row))
            self.assertEqual(sum(not any(row) for row in cube.weight_exponents), 1)
            times = cube.time_map(point)
            found, restored = atlas.locate(times)
            self.assertEqual((found.position, restored), (cube.position, point))
            cube.check_original()
        for position in (0, 72, 143):
            cube = atlas.chart(position)
            self.assertEqual(cube.current_jacobian(point), independent_jacobian(cube, point))
            self.assertEqual(len(cube.current_jacobian(point)), len(THREE))
            for corner in ((F(0),)*atlas.dimension, (F(1),)*atlas.dimension):
                self.assertEqual(len(cube.currents(corner)), len(THREE))
                self.assertEqual(len(cube.current_jacobian(corner)), len(THREE))

    def test_five_loop_complete_source_and_current_maps(self):
        atlas = ForestCubes(Graph(7, FIVE), full_sources(7), marked=(5, 6), check=0)
        self.assertEqual((atlas.tree_count, atlas.chart_count), (195, 23400))
        point = tuple(F(j+2, 2*j+5) for j in range(atlas.dimension))
        for position in (0, 11700, 23399):
            cube = atlas.chart(position)
            cube.check_original()
            self.assertEqual(cube.current_jacobian(point), independent_jacobian(cube, point))

    def test_exact_polynomial_series_and_achieved_interval(self):
        graph = Graph(3, ((0, 1), (0, 2), (1, 2)))
        atlas = ForestCubes(graph, ((1, 0), (0, 1), (-1, -1)), marked=(1, 2))
        # On the original simplex y_01,0=t_02 and y_02,1=t_12.
        polynomial = {(0,)*6: F(1), (1, 0, 0, 0, 0, 0): F(1),
                      (1, 0, 0, 1, 0, 0): F(-2)}
        expected = -F(1, 2)-F(1, 6)+2*F(1, 24)
        with patch("htk.forest_periods.CoreCurrentPeriods._source_numerator", side_effect=AssertionError("global numerator used")), \
             patch("htk.forest_periods._research_module", side_effect=AssertionError("period/Hepp backend used")):
            rows = list(atlas.series(polynomial=polynomial).terms(3))
            total = F(0)
            for row in rows:
                total += row.term
                self.assertEqual(total, row.partial_sum)
                self.assertLessEqual(row.lower, expected)
                self.assertGreaterEqual(row.upper, expected)
                self.assertLessEqual(abs(row.partial_sum-expected), row.error_bound)
            for previous, row in zip(rows, rows[1:]):
                self.assertLessEqual(previous.lower, row.lower)
                self.assertGreaterEqual(previous.upper, row.upper)
            result = atlas.series({(0, 0): 2, (2, 1): 1}).enclose(F(1, 2))
        self.assertLessEqual(result.upper-result.lower, F(1, 2))
        self.assertLessEqual(result.lower, -F(1, 15))
        self.assertGreaterEqual(result.upper, -F(1, 15))

    def test_original_edge_permutations_and_reversals(self):
        graph = Graph(3, ((0, 1), (0, 2), (1, 2)))
        source = ((1, 2), (3, -1), (-4, -1))
        original = next(ForestCubes(graph, source).series({(0, 0): 1}).terms(0))
        changed = graph.reverse_edges((0,)).permute_edges((2, 0, 1))
        transformed = next(ForestCubes(changed, source).series({(1, 0): 1}).terms(0))
        # The three-cycle is even; the one original current reversal gives -1.
        self.assertEqual(transformed.partial_sum, -original.partial_sum)
        self.assertEqual((transformed.lower, transformed.upper), (-original.upper, -original.lower))

    def test_singular_centre_interval_converges_without_pivots(self):
        centre = ((F(1), F(1)), (F(1), F(1)))
        previous = None
        for radius in (F(1, 2), F(1, 4), F(1, 8)):
            matrix = tuple(tuple(_Interval(value-radius, value+radius) for value in row) for row in centre)
            bound = _det_interval(matrix, centre, F(0), 1)
            for choices in product((0, 1), repeat=4):
                chosen = [[(matrix[i][j].lower, matrix[i][j].upper)[choices[2*i+j]] for j in range(2)] for i in range(2)]
                value = det(chosen)
                self.assertLessEqual(bound.lower, value)
                self.assertGreaterEqual(bound.upper, value)
            if previous is not None:
                self.assertLess(bound.upper-bound.lower, previous.upper-previous.lower)
            previous = bound

    def test_full_budget_preflights_and_zero_loop(self):
        graph = Graph(5, THREE)
        with patch("htk.forest_cubes.tree_flow", side_effect=AssertionError("allocated tree flows")):
            with self.assertRaises(ForestCubeBudgetExceeded):
                ForestCubes(graph, full_sources(5), max_charts=143)
        atlas = ForestCubes(graph, full_sources(5), max_evaluations=144, check=0)
        with patch.object(atlas, "chart", side_effect=AssertionError("evaluated a partial atlas")):
            with self.assertRaises(ForestCubeBudgetExceeded):
                list(atlas.series().terms(1))
        with self.assertRaises(ForestCubeBudgetExceeded):
            atlas.series({(1, 0): atlas.max_input_degree+1})
        with self.assertRaises(ValueError):
            atlas.series({(1, 0): 1}, polynomial={})
        with self.assertRaises(TypeError):
            ForestCubes(graph, full_sources(5), max_charts=True)
        zero = ForestCubes(Graph(2, ((1, 0),)), ((2, 3), (-2, -3)), marked=(0, 1))
        rows = list(zero.series({(0, 0): 2, (0, 1): 1}).terms(3))
        self.assertEqual([row.term for row in rows], [F(-12), F(0), F(0), F(0)])
        self.assertTrue(all(row.lower == row.upper == -12 for row in rows))
        self.assertEqual(atlas.series(polynomial={}).enclose(F(1, 100)).error_bound, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
