"""Exact controls for the general bulk dispatcher and half-integer series."""
from fractions import Fraction
from pathlib import Path
from math import prod
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from htk import brute
from htk.graph import Graph, banana, polygon, theta
from htk.topological_polynomial import BulkPolynomial, BulkPolynomialBudgetExceeded, rational
from htk.bulk_kernels import evaluate_bulk, _check_admissible
from htk.bulk_series import ExactBulkSeries, BulkSeriesBudgetExceeded
from htk.forest_cubes import GreedyCube
from htk.forest_cube_expansion import _q
from htk.spred import sp_reduce


def pi_bounds(terms=24):
    # Alternating-series bounds, combined with the exact Machin identity.
    def arctangent_bounds(denominator):
        value = sum((Fraction((-1)**k, (2*k+1)*denominator**(2*k+1))
                     for k in range(terms)), Fraction(0))
        other = value + Fraction((-1)**terms, (2*terms+1)*denominator**(2*terms+1))
        return min(value, other), max(value, other)
    lo5, hi5 = arctangent_bounds(5)
    lo239, hi239 = arctangent_bounds(239)
    return 16*lo5-4*hi239, 16*hi5-4*lo239


class BulkKernelChecks(unittest.TestCase):
    def test_half_integer_series_encloses_pi_with_a_narrow_rational_interval(self):
        density = BulkPolynomial(banana(3), [[]], H=0, T=1)
        series = ExactBulkSeries(density)
        lower, upper = series.interval(60)
        known_lower, known_upper = pi_bounds()
        self.assertLessEqual(lower, known_lower)
        self.assertGreaterEqual(upper, known_upper)
        self.assertLess(upper-lower, Fraction(1, 10000))
        self.assertEqual(series.statistics["complete_charts"], 6)
        self.assertTrue(series.statistics["squared_coordinates"])

    def test_half_integer_series_all_original_currents(self):
        graph, _ = theta(2, 2, 1)
        sources = [[1], [2], [-1]]
        density = BulkPolynomial(graph, sources, H=1, T=1)
        word = ((0, 0), (3, 0))
        series = ExactBulkSeries(density, word)
        atlas = series._atlas()
        # A pullback of the complete inserted density is checked against
        # the original grounded Laplacian in every dominant-tree chart.
        from itertools import permutations
        point = tuple(Fraction(i+1, i+3) for i in range(graph.ne-1))
        for tree_index, tree in enumerate(atlas.trees):
            for order in permutations(graph.cotree(tree)):
                chart = GreedyCube(atlas, 0, tree_index, order)
                denominator, numerator = series._chart_polynomial(chart)
                denominator = rational(denominator(*map(_q, point)))
                numerator = rational(numerator(*map(_q, point)))
                times = chart.time_map(tuple(x*x for x in point))
                jacobian = 2**series.dimension*prod(times)/prod(point)
                original = brute.density(graph, times, sources, T=1)
                currents = brute.currents(graph, times, sources)
                original = original.scaled(jacobian*prod(currents[e][a] for e, a in word))
                represented = brute.Density(numerator/denominator**int(series.power-Fraction(1,2)), denominator)
                self.assertTrue(original.same_as(represented))

    def test_quadratic_gaussian_identity_in_three_holomorphic_dimensions(self):
        graph, _ = theta(3, 3, 3)
        sources = [[(v+1)**(a+1) for a in range(3)] for v in range(graph.nv-1)]
        density = BulkPolynomial(graph, sources, H=3, T=1)
        self.assertEqual(density.method, "gaussian_minor")
        self.assertTrue(density.numerator)
        # Independent grounded-Laplacian Pfaffians use every source
        # component, without the polynomial minor/polarization algorithm.
        for times in (tuple(range(1, graph.ne+1)),
                      tuple(Fraction(e+3, e+1) for e in range(graph.ne))):
            self.assertTrue(density.check_at(times))
        # The original full series-parallel fibre integration supplies
        # an independent integral, including a mixed-component insertion.
        word = ((0, 1), (4, 2))
        actual = evaluate_bulk(graph, sources, H=3, T=1, word=word, route="periods")
        reference = sp_reduce(graph, D=4).integrate(sources, T=1,
                                                  powers={word[0]: 1, word[1]: 1})
        self.assertEqual(actual.as_sympy(), reference.sympy())

    def test_dispatcher_keeps_every_spinor_direction(self):
        graph = polygon(4)
        sources = [[1, 2, 3], [2, -1, 1], [-1, 3, 2]]
        for component in range(3):
            word = ((0, component), (2, (component+1)%3))
            actual = evaluate_bulk(graph, sources, H=3, T=0, word=word)
            reference = sp_reduce(graph, D=3).integrate(sources, powers={word[0]: 1, word[1]: 1})
            self.assertEqual(dict(actual.coefficients), {"1": reference} if reference else {})

    def test_new_complete_k33_derivative_is_nonzero(self):
        graph = Graph(6, [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)])
        sources = [[1], [-1], [0], [1], [0]]
        bare = evaluate_bulk(graph, sources, H=1, T=1)
        current = evaluate_bulk(graph, sources, H=1, T=1, word=((0,0),))
        self.assertEqual(dict(bare.coefficients), {})
        self.assertEqual(dict(current.coefficients), {"sqrtpi^4": Fraction(2243,4725),
                                                       "sqrtpi^8": Fraction(-121,2520)})

    def test_current_orientation_and_edge_order(self):
        edges = [(0,3),(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]
        graph = Graph(6, edges)
        sources = [[1], [-1], [0], [1], [0]]
        original = evaluate_bulk(graph, sources, H=1, T=1, word=((0,0),))
        reversed_graph = graph.reverse_edges([0])
        reversed_value = evaluate_bulk(reversed_graph, sources, H=1, T=1, word=((0,0),))
        self.assertEqual(dict(reversed_value.coefficients), {k: -v for k, v in original.coefficients.items()})
        permutation = [1,0,2,3,4,5,6,7,8]
        permuted_value = evaluate_bulk(graph.permute_edges(permutation), sources, H=1, T=1, word=((1,0),))
        self.assertEqual(dict(permuted_value.coefficients), {k: -v for k, v in original.coefficients.items()})

    def test_bulk_vanishing_retains_insertions(self):
        graph = polygon(3)
        self.assertEqual(dict(evaluate_bulk(graph, [[1], [2]], H=1, T=1, word=((0,0),)).coefficients), {})
        self.assertEqual(dict(evaluate_bulk(graph, None, H=0, T=2).coefficients), {})

    def test_exact_inputs_and_admissibility_fail_closed(self):
        graph = polygon(3)
        with self.assertRaises(TypeError):
            evaluate_bulk(graph, [[1.0], [2]], H=1, T=1)
        with self.assertRaises(TypeError):
            BulkPolynomial(graph, [[1.0], [2]], H=1, T=1)
        with self.assertRaises(ValueError):
            evaluate_bulk(graph, [[1], [2]], H=1, T=1, word=((0,1),))
        bad = Graph(5, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(3,4)])
        with self.assertRaises(ValueError):
            _check_admissible(bad, 2)
        for invalid in ({0: [1,0], 1: [0,1], 99: [3,4]},
                        [[1,0], [0,1], [-1,-1,7]],
                        {0: [1,0], 1: [0,1], 2: [-1,-1,7]},
                        {True: [1,0]}):
            for constructor in (evaluate_bulk, BulkPolynomial):
                with self.assertRaises((ValueError, TypeError)):
                    constructor(graph, invalid, H=2, T=0)
        sparse = evaluate_bulk(graph, {0: [1,2]}, H=2, T=0)
        complete = evaluate_bulk(graph, [[1,2], [0,0], [-1,-2]], H=2, T=0)
        self.assertEqual(dict(sparse.coefficients), dict(complete.coefficients))

    def test_incomplete_cube_atlas_is_never_returned(self):
        density = BulkPolynomial(banana(3), [[]], H=0, T=1)
        with self.assertRaises(BulkSeriesBudgetExceeded):
            ExactBulkSeries(density, max_charts=1).interval(0)
        with self.assertRaises(BulkSeriesBudgetExceeded):
            ExactBulkSeries(density, max_total_terms=1).interval(1)
        with self.assertRaises(BulkPolynomialBudgetExceeded):
            BulkPolynomial(banana(3), [[]], H=0, T=1, max_terms=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
