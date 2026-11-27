"""Complete polynomial composition, source transport and exact fallback controls.

Run: python3 code/tests/test_forest_reduction.py
"""
from fractions import Fraction as F
from itertools import product
from math import prod
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk import Graph, brute
from htk.forest_model import forest_graph
from htk.forest_contractions import ForestContractions
from htk.forest_ears import ForestEarMoments, EarBudgetExceeded
from htk.forest_periods import ExactPeriod, ExactForestSeries, MissingForestPeriod
from htk.forest_reduction import (ForestReduction, ForestReductionBudgetExceeded,
                                  reduce_product_forest)
from htk.linalg import perm_sign
from htk.spred import SPReduction
import htk.forest_reduction as implementation


HIGHER = Graph(5, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (2, 4), (3, 4)])
HIGHER_SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8))
FIVE = Graph(7, [(0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
                 (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
FIVE_SOURCES = ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (0, 0), (-1, -1))
FIVE_VALUE = ExactPeriod({'1': F(1187, 435456000), 'z3': F(-17, 9072000)})


def polynomial(graph, terms):
    result = {}
    for coefficient, word in terms:
        exponent = [0] * (2 * graph.ne)
        for (edge, component), power in word.items():
            exponent[2 * edge + component] = power
        exponent = tuple(exponent)
        result[exponent] = result.get(exponent, F(0)) + coefficient
    return {exponent: F(value) for exponent, value in result.items() if value}


def evaluate(polynomial, currents):
    currents = tuple(value for row in currents for value in row)
    return sum(coefficient * prod(value ** power for value, power in zip(currents, exponent))
               for exponent, coefficient in polynomial.items())


def passive_ear_fixture():
    # A complete sourced ear on a non-SP three-loop fork, with one
    # inverse zero-source single move on that ear's triangle.
    edges = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
             (5, 3), (6, 1), (6, 5), (6, 3)]
    sources = ((1, 2), (3, -1), (-2, 4), (5, 3), (-9, -5), (2, -3), (0, 0))
    return Graph(7, edges), sources, (3, 4)


class ForestReductionChecks(unittest.TestCase):
    def test_cubic_mixed_polynomial_and_complete_tensor_against_generic_density(self):
        graph, _ = forest_graph((("A", ("A", "B")),))
        sources = ((1, 2), (0, 0), (5, 3), (-6, -5))
        insertion = polynomial(graph, [(1, {}), (F(-2, 7), {(0, 0): 1}),
                                        (F(3, 5), {(1, 1): 1, (2, 0): 2})])
        direct = SPReduction(graph)
        expected = direct.integrate(sources, T=0, method='generic', check=1,
                                    insertion=lambda currents: evaluate(insertion, currents))
        with ForestReduction(graph, sources, (2, 3), check=1) as engine:
            with patch.object(implementation, 'CoreCurrentPeriods', side_effect=AssertionError('lost cubic route')):
                self.assertEqual(engine.integrate(insertion), expected)
                tensor = engine.tensor((0, 1), indices=('a', 'b'))
            self.assertEqual(tensor.indices, ('a', 'b'))
            self.assertEqual(tensor.edge_slots, (0, 1))
            self.assertEqual(tensor.rank, 2)
            self.assertEqual(set(tensor.entries), set(product((0, 1), repeat=2)))
            for coordinates, value in tensor.entries.items():
                word = {(0, coordinates[0]): 1, (1, coordinates[1]): 1}
                self.assertEqual(value, direct.integrate(sources, powers=word, method='generic', check=1))
            self.assertEqual(engine.series(insertion).interval(1), (expected, expected))
            with self.assertRaises(TypeError):
                tensor.entries[(0, 0)] = F(0)

    def test_generic_higher_valence_polynomial_reuses_one_full_density(self):
        insertion = polynomial(HIGHER, [(1, {}), (F(-2, 7), {(0, 0): 1}),
                                        (F(5, 11), {(0, 0): 1, (6, 1): 2})])
        expected = F(55625, 432) - F(2, 7) * F(-31083, 80) + F(5, 11) * F(4369801, 8400)
        calls, original = [], SPReduction.density
        def density(reduction, *args, **kwargs):
            calls.append(kwargs)
            return original(reduction, *args, **kwargs)
        with ForestReduction(HIGHER, HIGHER_SOURCES, (3, 4), check=1) as engine:
            with patch.object(SPReduction, 'density', density):
                self.assertEqual(engine.integrate(insertion), expected)
                self.assertEqual(engine.moment(), F(55625, 432))
                self.assertEqual(engine.moment({(0, 0): 1}), F(-31083, 80))
                self.assertEqual(engine.integrate(insertion), expected)
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0]['method'], 'generic')
            self.assertEqual(engine.sources, tuple(tuple(map(F, row)) for row in HIGHER_SOURCES))

    def test_passive_star_then_sourced_ear_recombines_every_degree(self):
        graph, sources, marked = passive_ear_fixture()
        insertion = polynomial(graph, [(1, {}), (F(-2, 7), {(8, 0): 1})])
        contractions = ForestContractions(graph, sources, marked, check=1)
        combined = {}
        for coefficient, word in ((1, {}), (F(-2, 7), {(8, 0): 1})):
            result = contractions.contract(polynomial=polynomial(graph, [(coefficient, word)]))
            for exponent, value in result.polynomial.items():
                combined[exponent] = combined.get(exponent, F(0)) + value
        with ForestEarMoments(contractions.reduced_graph, contractions.reduced_sources,
                              contractions.reduced_marked, check=1) as ears:
            expected = ears.integrate(combined, allow_series=False)
        self.assertEqual(expected, ExactPeriod({'1': F(-49356274469, 60963840),
                                               'z3': F(663612401, 1323000)}))
        with ForestReduction(graph, sources, marked, check=1) as engine:
            self.assertEqual(engine.integrate(insertion, allow_series=False), expected)
            self.assertEqual(engine.method, 'sourcefree_stars_then_sourced_ears')
            self.assertEqual(len(engine.last_stats['contraction']['homogeneous_groups']), 2)
            self.assertEqual(engine.last_stats['reduced_loops'], 3)
            self.assertEqual(engine.last_stats['ear_bundles'], 1)
            self.assertTrue(engine.last_stats['all_original_sources_retained'])

    def test_kirchhoff_complete_cancellation_and_constant_before_period_lookup(self):
        # At vertex1, y_(1,2),0 - y_(0,1),0 = source_(1,0) = 1.
        # The complete square plus two linear terms and three equals five.
        insertion = polynomial(FIVE, [(1, {(1, 0): 2}), (1, {(0, 0): 2}),
                                      (-2, {(1, 0): 1, (0, 0): 1}), (4, {})])
        zero = polynomial(FIVE, [(1, {(1, 0): 1}), (-1, {(0, 0): 1}), (-1, {})])
        with ForestReduction(FIVE, FIVE_SOURCES, (5, 6), check=1) as engine:
            normal = engine.normal_form(insertion)
            self.assertEqual(normal, {(0,) * (2 * FIVE.ne): F(5)})
            for shift in (1, 3):
                times = tuple(F(edge + shift, edge + shift + 2) for edge in range(FIVE.ne))
                currents = brute.currents(FIVE, times, engine.sources)
                self.assertEqual(evaluate(insertion, currents), evaluate(normal, currents))
            with patch.object(implementation, 'CoreCurrentPeriods', side_effect=AssertionError('needless core')):
                self.assertEqual(engine.integrate(insertion), 5 * FIVE_VALUE)
                self.assertEqual(engine.integrate(zero), 0)
                self.assertEqual(engine.moment(), FIVE_VALUE)
            self.assertIsNone(engine._core)

    def test_finite_certificate_transports_every_label_orientation_and_edge_sign(self):
        rng = random.Random(251127)
        for _ in range(3):
            labels, order = list(range(FIVE.nv)), list(range(FIVE.ne))
            rng.shuffle(labels)
            rng.shuffle(order)
            graph = FIVE.reverse_edges(range(0, FIVE.ne, 2)).permute_edges(order).relabel(labels)
            sources = [None] * FIVE.nv
            for vertex, row in enumerate(FIVE_SOURCES):
                sources[labels[vertex]] = row
            with ForestReduction(graph, sources, tuple(labels[v] for v in (5, 6)), check=0) as engine:
                self.assertEqual(engine.moment(allow_series=False), perm_sign(order) * FIVE_VALUE)
                self.assertEqual(engine.last_stats['backend'], 'certified_five_loop_source_polynomial')
                self.assertEqual(engine.last_stats['finite_certificate']['status'], 'verified')
        sources = [list(row) for row in FIVE_SOURCES]
        sources[0][0] = 1
        sources[6][0] -= 1
        value, status = implementation._certified_bare(FIVE, sources, (5, 6))
        self.assertIsNone(value)
        self.assertEqual(status, 'source_support_outside_certificate')

    def test_stale_certificate_has_complete_real_series_and_finite_only_rejection(self):
        # This is the genuine unavailable K3,3+edge denominator-power-seven
        # table. Disabling the separately certified whole-source polynomial
        # must retain its entire core numerator in the exact fallback.
        with patch.object(implementation, '_CERTIFICATE_SHA256', '0' * 64):
            with ForestReduction(FIVE, FIVE_SOURCES, (5, 6), check=1) as engine:
                fallback = engine.moment()
                self.assertIsInstance(fallback, ExactForestSeries)
                self.assertIn('q=7', fallback.reason)
                self.assertTrue(fallback.components)
                self.assertEqual(fallback.constant, 0)
                explicit = engine.series()
                self.assertEqual(fallback.edges, explicit.edges)
                self.assertEqual(fallback.components, explicit.components)
                with self.assertRaises(MissingForestPeriod):
                    engine.moment(allow_series=False)
        data, status = implementation._verified_certificate()
        self.assertEqual(status, 'verified')
        dependency = implementation._ROOT / next(iter(data['implementation_sha256']))
        original = Path.read_bytes
        def changed(path):
            return original(path) + (b'\n' if path == dependency else b'')
        with patch.object(implementation, '_CERTIFICATE_CACHE', None), patch.object(Path, 'read_bytes', changed):
            value, status = implementation._verified_certificate()
        self.assertIsNone(value)
        self.assertTrue(status.startswith('certificate_dependency_changed:'))

    def test_direct_cube_route_retains_original_data_without_any_core(self):
        from htk.forest_cubes import ExactCubeSeries
        graph, _ = forest_graph((("A", "B"),))
        sources = ((1, 2), (5, 3), (-6, -5))
        insertion = polynomial(graph, [(1, {}), (F(2, 7), {(0, 1): 1})])
        expected = SPReduction(graph).integrate(sources, method='generic', check=1,
                       insertion=lambda currents: evaluate(insertion, currents))
        with ForestReduction(graph, sources, (1, 2), check=1) as engine:
            with patch.object(implementation, 'CoreCurrentPeriods', side_effect=AssertionError('global numerator')):
                series = engine.cube_series(polynomial=insertion, max_evaluations=1000)
                self.assertIsInstance(series, ExactCubeSeries)
                self.assertEqual(series.atlas.G.edges, graph.edges)
                self.assertEqual(series.atlas.sources, engine.sources)
                self.assertEqual(series.polynomial, insertion)
                terms = list(series.terms(1))
            self.assertEqual(sum(row.term for row in terms), terms[-1].partial_sum)
            self.assertLessEqual(terms[-1].lower, expected)
            self.assertGreaterEqual(terms[-1].upper, expected)
            self.assertEqual(terms[-1].evaluations, 15)
            self.assertTrue(engine.last_stats['original_graph'])
            self.assertIsNone(engine._core)

    def test_optional_normalization_budget_retains_entire_input_and_prepare_is_atomic(self):
        insertion = polynomial(HIGHER, [(1, {}), (F(-2, 7), {(0, 0): 1})])
        expected = F(55625, 432) - F(2, 7) * F(-31083, 80)
        with ForestReduction(HIGHER, HIGHER_SOURCES, (3, 4), max_terms=4, check=0) as engine:
            with self.assertRaises(ForestReductionBudgetExceeded):
                engine.normal_form(polynomial(HIGHER, [(1, {(0, 0): 10000})]))
        with ForestReduction(HIGHER, HIGHER_SOURCES, (3, 4), check=0) as engine:
            with patch.object(engine, 'normal_form', side_effect=ForestReductionBudgetExceeded('test')):
                self.assertEqual(engine.integrate(insertion), expected)
            self.assertEqual(engine.last_stats['kirchhoff_normalization'], 'deferred_at_budget')
        graph, sources, marked = passive_ear_fixture()
        with ForestReduction(graph, sources, marked, check=0) as engine:
            with patch.object(implementation, 'ForestEarMoments', side_effect=EarBudgetExceeded('test')):
                with self.assertRaises(EarBudgetExceeded):
                    engine._prepare()
            self.assertIsNone(engine._core)
            self.assertIsNone(engine._contractions)
            engine._prepare()
            self.assertEqual(engine.method, 'sourcefree_stars_then_sourced_ears')

    def test_exact_inputs_complete_tensor_budget_and_worker_cleanup(self):
        self.assertEqual(reduce_product_forest(HIGHER, HIGHER_SOURCES, (3, 4)), F(55625, 432))
        with ForestReduction(HIGHER, HIGHER_SOURCES, (3, 4), max_tensor_entries=2) as engine:
            with patch.object(engine, 'moment', side_effect=AssertionError('partial tensor computed')):
                with self.assertRaises(ForestReductionBudgetExceeded):
                    engine.tensor((0, 1), indices=('a', 'b'))
            for word in ({(0, 0): -1}, {(0, 2): 1}, {(0, 0): F(3, 2)}):
                with self.assertRaises((TypeError, ValueError)):
                    engine.moment(word)
            with self.assertRaises(TypeError):
                engine.integrate({(0,) * (2 * HIGHER.ne): 0.5})
        graph, sources, marked = passive_ear_fixture()
        engine = ForestReduction(graph, sources, marked, check=0)
        with engine:
            engine.moment(allow_series=False)
            self.assertIsNotNone(engine._ears._worker)
        self.assertIsNone(engine._ears._worker)


if __name__ == '__main__':
    unittest.main()
