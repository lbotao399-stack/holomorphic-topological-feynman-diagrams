"""Exact terminal and whole-expression controls; no period table generation.

The four K5 values below were obtained by the independent symmetry routes
stored in check_k5_examples.log before K5 tables became direct terminals.
The cancellation control disables every finite terminal: it must use the
Fourier identity with its full rational normalization, across two graphs.
"""
from fractions import Fraction as F
from math import factorial, prod
import unittest
from unittest.mock import patch

import engine
from combination import (CombinationBudgetExceeded, IncompleteCombination,
                         reduce_period_combination)


K5 = tuple((i, j) for i in range(5) for j in range(i + 1, 5))
ARCHIVED = (
    ((5, 5, 6, 1, 6, 4, 5, 6, 5, 5), 8,
     {(): F(105943, 1569209241600), ('z3',): F(-17, 302702400)}),
    ((4, 4, 6, 7, 6, 1, 5, 5, 5, 5), 8,
     {(): F(-62311, 224172748800), ('z3',): F(1, 4324320)}),
    ((6, 5, 5, 4, 3, 1, 7, 7, 7, 3), 8,
     {(): F(-2759219, 2615348736000), ('z3',): F(443, 504504000)}),
    ((4, 2, 5, 5, 5, 6, 5, 2, 5, 3), 7,
     {(): F(3251321, 64665216000), ('z3',): F(52643, 449064000),
      ('z5',): F(-23, 124740)}),
)
FOURIER_GRAPH = ((0, 1), (1, 2), (0, 3), (3, 4), (1, 4), (2, 4),
                 (0, 5), (1, 5), (2, 5), (3, 5), (4, 5))
FOURIER_POWERS = (1, 6, 5, 4, 2, 4, 6, 5, 2, 3, 4)


class ExactConformalCombination(unittest.TestCase):
    def test_new_terminals_match_archived_independent_routes(self):
        for powers, q, expected in ARCHIVED:
            with self.subTest(powers=powers, q=q):
                direct = engine.eval_I(K5, powers, q)
                self.assertEqual(direct, (expected, 'K5'))
                completed = engine.period(K5, powers, q, budget=1)
                self.assertIsNotNone(completed)
                self.assertEqual(completed[0], expected)
                self.assertIn('K5-e', completed[1][-1])

    def test_previously_missing_completed_orbit(self):
        key = (7, ((0, 1, -2), (0, 2, 1), (0, 3, 4), (0, 4, 2),
                   (0, 5, 4), (0, 6, 5), (1, 2, 6), (1, 3, 3),
                   (1, 5, 3), (1, 6, 4), (2, 3, 3), (2, 4, 4),
                   (3, 5, 2), (3, 6, 2), (4, 5, 5), (4, 6, 3)))
        vertices, weights = engine.unkey(key)
        self.assertTrue(engine.balanced(vertices, weights, 7))
        actual = engine.solve(vertices, weights, 7, budget=40)
        self.assertIsNotNone(actual)
        self.assertEqual(actual[0], {
            (): F(1269429229, 238381852262400),
            ('z3',): F(306217, 21499084800),
            ('z5',): F(-49771, 2299207680),
        })

    def test_product_labels_and_read_only_missing_tables(self):
        self.assertEqual(engine.from_cores({'1': F(3), 'z3*z3': F(2), 'z5*z3': F(-1)}),
                         {(): F(3), ('z3', 'z3'): F(2), ('z3', 'z5'): F(-1)})
        self.assertEqual(engine.from_cores({'z3*z5': F(2), 'z5*z3': F(3)}),
                         {('z3', 'z5'): F(5)})
        self.assertEqual(engine.from_cores({'z3*z5': F(2), 'z5*z3': F(-2)}), {})
        # Load the lazy bridge before patching the relation-system builder.
        engine.eval_I(K5, ARCHIVED[-1][0], 7)
        import per5
        with patch.object(per5.msys.MultiSystem, 'build',
                          side_effect=AssertionError('unexpected table construction')):
            self.assertIsNone(engine.eval_I(K5, (7,) * 10, 9))
            self.assertIsNone(engine.eval_I(engine.cores.PRISM, (5,) * 9, 12))
            self.assertIsNone(engine.eval_I(engine.cores.K33, (5,) * 9, 12))
        with patch.object(engine.os.path, 'isfile', return_value=False), \
                patch.object(engine, '_orig_table', side_effect=AssertionError('unexpected table access')):
            with self.assertRaises(engine.NotTabulated):
                engine._guarded_table('prism', 7)

    def test_finite_whole_combination(self):
        coefficients = (F(2), F(-3), F(5, 7))
        terms = [(K5, row[0], coefficient)
                 for row, coefficient in zip(ARCHIVED[:3], coefficients)]
        expected = {}
        for row, coefficient in zip(ARCHIVED[:3], coefficients):
            for label, value in row[2].items():
                text = '*'.join(label) or '1'
                expected[text] = expected.get(text, F(0)) + coefficient * value
        result = reduce_period_combination(terms, 8, max_states=10)
        self.assertTrue(result.complete)
        self.assertEqual(dict(result.value), expected)
        self.assertFalse(result.residuals)

    def test_nontrivial_fourier_cancellation_without_any_finite_lookup(self):
        terms = [(K5, ARCHIVED[-1][0], F(1)),
                 (FOURIER_GRAPH, FOURIER_POWERS, F(-16, 5))]
        with patch.object(engine, 'terminal', return_value=None):
            retained = reduce_period_combination(terms, 7, max_states=0)
            reduced = reduce_period_combination(terms, 7, max_states=40, max_nodes=500)
        self.assertFalse(retained.complete)
        self.assertEqual(len(retained.residuals), 2)
        with self.assertRaises(IncompleteCombination):
            _ = retained.value
        self.assertTrue(reduced.complete)
        self.assertEqual(dict(reduced.value), {})
        self.assertEqual(reduced.statistics['terminal_values'], 0)
        self.assertEqual(reduced.statistics['cancelled_unresolved_components'], 1)
        self.assertGreater(reduced.statistics['merged_states'], 0)

    def test_search_budget_preserves_the_complete_residual_coefficient(self):
        powers = ARCHIVED[-1][0]
        coefficient = F(-13, 17)
        with patch.object(engine, 'terminal', return_value=None):
            result = reduce_period_combination([(K5, powers, coefficient)], 7,
                                               max_states=1, max_nodes=1)
        self.assertFalse(result.complete)
        self.assertEqual(dict(result.finite_coefficients), {})
        self.assertEqual(len(result.residuals), 1)
        expected = coefficient * F(prod(factorial(6 - power) for power in powers), factorial(6))
        self.assertEqual(result.residuals[0].coefficient, expected)
        self.assertGreater(result.statistics['node_budget_skips'], 0)

    def test_exact_input_validation_and_input_budgets(self):
        row = (K5, ARCHIVED[-1][0], F(1))
        for options in ({'max_states': -1}, {'max_nodes': 0}, {'max_edges': True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                reduce_period_combination([row], 7, **options)
        with self.assertRaises(ValueError):
            reduce_period_combination([(K5, row[1], 0.5)], 7)
        with self.assertRaises(ValueError):
            reduce_period_combination([(K5, (1,) * 10, F(1))], 7)
        with self.assertRaises(CombinationBudgetExceeded):
            reduce_period_combination([row], 7, max_edges=9)
        with self.assertRaises(CombinationBudgetExceeded):
            reduce_period_combination([row, row], 7, max_input_terms=1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
