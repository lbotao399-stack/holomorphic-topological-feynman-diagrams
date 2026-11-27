"""Focused exact controls for coordinate planning and shared row proofs."""
from fractions import Fraction
from pathlib import Path
import json
import random
import sys
import unittest

from source_family import HERE, ROOT, BASIS, EDGES, five_loop_source_kernel, sources_at
from htk import Graph, fan
from htk.forest_periods import CoreCurrentPeriods, ForestPeriodBudgetExceeded
from tree_minor import prepare_source_numerator
from batch_witness import eliminate_batch_with_witness
from family_target import FamilyTargetBudgetExceeded


class TreeMinorChecks(unittest.TestCase):
    def compare(self, graph, sources):
        original = CoreCurrentPeriods(graph, sources, check=1, max_terms=100000)
        planned = CoreCurrentPeriods(graph, sources, check=2, max_terms=100000)
        try:
            original._source_numerator()
            statistics = prepare_source_numerator(planned, max_tree_subsets=100, max_plans=8)
            self.assertEqual(original._bare, planned._bare)
            self.assertIsNone(planned._worker)
            self.assertTrue(statistics['exhaustive_tree_scan'])
            self.assertEqual(statistics['independent_density_checks'], 2)
        finally:
            original.close()
            planned.close()

    def test_full_rational_sources_and_orientation(self):
        graph = fan(3)
        sources = ((Fraction(1, 2), 2), (-1, Fraction(1, 3)), (3, -2),
                   (Fraction(2, 5), 1), (Fraction(-29, 10), Fraction(-4, 3)))
        self.compare(graph, sources)
        self.compare(graph.reverse_edges([0, 3]).permute_edges((6, 4, 2, 0, 5, 3, 1)), sources)

    def test_parallel_edges_keep_individual_positions(self):
        graph = Graph(4, ((0, 1), (0, 1), (1, 2), (2, 3), (0, 3)))
        self.compare(graph, ((1, 2), (3, -1), (-2, 4), (-2, -5)))

    def test_rank_zero_is_exact_zero(self):
        core = CoreCurrentPeriods(fan(3), ((0, 0),) * 5, check=2)
        try:
            statistics = prepare_source_numerator(core)
            self.assertFalse(core._bare)
            self.assertTrue(statistics['constant_rank_zero'])
        finally:
            core.close()

    def test_exhausted_budget_never_installs_partial_cache(self):
        core = CoreCurrentPeriods(fan(3), ((1, 2), (3, -1), (-2, 4), (2, 1), (-4, -6)),
                                  check=0, max_terms=1)
        try:
            with self.assertRaises(ForestPeriodBudgetExceeded):
                prepare_source_numerator(core, max_tree_subsets=2, max_plans=2)
            self.assertIsNone(core._bare)
            self.assertIsNone(core._context)
        finally:
            core.close()

    def test_both_structural_source_zeros(self):
        for sources, pointwise in [(((0,0),(0,0),(0,0),(1,2),(0,0),(3,-1),(-4,-1)), True),
                                   (((0,0),(1,0),(0,0),(0,0),(0,0),(0,1),(-1,-1)), False)]:
            core = CoreCurrentPeriods(Graph(7, EDGES), sources, check=1, max_terms=750000)
            try:
                prepare_source_numerator(core)
                if pointwise:
                    self.assertFalse(core._bare)
                else:
                    values = core._core_polynomials({(0,) * 22: Fraction(1)})
                    self.assertFalse(any(values.values()))
            finally:
                core.close()

    def test_tree_and_allocation_guards(self):
        core = CoreCurrentPeriods(fan(3), ((0, 0),) * 5)
        try:
            for bad in (0, -1, 2.0, True):
                with self.assertRaises(ValueError):
                    prepare_source_numerator(core, max_plans=bad)
            with self.assertRaises(ValueError):
                prepare_source_numerator(core, tree_candidates=((0, 0, 1, 2),))
            with self.assertRaises(ValueError):
                prepare_source_numerator(core, tree_candidates=())
            with self.assertRaises(ForestPeriodBudgetExceeded):
                prepare_source_numerator(core, tree_candidates=(tuple(core.G.spanning_tree()),) * 2,
                                         max_tree_subsets=1)
            tree_core = CoreCurrentPeriods(Graph(2, ((0, 1),)), ((1, 2), (-1, -2)))
            try:
                with self.assertRaises(ValueError):
                    prepare_source_numerator(tree_core)
            finally:
                tree_core.close()
            core.max_matrix_entries = 1
            with self.assertRaises(ForestPeriodBudgetExceeded):
                prepare_source_numerator(core)
        finally:
            core.close()


class BatchWitnessChecks(unittest.TestCase):
    def test_complete_and_budget_stopped_targets(self):
        randomizer = random.Random(5112026)
        for trial in range(8):
            exact = [Fraction(randomizer.randrange(-9, 10), randomizer.randrange(1, 8))
                     for _ in range(7)]
            rows = []
            for _ in range(12):
                row = {i: Fraction(randomizer.randrange(-3, 4)) for i in randomizer.sample(range(7), 3)}
                row = {i: value for i, value in row.items() if value}
                rows.append((row, {'1': sum((value * exact[i] for i, value in row.items()), Fraction())}))
            targets = [{i: Fraction(randomizer.randrange(-5, 6), randomizer.randrange(1, 5))
                        for i in range(7)} for _ in range(4)] + [{}]
            representatives = [(i + 1,) for i in range(7)]
            for pivots in (0, 1, 7):
                results = eliminate_batch_with_witness(rows, targets, representatives, ((0, 1),), 7,
                                                       max_pivots=pivots, max_history_entries=1000)
                for target, result in zip(targets, results):
                    expected = sum((value * exact[i] for i, value in target.items()), Fraction())
                    observed = result.finite_coefficients.get('1', 0)
                    observed += sum((term.coefficient * exact[term.powers[0] - 1]
                                     for term in result.residuals), Fraction())
                    self.assertEqual(observed, expected)
                    if pivots == 0:
                        self.assertFalse(result.row_weights)
                    self.assertFalse(result.statistics['stopped_by'] is None and pivots == 0 and target)
                self.assertTrue(results[-1].complete)

    def test_history_and_initial_storage_guards(self):
        rows = [({0: 1, 1: 1}, {'1': 3}), ({0: 2, 1: -1}, {})]
        targets = [{0: 1}, {1: 1}]
        with self.assertRaises(FamilyTargetBudgetExceeded):
            eliminate_batch_with_witness(rows, targets, [(1,), (2,)], ((0, 1),), 7,
                                          max_active_coefficients=1)
        results = eliminate_batch_with_witness(rows, targets, [(1,), (2,)], ((0, 1),), 7,
                                               max_history_entries=0)
        self.assertTrue(all(result.statistics['stopped_by'] == 'max_history_entries' for result in results))
        self.assertTrue(all(not result.row_weights for result in results))


class CovariantFormulaChecks(unittest.TestCase):
    def test_full_source_evaluator_support_and_degenerate_sources(self):
        if not (HERE / 'source_certificate.json').exists():
            self.skipTest('certificate generation is required')
        self.assertEqual(five_loop_source_kernel(sources_at(1,1)),
                         {'1': Fraction(1187,435456000), 'z3': Fraction(-17,9072000)})
        self.assertEqual(five_loop_source_kernel(((0,0),)*7), 0)
        original = sources_at(3,1)
        swapped = (*original[:5], original[6], original[5])
        self.assertEqual(five_loop_source_kernel(swapped), -five_loop_source_kernel(original))
        bad = [list(row) for row in original]
        bad[0][0] = 1
        bad[6][0] -= 1
        with self.assertRaises(ValueError):
            five_loop_source_kernel(bad)
        bad = [list(row) for row in original]
        bad[1][0] = 1.0
        with self.assertRaises((TypeError,ValueError)):
            five_loop_source_kernel(bad)


    def test_ten_real_coefficients_are_nonzero_without_zeta_independence(self):
        path = HERE / 'source_certificate.json'
        if not path.exists():
            self.skipTest('certificate generation is required')
        certificate = json.loads(path.read_text())
        bounds = {}
        for power in (3,5):
            partial = sum((Fraction(1, k ** power) for k in range(1,65)), Fraction())
            bounds['z%d' % power] = (partial + Fraction(1,(power-1)*65**(power-1)),
                                     partial + Fraction(1,(power-1)*64**(power-1)))
        for index, values in enumerate(certificate['finite_coefficients_in_basis_order']):
            lower = upper = Fraction(values.get('1',0))
            for label, (first,last) in bounds.items():
                coefficient = Fraction(values.get(label,0))
                lower += coefficient * (first if coefficient >= 0 else last)
                upper += coefficient * (last if coefficient >= 0 else first)
            if BASIS[index] == (0,5):
                self.assertEqual((lower,upper),(0,0))
            else:
                self.assertTrue(lower > 0 or upper < 0)

    def test_zeta5_factorization_and_source_specializations(self):
        path = HERE / 'source_certificate.json'
        if not path.exists():
            self.skipTest('run certify_sources.py --write before the finite formula check')
        certificate = json.loads(path.read_text())
        import sympy
        first, second, third = sympy.symbols('bracket_13 bracket_3_difference bracket_1_difference')
        zeta5 = sum(sympy.Rational(values.get('z5', '0')) * first ** (5-i-j) * second ** i * third ** j
                    for (i,j),values in zip(BASIS, certificate['finite_coefficients_in_basis_order']))
        factorized = (second+third)*(18*first**4-33*first**2*second**2
                         -39*first**2*second*third-24*first**2*third**2
                         +2*second**3*third+3*second**2*third**2+second*third**3)/120960
        self.assertEqual(sympy.expand(zeta5-factorized), 0)
        self.assertEqual(zeta5.subs({first:1,second:-3,third:1}), sympy.Rational(1,280))
        self.assertEqual(zeta5.subs(third,-second), 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
