"""Independent table controls and bounded exact-target regression checks."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symmetry'))
from family_target import (FamilyTarget, FamilyTargetBudgetExceeded,
                           IncompleteFamilyTarget, eliminate_target)


def table(q):
    result = json.loads((HERE.parent / ('fiveloop/tables/K33+e_q%d.json' % q)).read_text())
    return tuple(map(tuple, result['edges'])), {
        tuple(powers): {label: Fraction(value) for label, value in values.items() if Fraction(value)}
        for powers, values in result['values']}


def check_rows(test, system, values):
    for row, expected in system.rows:
        total = Counter()
        for column, coefficient in row.items():
            for label, value in values[system.reps[column][1]].items():
                total[label] += coefficient * value
        test.assertEqual({label: value for label, value in total.items() if value}, expected)


class TargetChecks(unittest.TestCase):
    def test_every_budget_preserves_the_entire_expression(self):
        rows = [({0: 1, 1: 2}, {'1': 7, 'z3': 3}),
                ({1: 1, 2: -1}, {'z5': 5}),
                ({2: 1}, {'1': 11})]
        powers = ((1, 1), (1, 2), (2, 1))
        values = {powers[0]: {'1': -15, 'z3': 3, 'z5': -10},
                  powers[1]: {'1': 11, 'z5': 5}, powers[2]: {'1': 11}}
        for budget in (0, 1, 2, 3):
            result = eliminate_target(rows, {0: 3, 1: -4, 2: 2}, powers,
                                      ((0, 1), (1, 2)), 3, max_pivots=budget)
            total = Counter(result.finite_coefficients)
            for term in result.residuals:
                for label, value in values[term.powers].items():
                    total[label] += term.coefficient * value
            self.assertEqual({k: v for k, v in total.items() if v},
                             {'1': -67, 'z3': 9, 'z5': -50})
            if result.residuals:
                with self.assertRaises(IncompleteFamilyTarget):
                    result.value
        self.assertTrue(result.complete)
        self.assertEqual(dict(result.value), {'1': -67, 'z3': 9, 'z5': -50})

    def test_inconsistent_rows_raise(self):
        with self.assertRaises(ArithmeticError):
            eliminate_target([({0: 1}, {'1': 1}), ({0: 1}, {'1': 2})],
                             {0: 1}, [(1,)], ((0, 1),), 3)

    def test_constant_labels_add_after_canonicalization(self):
        result = eliminate_target([({0: 1}, {'z5*z3': 2, 'z3*z5': -1})],
                                  {0: 1}, [(1,)], ((0, 1),), 3)
        self.assertEqual(dict(result.value), {'z3*z5': 1})

    def test_parallel_input_is_rejected_before_orbit_canonicalization(self):
        # gutil.edge_automorphisms is for simple input cores.  Without this
        # check the old repeated-edge map could change (2,3,3) to (3,3,3).
        for edges in (((0, 1), (0, 1), (0, 1)),
                      ((0, 1), (1, 0), (0, 1))):
            with self.assertRaisesRegex(ValueError, 'input core must be simple'):
                FamilyTarget(edges, [((2, 3, 3), 1)], 4)

    def test_exact_input_and_storage_budgets(self):
        edges, values = table(3)
        first, second = list(values)[:2]
        with self.assertRaises(ValueError):
            FamilyTarget(edges, [(first, 0.5)], 3)
        with self.assertRaises(FamilyTargetBudgetExceeded):
            FamilyTarget(edges, [(first, 1), (second, 1)], 3, max_variables=1)
        target = FamilyTarget(edges, [(first, 1), (first, -1)], 3)
        self.assertEqual(dict(target.reduce().value), {})
        with self.assertRaises(FamilyTargetBudgetExceeded):
            eliminate_target([({0: 1, 1: 1}, {})], {0: 1}, [(1,), (2,)],
                             ((0, 1),), 3, max_active_coefficients=1)

    def test_all_lower_q_local_rows_against_complete_tables(self):
        expected_counts = {3: {'D': 8, 'S': 8, 'Y': 6, 'X': 1},
                           4: {'D': 518, 'S': 135, 'Y': 129, 'X': 89}}
        for q in (3, 4):
            edges, values = table(q)
            target = FamilyTarget(edges, [(powers, 1) for powers in values], q)
            target.close_local()
            self.assertEqual({key: target.system.counts[key] for key in ('D', 'S', 'Y', 'X')},
                             expected_counts[q])
            self.assertEqual(len(target.system.rows), 23 if q == 3 else 871)
            check_rows(self, target.system, values)
            old = len(target.system.rows), len(target.system.reps), target.local_processed
            target.close_local()
            self.assertEqual((len(target.system.rows), len(target.system.reps), target.local_processed), old)

    def test_all_double_dimension_rows_against_independent_tables(self):
        edges, lower = table(3)
        higher_edges, higher = table(5)
        self.assertEqual(edges, higher_edges)
        target = FamilyTarget(edges, [], 5)
        target.add_dimension_rows(lower.items(), 3)
        self.assertEqual(target.system.counts['dimension'], 18)
        check_rows(self, target.system, higher)


if __name__ == '__main__':
    unittest.main(verbosity=2)
