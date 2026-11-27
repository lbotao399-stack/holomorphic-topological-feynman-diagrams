"""Exact proof, budget, and malformed-annihilator regression controls."""
from dataclasses import replace
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'symmetry'))
from annihilator_target import (AnnihilatorTarget, UnavailableAnnihilators,
                                 load_annihilators, verify_annihilators)
from family_target import IncompleteFamilyTarget
from relation_witness import eliminate_with_witness, verify_row_witness


class WitnessChecks(unittest.TestCase):
    def test_finite_and_partial_witnesses_reproduce_the_entire_target(self):
        rows = [({0: 1, 1: 1}, {'1': 3}), ({1: 1, 2: 1}, {'1': 5}),
                ({0: 1, 2: 1}, {'1': 4})]
        target = {0: 2, 1: -3, 2: 5}
        representatives = ((1,), (2,), (3,))
        for budgets in ({}, {'max_pivots': 0}, {'max_pivots': 1},
                        {'max_operations': 0}, {'max_history_entries': 0}):
            result = eliminate_with_witness(rows, target, representatives, ((0, 1),), 3, **budgets)
            value = result.finite_coefficients.get('1', 0)
            value += sum(term.coefficient * term.powers[0] for term in result.residuals)
            self.assertEqual(value, 11)
            residual = {representatives.index(term.powers): term.coefficient for term in result.residuals}
            self.assertTrue(verify_row_witness(rows, target, residual,
                                              result.finite_coefficients, result.row_weights))
            if result.residuals:
                with self.assertRaises(IncompleteFamilyTarget):
                    result.value
        complete = eliminate_with_witness(rows, target, representatives, ((0, 1),), 3)
        self.assertEqual(dict(complete.value), {'1': 11})
        tampered = dict(complete.row_weights)
        tampered[next(iter(tampered))] += 1
        with self.assertRaises(ArithmeticError):
            verify_row_witness(rows, target, {}, complete.finite_coefficients, tampered)

    def test_fill_budget_stops_before_a_whole_pivot(self):
        rows = []
        for start in range(7):
            columns = [(start + shift) % 7 for shift in range(4)]
            rows.append(({column: 1 for column in columns},
                         {'1': sum(column + 1 for column in columns)}))
        target = {column: (-1) ** column for column in range(7)}
        result = eliminate_with_witness(rows, target, [(index + 1,) for index in range(7)],
                                        ((0, 1),), 3, max_active_coefficients=28)
        self.assertEqual(result.statistics['stopped_by'], 'max_active_coefficients')
        self.assertEqual(result.statistics['pivots'], 0)
        self.assertEqual(dict(result.row_weights), {})
        self.assertEqual({term.powers[0] - 1: term.coefficient for term in result.residuals}, target)

    def test_inconsistent_finite_rows_raise(self):
        with self.assertRaises(ArithmeticError):
            eliminate_with_witness([({0: 1}, {'1': 1}), ({0: 1}, {'1': 2})],
                                    {0: 1}, [(1,)], ((0, 1),), 3)

    def test_invalid_polynomial_or_degree_cannot_be_certified(self):
        identity = json.loads((HERE.parent / 'forest-target-reduction/dimension_certificate.json').read_text())
        edges = tuple(map(tuple, identity['scalar_core']['edges']))
        generators = load_annihilators(HERE / 'annihilators.json', edges)
        first = generators[0]
        vector = [dict(polynomial) for polynomial in first.vector_fields]
        index = next(index for index, polynomial in enumerate(vector) if polynomial)
        powers = next(iter(vector[index]))
        vector[index][powers] += 1
        with self.assertRaises(ArithmeticError):
            verify_annihilators(edges, [replace(first, vector_fields=tuple(vector))])
        with self.assertRaises(ValueError):
            verify_annihilators(edges, [replace(first, degree=first.degree + 1)])
        with self.assertRaises(UnavailableAnnihilators):
            load_annihilators(HERE / 'annihilators.json', edges[:-1])

    def test_candidate_budget_preserves_every_archived_target_coefficient(self):
        identity = json.loads((HERE.parent / 'forest-target-reduction/dimension_certificate.json').read_text())
        edges = tuple(map(tuple, identity['scalar_core']['edges']))
        terms = [(tuple(powers), Fraction(coefficient))
                 for powers, coefficient in identity['all_residual_scalar_periods']]
        generators = load_annihilators(HERE / 'annihilators.json', edges)
        target = AnnihilatorTarget(edges, terms, 7, generators)
        target.generate(max_candidates=0)
        result = target.reduce(max_pivots=0)
        self.assertEqual(result.statistics['generation_stopped_by'], 'max_candidates')
        self.assertEqual({term.powers: term.coefficient for term in result.residuals}, dict(terms))
        self.assertEqual(dict(result.row_weights), {})
        self.assertEqual(dict(result.finite_coefficients), {})
        with self.assertRaises(IncompleteFamilyTarget):
            result.value


if __name__ == '__main__':
    unittest.main(verbosity=2)
