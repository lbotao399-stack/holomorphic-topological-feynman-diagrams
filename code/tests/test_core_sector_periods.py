"""Independent density, calculus, exact-integral and certificate controls."""
from dataclasses import replace
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path
import runpy
import sys
import unittest

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from htk import brute
from htk.graph import Graph
from htk.core_sector_periods import (CoreSector, CubePeriodRelations,
    SectorBudgetExceeded, ReductionNotFound, check_convergence)


def evaluate(poly, point):
    return sum((c * sp.prod(x ** a for x, a in zip(point, powers))
                for powers, c in poly.items()), F())


class CoreSectorTests(unittest.TestCase):
    def test_original_kirchhoff_determinant_and_complete_jacobian(self):
        cases = [
            (Graph(3, [(0, 1), (1, 2), (0, 2)]), (1, 1, 2), 4),
            (Graph(2, [(0, 1)] * 3), (2, 2, 2), 3),
            (Graph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
             (1,) * 6, 2),
        ]
        checked = 0
        for graph, alpha, q in cases:
            orders = [tuple(range(graph.ne)), tuple(reversed(range(graph.ne)))]
            for order in orders:
                sector = CoreSector(graph, order)
                form = sector.pullback(alpha, q)
                x = sp.symbols('x:' + str(graph.ne - 1))
                parameter = [None] * graph.ne
                for i, edge in enumerate(order):
                    parameter[edge] = sp.prod(x[i:])
                jacobian = sp.Matrix([parameter[edge] for edge in order[:-1]]).jacobian(x).det()
                for offset in (1, 4):
                    point = tuple(F(i + offset, i + offset + 2) for i in range(graph.ne - 1))
                    substitutions = dict(zip(x, point))
                    t = [F(value.subs(substitutions)) if hasattr(value, 'subs') else F(value)
                         for value in parameter]
                    original_U = brute.symanzik_U(graph, t)
                    factor = sp.prod(value ** loops for value, loops in zip(point, sector.prefix_loops))
                    self.assertEqual(original_U,
                                     factor * evaluate(sector.relations.denominator, point))
                    original = sp.prod(value ** (a - 1) for value, a in zip(t, alpha))
                    original *= abs(F(jacobian.subs(substitutions))) / original_U ** q
                    transformed = evaluate(form.numerator, point)
                    transformed /= evaluate(sector.relations.denominator, point) ** form.power
                    self.assertEqual(original, transformed)
                    checked += 1
        self.assertEqual(checked, 12)

    def test_all_subset_convergence_and_all_shifted_weights(self):
        graph = Graph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
        sector = CoreSector(graph, range(6))
        self.assertEqual(check_convergence(graph, (1,) * 6, 2),
                         {'proper_subsets': 62, 'minimum_integer_margin': 1})
        shifts = sector.shifted_weights((1,) * 6, 2)
        self.assertEqual(len(shifts), 16)
        self.assertEqual(sector.shifted_weights(iter((1,) * 6), 2), shifts)
        for alpha in shifts:
            self.assertGreater(check_convergence(graph, alpha, 3)['minimum_integer_margin'], 0)
        # Independent integrand check of U multiplication at a rational point.
        t = [F(i + 1, i + 3) for i in range(6)]
        U = brute.symanzik_U(graph, t)
        self.assertEqual(sum(sp.prod(value ** (a - 1) for value, a in zip(t, alpha))
                             / U ** 3 for alpha in shifts), 1 / U ** 2)
        with self.assertRaises(ValueError):
            check_convergence(Graph(2, [(0, 1)] * 3), (1, 1, 2), 2)

    def test_complete_boundary_derivative_against_symbolic_calculus(self):
        x, y = sp.symbols('x y')
        cube = CubePeriodRelations({(0, 0): 2, (1, 0): F(3, 5), (1, 1): 1})
        denominator = 2 + F(3, 5) * x + x * y
        for power in (0, 1, 4):
            for coordinate, variable in enumerate((x, y)):
                form = cube.boundary_derivative((2, 3), power, coordinate)
                constructed = evaluate(form.numerator, (x, y)) / denominator ** form.power
                direct = sp.diff(variable * (1 - variable) * x ** 2 * y ** 3
                                 / denominator ** power, variable)
                self.assertEqual(sp.cancel(constructed - direct), 0)

    def test_nonrational_one_variable_reduction_and_full_witness(self):
        cube = CubePeriodRelations({(0,): 1, (1,): 1})
        target = cube.monomial((7,), 5)
        masters = [cube.monomial((0,)), cube.monomial((0,), 1)]
        certificate = cube.search(target, masters, max_primitive_power=4,
                                  max_primitive_degree=8)
        self.assertEqual(certificate.coefficients, (F(1553, 64), F(-35)))
        x = sp.symbols('x', positive=True)
        original = sp.integrate(x ** 7 / (1 + x) ** 5, (x, 0, 1))
        self.assertEqual(sp.simplify(original - (sp.Rational(1553, 64) - 35 * sp.log(2))), 0)
        self.assertTrue(cube.verify(certificate))
        record = json.loads(json.dumps(certificate.record()))
        self.assertEqual(cube.read_certificate(record).coefficients, certificate.coefficients)

    def test_complete_three_edge_banana_by_all_six_sectors(self):
        graph = Graph(2, [(0, 1)] * 3)
        total = F()
        for order in permutations(range(3)):
            sector = CoreSector(graph, order)
            cube = sector.relations
            certificate = cube.search(sector.pullback((2, 2, 2), 3),
                                      [cube.monomial((0, 0))],
                                      max_primitive_power=2, max_primitive_degree=4)
            self.assertEqual(certificate.coefficients, (F(1, 12),))
            total += certificate.coefficients[0]
        # Independent projective inversion gives a triangle Dirichlet integral.
        self.assertEqual(total, F(1, 2))

    def test_two_variable_nonrational_reduction(self):
        cube = CubePeriodRelations({(0, 0): 1, (1, 1): 1})
        certificate = cube.search(cube.monomial((2, 2), 3),
            [cube.monomial((0, 0)), cube.monomial((0, 0), 1), cube.monomial((1, 0), 1)],
            max_primitive_power=2, max_primitive_degree=4)
        self.assertEqual(certificate.coefficients, (F(-1, 2), F(1), F(-3, 4)))
        # Independent analytic integration:
        # J1=pi^2/12, J2=log2, J3=1/4+log2/2;
        # x^2*y^2/(1+xy)^3 = U^-1 -2U^-2 +U^-3.
        x, y = sp.symbols('x y', positive=True)
        for power, expected in ((2, sp.log(2)), (3, sp.Rational(1, 4) + sp.log(2) / 2)):
            first = sp.integrate(1 / (1 + x * y) ** power, (y, 0, 1))
            self.assertEqual(sp.simplify(sp.integrate(first, (x, 0, 1)) - expected), 0)
        self.assertTrue(cube.verify(certificate))

    def test_unresolved_search_retains_whole_residual(self):
        cube = CubePeriodRelations({(0,): 1, (1,): 1})
        with self.assertRaises(ReductionNotFound) as caught:
            cube.search(cube.monomial((0,), 1), [],
                        max_primitive_power=0, max_primitive_degree=0)
        self.assertEqual(dict(caught.exception.residual.numerator), {(1,): F(1), (2,): F(2)})
        self.assertEqual(caught.exception.residual.power, 1)
        self.assertEqual(caught.exception.statistics['generators'], 1)

    def test_corrupted_witness_is_rejected(self):
        cube = CubePeriodRelations({(0,): 1, (1,): 1})
        certificate = cube.search(cube.monomial((0,), 2), [cube.monomial((0,))],
                                  max_primitive_power=1, max_primitive_degree=2)
        bad = replace(certificate, coefficients=(certificate.coefficients[0] + 1,))
        self.assertFalse(cube.verify(bad))
        record = certificate.record()
        record['primitives'] = []
        with self.assertRaises(ArithmeticError):
            cube.read_certificate(record)

    def test_certificate_replay_binds_to_requested_target(self):
        directory = Path(__file__).resolve().parents[2] / 'explore/circuit-general/core-sector-reduction'
        replay = runpy.run_path(str(directory / 'certify.py'))
        data = json.loads((directory / 'reductions.json').read_text())
        record = data['cases'][0]['witness']
        x = sp.symbols('x')
        # This is a valid full identity with the same integrated value, but a
        # different target density. A target-specific certificate must reject it.
        numerator = sp.Poly(x ** 7 + (1 - 2 * x) * (1 + x) ** 5, x)
        record['target']['numerator'] = [[list(powers), str(coefficient)]
                                         for powers, coefficient in numerator.terms()]
        record['primitives'].append(['1', 0, [0], 0])
        cube = CubePeriodRelations({(0,): 1, (1,): 1})
        self.assertTrue(cube.verify(cube.read_certificate(record)))
        with self.assertRaisesRegex(ArithmeticError, 'intended target'):
            replay['verify'](data)

    def test_declared_budgets_and_exact_input_rejection(self):
        cube = CubePeriodRelations({(0,): 1, (1,): 1})
        with self.assertRaises(SectorBudgetExceeded):
            cube.search(cube.monomial((0,), 1), [], max_primitive_power=1,
                        max_primitive_degree=1, max_generators=1)
        # No pivot overlaps this target: the input still has to respect the
        # elimination budget before any row operation is attempted.
        with self.assertRaises(SectorBudgetExceeded):
            cube.search(cube.form({(a,): 1 for a in range(100, 106)}, 0), [],
                        max_primitive_power=0, max_primitive_degree=0,
                        max_elimination_entries=4)
        with self.assertRaises(SectorBudgetExceeded):
            CoreSector(Graph(3, [(0, 1), (1, 2), (0, 2)]), (0, 1, 2), max_tree_candidates=1)
        with self.assertRaises(SectorBudgetExceeded):
            check_convergence(Graph(2, [(0, 1)] * 3), (2, 2, 2), 3, max_subsets=1)
        with self.assertRaises(ValueError):
            CubePeriodRelations({(0,): 1, (1,): -1})
        with self.assertRaises(TypeError):
            cube.monomial((1,), coefficient=0.5)
        with self.assertRaises(ValueError):
            cube.monomial((1,), power=1.5)
        with self.assertRaises(ValueError):
            CoreSector(Graph(3, [(0, 1), (1, 2), (0, 2)]), (0, 1, 1))
        with self.assertRaises(ValueError):
            check_convergence(Graph(1, []), (), 1)
        with self.assertRaises(SectorBudgetExceeded):
            CoreSector(Graph(3, [(0, 1), (1, 2), (0, 2)]), (0, 1, 2), max_terms=1)
        limited = CubePeriodRelations({(0,): 1}, max_terms=3)
        first = limited.form({(0,): 1, (1,): 1, (2,): 1}, 0)
        second = limited.form({(0,): 1, (3,): 1, (4,): 1}, 0)
        with self.assertRaises(SectorBudgetExceeded):
            limited.search(first, [first, second], max_primitive_power=0,
                           max_primitive_degree=0)
        # Constant denominators do not need a list with one entry for every
        # integer pole order: exact binary powering keeps this sparse.
        constant = CubePeriodRelations({(0,): 1})
        self.assertEqual(constant.lift(constant.monomial((0,)), 10 ** 9), {(0,): F(1)})


if __name__ == '__main__':
    unittest.main()
