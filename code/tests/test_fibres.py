"""Independent exact checks of the chain-fibre integral.

Run ``python3 code/tests/test_fibres.py``.  The main oracle reconstructs the
old whole-graph numerator and integrates every monomial with the elementary
factorial Dirichlet formula.  It uses neither the new cubature weights nor
the contracted determinant.  The three-loop H=3 replay takes several seconds.
"""
import math
import os
import sys
import unittest
from fractions import Fraction as F

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "cores"))

import numer
import htk
from htk import brute
from htk.dirichlet import CubatureLimitError
from htk.fibres import integrate_fibre
from htk.graph import polygon, theta
from htk.series import SeriesReduction
from test_series import subdivide, complete_edges, rational_sources


def old_numerator_fibre(G, lam, reduction, resistances, powers=None):
    """Whole-graph reconstruction followed by the factorial monomial formula."""
    insertions = [key for key, power in (powers or {}).items() for _ in range(power)]
    polynomial, q = numer.numerator(G, lam, insertions=insertions)
    centre_t = [F(0)] * G.ne
    for chain, R in zip(reduction.chains, resistances):
        for edge in chain.edges:
            centre_t[edge] = F(R, chain.length)
    U = brute.symanzik_U(G, centre_t)
    value = F(0)
    for exponents, coefficient in polynomial.items():
        contribution = coefficient
        for chain, R in zip(reduction.chains, resistances):
            degree = sum(exponents[edge] for edge in chain.edges)
            numerator = math.prod(math.factorial(exponents[edge]) for edge in chain.edges)
            denominator = math.factorial(degree + chain.length - 1)
            contribution *= F(numerator, denominator) * F(R) ** (degree + chain.length - 1)
        value += contribution
    return value / U ** q


class FibreTests(unittest.TestCase):
    def test_full_theta_numerator_and_insertion(self):
        G = theta(1, 3, 3)[0]
        lam = rational_sources(G, 3, 217)
        reduction = SeriesReduction(G, lam)
        resistances = [1, 2, 3]
        for powers in (None, {(0, 0): 1, (G.ne - 1, 2): 1}):
            with self.subTest(powers=powers):
                expected = old_numerator_fibre(G, lam, reduction, resistances, powers)
                actual = integrate_fibre(reduction, resistances, powers)
                self.assertEqual(actual.value, expected)
                self.assertNotEqual(actual.value, 0)
                self.assertGreater(actual.nodes, 0)

    def test_generic_three_loop_tetrahedron_full_numerator(self):
        G = subdivide(4, complete_edges(4), [3, 2, 1, 1, 2, 1])
        lam = rational_sources(G, 3, 217)
        reduction = SeriesReduction(G, lam)
        resistances = [1, 2, 3, 2, 1, 3]
        expected = old_numerator_fibre(G, lam, reduction, resistances)
        actual = integrate_fibre(reduction, resistances)
        self.assertEqual(actual.value, expected)
        self.assertNotEqual(actual.value, 0)

    def test_polygon_is_a_completed_period(self):
        G = polygon(4)
        lam = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        reduction = SeriesReduction(G, lam)
        self.assertEqual(len(reduction.chains), 1)
        for powers in (None, {(0, 0): 2, (2, 1): 1}):
            expected = htk.kernel_exact(G, lam, powers=powers)
            actual = integrate_fibre(reduction, [1], powers)
            self.assertEqual(actual.value, expected)

    def test_gauge_reversal_and_permutation(self):
        G = subdivide(4, complete_edges(4), [2, 1, 1, 1, 1, 1])
        lam = rational_sources(G, 2, 317)
        order = [3, 1, 6, 0, 5, 2, 4]
        changed = G.permute_edges(order).reverse_edges([0, 3, 5])
        for graph in (G, changed):
            reduction = SeriesReduction(graph, lam)
            resistances = [1 + i % 3 for i in range(len(reduction.chains))]
            expected = old_numerator_fibre(graph, lam, reduction, resistances)
            actual = integrate_fibre(reduction, resistances).value
            self.assertEqual(actual, expected)
            alternate = SeriesReduction(graph, lam, gauge=len(reduction.chains) - 1)
            self.assertEqual(integrate_fibre(alternate, resistances).value, actual)

    def test_zero_and_input_limits(self):
        G = theta(1, 3, 3)[0]
        reduction = SeriesReduction(G, [[0, 0, 0] for _ in range(G.nv - 1)])
        result = integrate_fibre(reduction, [1, 2, 3])
        self.assertEqual(result.value, 0)
        self.assertEqual(result.nodes, 0)
        self.assertTrue(result.zero_from_sources)
        with self.assertRaises(ValueError):
            integrate_fibre(reduction, [1, 0, 3])
        reduction = SeriesReduction(G, rational_sources(G, 3, 217))
        with self.assertRaises(CubatureLimitError):
            integrate_fibre(reduction, [1, 2, 3], max_nodes=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
