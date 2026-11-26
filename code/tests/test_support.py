"""Independent exact and modular checks of the numerator support bounds.

Run: python3 code/tests/test_support.py
The unrestricted interpolation and the existing Fraction numerator supply
independent coefficient references; no period table is used.
"""
import itertools
import os
import random
import sys
import types
import unittest
from unittest.mock import patch

import numpy as np
import sympy as sp

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "cores"))
sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "assembly"))

from htk.graph import Graph, fan, polygon
from htk import kernel_exact
from htk.support import cycle_column_rank, numerator_support
from numerator import BareNumerator, LowerSet
import numer


def nonzero_sources(G, H):
    rng = random.Random(47 + 13 * G.ne + H)
    choices = tuple(range(-7, 0)) + tuple(range(1, 8))
    rows = [[rng.choice(choices) for _ in range(H)] for _ in range(G.nv - 1)]
    return rows + [[-sum(row[i] for row in rows) for i in range(H)]]


def subdivide_complete_four(lengths):
    edges = []
    vertex = 4
    for (u, v), length in zip(itertools.combinations(range(4), 2), lengths):
        previous = u
        for _ in range(length - 1):
            edges.append((previous, vertex))
            previous = vertex
            vertex += 1
        edges.append((previous, v))
    return Graph(vertex, edges)


def coefficient_dict(numerator, prime):
    coefficients = numerator.coefficients_mod(prime)
    return {tuple(int(x) for x in exponent): int(value)
            for exponent, value in zip(numerator.exponents(), coefficients) if value}


class SupportChecks(unittest.TestCase):
    def test_graph_rank_against_exact_cycle_columns(self):
        graphs = [polygon(3), fan(2),
                  Graph(4, [(0, 1), (1, 0), (1, 2), (2, 3), (3, 0), (1, 3)]),
                  Graph(4, [(0, 1), (1, 2), (2, 0), (2, 3)])]
        for G in graphs:
            C = sp.Matrix(G.cycle_matrix())
            for mask in range(1 << G.ne):
                selected = tuple(e for e in range(G.ne) if mask & (1 << e))
                self.assertEqual(cycle_column_rank(G, selected),
                                 C.extract(range(G.loops), selected).rank())
        with self.assertRaises(ValueError):
            cycle_column_rank(graphs[0], [graphs[0].ne])
        with self.assertRaises(ValueError):
            cycle_column_rank(Graph(3, [(0, 1)]), [0])

    def test_flats_match_linear_algebra_and_budgets(self):
        graphs = [fan(3), Graph(4, [(0, 1), (1, 0), (1, 2), (2, 3), (3, 0), (1, 3)])]
        for G in graphs:
            C = sp.Matrix(G.cycle_matrix())
            for gauge in (0, G.ne - 1):
                ground = [e for e in range(G.ne) if e != gauge]
                ranks = {}
                for bits in itertools.product((False, True), repeat=len(ground)):
                    subset = frozenset(e for e, present in zip(ground, bits) if present)
                    ranks[subset] = C.extract(range(G.loops), sorted(subset)).rank()
                flats = {subset: rank for subset, rank in ranks.items()
                         if all(ranks[subset | {e}] > rank for e in ground if e not in subset)}
                expected = {tuple(sorted(subset)): G.loops * rank - 1
                            for subset, rank in flats.items() if subset
                            and G.loops * rank - 1 < min((G.loops - 1) * len(subset),
                                                        G.loops * G.loops - 1)}
                full = numerator_support(G, gauge, max_flats=1024)
                self.assertTrue(full.complete)
                self.assertFalse(full.zero)
                self.assertEqual(full.flat_count, len(flats))
                self.assertEqual(dict(full.constraints), expected)
                for budget in (1, 2, 3, len(flats)):
                    partial = numerator_support(G, gauge, max_flats=budget)
                    self.assertLessEqual(partial.flat_count, budget)
                    self.assertEqual(partial.complete, budget >= len(flats))
                    self.assertTrue(set(partial.constraints).issubset(set(full.constraints)))
        bridge = Graph(4, [(0, 1), (1, 2), (2, 0), (2, 3)])
        for gauge in range(bridge.ne):
            self.assertTrue(numerator_support(bridge, gauge).zero)

    def test_grid_order_downwardness_and_newton_transform(self):
        prime = 1000003
        for constraints in ((), (((0, 2), 2), ((1, 2, 3), 3))):
            grid = LowerSet(4, 3, 6, constraints=constraints)
            expected = [exponent for exponent in itertools.product(range(4), repeat=4)
                        if sum(exponent) <= 6
                        and all(sum(exponent[j] for j in selected) <= bound
                                for selected, bound in constraints)]
            expected.sort(key=lambda exponent: sum(x * 4 ** j for j, x in enumerate(exponent)))
            self.assertEqual([tuple(row) for row in grid.dig], expected)
            present = set(expected)
            for exponent in expected:
                for j, value in enumerate(exponent):
                    if value:
                        lower = list(exponent)
                        lower[j] -= 1
                        self.assertIn(tuple(lower), present)
            coefficients = [(17 * i * i + 11 * i + 3) % prime for i in range(grid.M)]
            samples = []
            for point in grid.dig:
                value = 0
                for exponent, coefficient in zip(expected, coefficients):
                    term = coefficient
                    for x, power in zip(point, exponent):
                        term = term * pow(int(x) + 1, power, prime) % prime
                    value = (value + term) % prime
                samples.append(value)
            result = grid.newton_to_monomial(np.array(samples, dtype=np.int64), prime)
            self.assertEqual(list(result), coefficients)

    def test_bounded_and_unrestricted_numerators(self):
        graphs = [subdivide_complete_four((3, 2, 1, 1, 2, 1)),
                  Graph(6, [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5),
                            (2, 3), (2, 5), (3, 4), (4, 5)])]
        for G, H in zip(graphs, (3, 2)):
            sources = nonzero_sources(G, H)
            old = BareNumerator(G, sources, rank_bounds=False)
            new = BareNumerator(G, sources)
            partial = BareNumerator(G, sources, max_flats=2)
            self.assertLessEqual(new.lowerset.M, partial.lowerset.M)
            self.assertLessEqual(partial.lowerset.M, old.lowerset.M)
            self.assertLess(new.lowerset.M, old.lowerset.M)
            for prime in (2147483629, 2147483587):
                old_coefficients = coefficient_dict(old, prime)
                self.assertTrue(old_coefficients)
                self.assertEqual(coefficient_dict(new, prime), old_coefficients)
                self.assertEqual(coefficient_dict(partial, prime), old_coefficients)
            print("support interpolation:", dict(loops=G.loops, H=H,
                                                   unrestricted=old.lowerset.M,
                                                   bounded=new.lowerset.M,
                                                   nonzero=len(old_coefficients)), flush=True)

    def test_every_subset_on_old_fraction_polynomials(self):
        for G in (fan(2), subdivide_complete_four((1, 1, 1, 2, 1, 1))):
            polynomial = numer.bare_numerator(G, nonzero_sources(G, 2))
            self.assertTrue(polynomial)
            C = sp.Matrix(G.cycle_matrix())
            for mask in range(1, 1 << G.ne):
                selected = [e for e in range(G.ne) if mask & (1 << e)]
                bound = G.loops * C.extract(range(G.loops), selected).rank() - 1
                self.assertTrue(all(sum(exponent[e] for e in selected) <= bound
                                    for exponent, coefficient in polynomial.items() if coefficient))

    def test_singular_samples_fail_explicitly(self):
        G = polygon(3)
        numerator = BareNumerator(G, nonzero_sources(G, 2))
        for point in ([[1, 2, 4]], [[7, 1, 1]]):
            for sampler in (numerator.samples_mod, numerator.samples_fast):
                with self.assertRaises(ArithmeticError):
                    sampler(np.array(point, dtype=np.int64), 7)
        # Force the accelerated path even in environments without numba.
        fake = types.SimpleNamespace(samples=lambda times, *args: np.zeros(len(times), dtype=np.int64))
        with patch.dict(sys.modules, {"fastsample": fake}):
            with self.assertRaises(ArithmeticError):
                numerator.samples_fast(np.array([[7, 1, 1]], dtype=np.int64), 7)

    def test_five_loop_client_two_loop_control(self):
        # Exercises the actual rank-grid client and collision-free row grouping.
        sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "fiveloop"))
        import fk5
        from htk.graph import theta
        G = theta(1, 2, 2)[0]
        for sources in ([[2, -3], [5, 7], [-4, 11]], [[0, 0]] * (G.nv - 1)):
            edges, q, coefficients = fk5.core_numerator(G.nv, G.edges, sources, verbose=False)
            direct = kernel_exact(G, sources)
            expected = {"1": direct} if direct else {}
            self.assertEqual(fk5.kernel_from_core(edges, q, coefficients), expected)
            self.assertEqual(fk5.kernel(G.nv, G.edges, sources, verbose=False), expected)
        bridge = Graph(4, [(0, 1), (1, 2), (2, 0), (2, 3)])
        zero = BareNumerator(bridge, nonzero_sources(bridge, 3))
        self.assertEqual(zero.lowerset.M, 0)
        self.assertEqual(coefficient_dict(zero, 2147483629), {})

    def test_memory_guards_and_empty_sets(self):
        with self.assertRaises(MemoryError):
            LowerSet(8, 4, 24, max_points=20)
        with self.assertRaises(ValueError):
            LowerSet(64, 1, 64)
        sparse = LowerSet(12, 3, 0)
        self.assertIsNone(sparse._pos)
        self.assertEqual(sparse.M, 1)
        self.assertEqual(LowerSet(1, 130, 130).dig.max(), 130)
        self.assertEqual(LowerSet(0, 0, 0).dig.shape, (1, 0))
        self.assertEqual(LowerSet(12, 3, 24, constraints=(((), -1),), max_points=1).M, 0)


if __name__ == "__main__":
    unittest.main()
