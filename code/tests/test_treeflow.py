"""Exact one-tree initialization checks, including the modular numerator clients.

Run: python3 code/tests/test_treeflow.py
Timings: python3 code/tests/test_treeflow.py --benchmark
"""
import itertools
import math
import os
import random
import sys
import time
import unittest
from fractions import Fraction as F
from unittest.mock import patch

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "cores"))
sys.path.insert(0, os.path.join(ROOT, "explore", "circuit-general", "assembly"))

from htk import brute
from htk.graph import Graph, fan, polygon
from htk.linalg import det, solve
from htk.treeflow import tree_flow
import fastker
import numer
import numerator as assembly_numerator


def electrical_from_tree(G, lam, times, tree=None):
    """Independent exact loop correction to the particular current."""
    initial = tree_flow(G, lam, tree)
    H = len(initial[0]) if initial else 0
    cycles = G.cycle_matrix()
    matrix = [[sum(F(a * b) * t for a, b, t in zip(row, col, times))
               for col in cycles] for row in cycles]
    rhs = [[-sum(F(row[e]) * times[e] * initial[e][i] for e in range(G.ne))
            for i in range(H)] for row in cycles]
    loops = solve(matrix, rhs) if cycles else []
    return [[initial[e][i] + sum(F(cycles[r][e]) * loops[r][i]
                                for r in range(len(cycles)))
             for i in range(H)] for e in range(G.ne)]


class TreeFlowChecks(unittest.TestCase):
    def test_arbitrary_sources_and_multigraphs(self):
        rng = random.Random(20261126)
        graphs = [polygon(3), fan(4),
                  Graph(4, [(0, 1), (1, 0), (2, 1), (2, 3), (3, 0), (1, 3)])]
        for G in graphs:
            for G in (G, G.reverse_edges(range(0, G.ne, 2))):
                for H in (1, 2, 4):
                    lam = [[F(rng.randrange(-5, 6), rng.randrange(1, 5))
                            for _ in range(H)] for _ in range(G.nv - 1)]
                    full = lam + [[-sum(row[i] for row in lam) for i in range(H)]]
                    times = [F(rng.randrange(1, 9), rng.randrange(1, 5)) for _ in G.edges]
                    expected = brute.currents(G, times, lam)
                    for source_input in (lam, full, dict(enumerate(full))):
                        initial = tree_flow(G, source_input)
                        B = G.incidence()
                        self.assertEqual([[sum(F(B[e][v]) * initial[e][i]
                                               for e in range(G.ne)) for i in range(H)]
                                          for v in range(G.nv - 1)], lam)
                        self.assertEqual(electrical_from_tree(G, source_input, times), expected)

    def test_every_tree_and_independent_linear_solve(self):
        G = Graph(4, [(0, 1), (1, 0), (2, 1), (2, 3), (3, 0), (1, 3)])
        lam = [[F(1, 2), F(-2, 3)], [F(-7, 5), F(3)], [F(4), F(8, 7)]]
        B = G.incidence()
        times = [F(e + 2, e + 1) for e in range(G.ne)]
        expected = brute.currents(G, times, lam)
        for tree in itertools.combinations(range(G.ne), G.nv - 1):
            tree_incidence = [[B[e][v] for e in tree] for v in range(G.nv - 1)]
            if not det(tree_incidence):
                continue
            linear_solution = solve(tree_incidence, lam)
            initial = tree_flow(G, lam, tree)
            self.assertEqual([initial[e] for e in tree], linear_solution)
            self.assertEqual(electrical_from_tree(G, lam, times, tree), expected)

    def test_invalid_trees_and_sources(self):
        G = Graph(4, [(0, 1), (1, 2), (2, 0), (2, 3)])
        lam = [[1], [2], [-3]]
        for tree in ([0, 1], [0, 0, 3], [0, 1, 4], [0, 1, 2]):
            with self.assertRaises(ValueError):
                tree_flow(G, lam, tree)
        with self.assertRaises(ValueError):
            tree_flow(Graph(3, [(0, 1)]), [[1], [2]])
        with self.assertRaises(ValueError):
            tree_flow(G, [[1], [2], [-3], [1]])
        with self.assertRaises(ValueError):
            tree_flow(G, [[1], [2, 3], [-3]])

    def test_modular_clients_against_exact_density(self):
        rng = random.Random(91)
        graphs = [fan(2), fan(4),
                  Graph(4, [(0, 1), (1, 0), (1, 2), (2, 3), (3, 0)]),
                  Graph(6, [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4),
                            (2, 3), (2, 5), (3, 4), (4, 5)])]
        prime = 1000003
        for G in graphs:
            for G in (G, G.reverse_edges([0, G.ne - 1])):
                for zero in (False, True):
                    lam = [[0 if zero else rng.randrange(-3, 4) for _ in range(2)]
                           for _ in range(G.nv - 1)]
                    full = lam + [[-sum(row[i] for row in lam) for i in range(2)]]
                    # A regression guard: neither client may enumerate trees.
                    with patch.object(numer, "tree_current_polys", side_effect=AssertionError("tree enumeration")):
                        fast = fastker.FastKernel(G.nv, G.edges, lam)
                        assembly = assembly_numerator.BareNumerator(G, full)
                    self.assertTrue(np.array_equal(fast.y0, assembly.y0))
                    grid = np.array([[rng.randrange(0, 5) for _ in range(G.ne - 1)]
                                     for _ in range(4)], dtype=np.int64)
                    values = fast.samples_mod(grid, prime)
                    for point, value in zip(grid, values):
                        times = [F(1)] + [F(int(x) + 1) for x in point]
                        U = brute.symanzik_U(G, times)
                        expected = brute.density(G, times, lam).rational() * U ** (G.loops + 2)
                        residue = expected.numerator * pow(expected.denominator, -1, prime) % prime
                        self.assertEqual(int(value), residue)
                    if zero:
                        self.assertFalse(fast.y0.any())
                        self.assertFalse(values.any())


def benchmark():
    """Compare the removed enumeration with the one-tree construction.

    Timings are local observations.  Graph sizes and the number of candidate
    tree subsets are exact.  This benchmarks initialization, not integration.
    """
    for n in (5, 7, 9):
        G = fan(n)
        lam = [[(3 * v + 1) % 7 - 3, (5 * v + 2) % 11 - 5]
               for v in range(G.nv - 1)]
        full = lam + [[-sum(row[i] for row in lam) for i in range(2)]]
        start = time.perf_counter()
        all_currents = numer.tree_current_polys(G, full)
        cotrees = {cot for currents in all_currents.values()
                   for polynomial in currents for cot in polynomial}
        cotree = min(cotrees)
        previous = np.array([[all_currents[e][i].get(cotree, 0) for i in range(2)]
                             for e in range(G.ne)], dtype=np.int64)
        old_seconds = time.perf_counter() - start
        start = time.perf_counter()
        for _ in range(100):
            direct = np.array(tree_flow(G, lam), dtype=np.int64)
        new_seconds = (time.perf_counter() - start) / 100
        assert np.array_equal(previous, direct)
        print(dict(loops=n, vertices=G.nv, edges=G.ne,
                   candidate_tree_subsets=math.comb(G.ne, G.nv - 1),
                   old_flow_seconds=round(old_seconds, 6),
                   new_flow_seconds=round(new_seconds, 6)), flush=True)
    for n in (12, 20, 32):
        G = fan(n)
        lam = [[(3 * v + 1) % 7 - 3, (5 * v + 2) % 11 - 5]
               for v in range(G.nv - 1)]
        start = time.perf_counter()
        fastker.FastKernel(G.nv, G.edges, lam)
        print(dict(loops=n, vertices=G.nv, edges=G.ne,
                   constructor_seconds=round(time.perf_counter() - start, 6)), flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["--benchmark"]:
        benchmark()
    else:
        unittest.main()
