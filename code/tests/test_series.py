"""Independent rational checks of series determinant contraction.

Run ``python3 code/tests/test_series.py``.  The oracle solves the full grounded
vertex Laplacian and evaluates the original, uncontracted current determinant.
It does not use chain reduction or the loop solver in ``SeriesReduction``.
"""
import math
import os
import random
import sys
import unittest
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from htk import brute
from htk.graph import Graph, polygon, theta, fan
from htk.series import SeriesReduction, series_chains, _contract_constant_columns


def subdivide(nv, edges, lengths):
    result = []
    for (u, v), length in zip(edges, lengths):
        path = [u] + list(range(nv, nv + length - 1)) + [v]
        nv += length - 1
        result.extend(zip(path, path[1:]))
    return Graph(nv, result)


def complete_edges(nv):
    return [(u, v) for u in range(nv) for v in range(u + 1, nv)]


def rational_sources(G, H, seed):
    rng = random.Random(seed)
    return [[F(rng.randint(-4, 4), rng.randint(1, 3)) for _ in range(H)]
            for _ in range(G.nv - 1)]


def rational_times(G, seed):
    rng = random.Random(seed)
    return [F(rng.randint(1, 11), rng.randint(1, 5)) for _ in G.edges]


def cases():
    yield polygon(4), 3
    yield Graph(2, [(0, 1), (1, 0)]), 1
    yield theta(1, 3, 3)[0], 3
    yield theta(3, 3, 3)[0], 4
    yield fan(4), 2
    yield subdivide(4, complete_edges(4), [3, 2, 1, 1, 2, 1]), 3
    yield Graph(6, [(u, v) for u in range(3) for v in range(3, 6)]), 2
    cube = [(0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]
    yield subdivide(8, cube, [2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1]), 3
    yield subdivide(5, complete_edges(5), [2 if i in (0, 7, 9) else 1 for i in range(10)]), 2


class SeriesTests(unittest.TestCase):
    def test_rank_after_an_early_dependent_column(self):
        rank, multiplier, rows = _contract_constant_columns(
            [(0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0)], 4)
        self.assertEqual(rank, 2)
        self.assertEqual(multiplier, 0)
        self.assertEqual(rows, ())

    def test_currents_and_original_density(self):
        for index, (G, H) in enumerate(cases()):
            lam = rational_sources(G, H, 171 + index)
            t = rational_times(G, 771 + index)
            reference_currents = brute.currents(G, t, lam)
            reference = brute.density(G, t, lam).rational()
            for gauge in {0, len(series_chains(G)) - 1}:
                with self.subTest(loops=G.loops, H=H, gauge=gauge):
                    reduction = SeriesReduction(G, lam, gauge=gauge)
                    self.assertEqual(reduction.currents(t), reference_currents)
                    self.assertEqual(reduction.density(t), reference)
                    self.assertEqual(reduction.determinant_size, len(reduction.chains) - 1)
                    resistances = reduction.core_resistances(t)
                    fractions = reduction.fractions(t)
                    self.assertEqual(reduction.times(resistances, fractions), t)
                    jacobian = math.prod(R ** (chain.length - 1)
                                         for chain, R in zip(reduction.chains, resistances))
                    self.assertEqual(reduction.fibre_density(resistances, fractions),
                                     reference * jacobian)
                    powers = {(0, 0): 2, (G.ne - 1, H - 1): 1}
                    insertion = math.prod(reference_currents[e][i] ** p
                                          for (e, i), p in powers.items())
                    self.assertEqual(reduction.density(t, powers), reference * insertion)
                    self.assertEqual(reduction.fibre_density(resistances, fractions, powers),
                                     reference * insertion * jacobian)

    def test_eight_and_twelve_loop_pointwise_controls(self):
        # A subdivided wheel at H=3 and a once-subdivided wheel at H=2.
        # These check the exact determinant identity, not completed periods.
        controls = [(8, 3, [2] * 9 + [1] * 7, 15),
                    (12, 2, [2] + [1] * 23, 23)]
        for loops, H, lengths, contracted_size in controls:
            edges = ([(0, j) for j in range(1, loops + 1)]
                     + [(j, 1 + j % loops) for j in range(1, loops + 1)])
            G = subdivide(loops + 1, edges, lengths)
            rng = random.Random(108 + loops + H)
            lam = [[rng.randint(-3, 3) for _ in range(H)] for _ in range(G.nv - 1)]
            full = lam + [[-sum(row[i] for row in lam) for i in range(H)]]
            times = [F(rng.randint(1, 5)) for _ in G.edges]
            with self.subTest(loops=loops, H=H):
                self.assertEqual(G.loops, loops)
                reduction = SeriesReduction(G, full)
                self.assertEqual(G.ne - 1, 24)
                self.assertEqual(reduction.determinant_size, contracted_size)
                self.assertEqual(reduction.currents(times), brute.currents(G, times, full))
                value = reduction.density(times)
                self.assertNotEqual(value, 0)
                self.assertEqual(value, brute.density(G, times, full).rational())

    def test_original_edge_orientations_and_orders(self):
        rng = random.Random(891)
        G = subdivide(4, complete_edges(4), [3, 2, 1, 1, 2, 1])
        lam = rational_sources(G, 3, 131)
        t = rational_times(G, 231)
        for _ in range(8):
            order = list(range(G.ne))
            rng.shuffle(order)
            reordered = G.permute_edges(order).reverse_edges(rng.sample(range(G.ne), 4))
            reordered_t = [t[e] for e in order]
            reduction = SeriesReduction(reordered, lam)
            self.assertEqual(reduction.currents(reordered_t), brute.currents(reordered, reordered_t, lam))
            self.assertEqual(reduction.density(reordered_t),
                             brute.density(reordered, reordered_t, lam).rational())

    def test_zero_columns_and_source_conventions(self):
        G = theta(1, 3, 3)[0]
        lam = rational_sources(G, 3, 6)
        reduction = SeriesReduction(G, lam)
        vertex = next(chain.vertices[1] for chain in reduction.chains if chain.length > 1)
        # If the selected vertex is the root, adjust the first non-root source
        # so that the balanced root source vanishes instead.
        if vertex == G.root:
            lam[0] = [-sum(lam[v][i] for v in range(1, G.root)) for i in range(3)]
        else:
            lam[vertex] = [F(0)] * 3
        reduction = SeriesReduction(G, lam)
        t = rational_times(G, 91)
        self.assertTrue(reduction.zero)
        self.assertEqual(reduction.density(t), 0)
        self.assertEqual(brute.density(G, t, lam).rational(), 0)
        zero = SeriesReduction(G, [[0] * 3 for _ in range(G.root)])
        self.assertTrue(zero.zero)
        self.assertEqual(zero.currents(t), [[0] * 3 for _ in G.edges])

        lam = rational_sources(G, 3, 11)
        full = lam + [[-sum(row[i] for row in lam) for i in range(3)]]
        value = SeriesReduction(G, lam).density(t)
        self.assertEqual(SeriesReduction(G, full).density(t), value)
        self.assertEqual(SeriesReduction(G, dict(enumerate(full))).density(t), value)

    def test_polynomial_degree_and_boundary_values(self):
        G = theta(1, 3, 3)[0]
        reduction = SeriesReduction(G, rational_sources(G, 3, 55))
        resistances = [F(j + 1) for j in range(len(reduction.chains))]
        powers = {(0, 0): 2, (G.ne - 1, 2): 1}
        bound = reduction.fibre_degree(powers)
        values = []
        for step in range(bound + 2):
            # A line inside the product of simplices; finite differences above
            # the safe degree vanish exactly.  All sources are nonzero.
            fractions = []
            for chain in reduction.chains:
                row = [F(1, chain.length)] * chain.length
                if chain.length > 1:
                    shift = F(step, 2 * chain.length * (bound + 2))
                    row[0] += shift
                    row[-1] -= shift
                fractions.append(row)
            t = reduction.times(resistances, fractions)
            oracle = brute.density(G, t, reduction.lam).rational()
            y = brute.currents(G, t, reduction.lam)
            oracle *= math.prod(y[e][i] ** p for (e, i), p in powers.items())
            oracle *= math.prod(R ** (c.length - 1) for R, c in zip(resistances, reduction.chains))
            value = reduction.fibre_density(resistances, fractions, powers)
            self.assertEqual(value, oracle)
            values.append(value)
        for _ in range(bound + 1):
            values = [b - a for a, b in zip(values, values[1:])]
        self.assertEqual(values, [0])
        fractions = [tuple([F(1)] + [F(0)] * (chain.length - 1)) for chain in reduction.chains]
        boundary_t = reduction.times(resistances, fractions)
        y = reduction.currents(boundary_t)
        # The loop equations remain meaningful with individual zero times.
        self.assertTrue(brute.check_kirchhoff(G, boundary_t, reduction.lam, y))
        self.assertIsInstance(reduction.fibre_density(resistances, fractions, powers), F)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            series_chains(Graph(4, [(0, 1), (2, 3)]))
        G = polygon(3)
        with self.assertRaises(ValueError):
            SeriesReduction(G, [[1], [2]])
        reduction = SeriesReduction(G, [[1, 0], [0, 1]])
        for t in ([1, 2], [0, 0, 0], [-1, 1, 1]):
            with self.assertRaises(ValueError):
                reduction.density(t)
        for powers in ({(0, 0): -1}, {(3, 0): 1}, {(0, 2): 1}, {(0, 0): F(1, 2)}):
            with self.assertRaises(ValueError):
                reduction.density([1, 1, 1], powers)
        with self.assertRaises(ValueError):
            reduction.times([1], [(F(1, 2), F(1, 2), F(1, 2))])


if __name__ == "__main__":
    unittest.main(verbosity=2)
