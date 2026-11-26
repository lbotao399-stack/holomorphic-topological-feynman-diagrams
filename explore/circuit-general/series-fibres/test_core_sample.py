"""Independent exact validation of the compiled reciprocal core sampler.

The expected values use rational chain integration (Grundmann--Moeller),
full original-graph U, and the stored all-ordinate six-loop certificate.
They do not use the finite-field Gauss nodes or the compiled determinant.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "code" / "tests"))

from htk import brute
from htk.fibres import integrate_fibre
from htk.graph import Graph, theta, fan
from htk.series import SeriesReduction
from test_series import subdivide, complete_edges, rational_sources
from core_sample import CoreSampler, SingularSampleError, splitting_prime, _gauss, _mod


def rational_reciprocal(reduction, conductances):
    resistances = [1 / F(value) for value in conductances]
    fibre = integrate_fibre(reduction, resistances).value
    fractions = [[F(1, chain.length)] * chain.length for chain in reduction.chains]
    U = brute.symanzik_U(reduction.G, reduction.times(resistances, fractions))
    multiplier = F(1)
    for chain, R, g in zip(reduction.chains, resistances, conductances):
        multiplier *= F(g) ** (reduction.n - 1) / R ** (chain.length - 1)
    return fibre * U ** (reduction.n + reduction.H) * multiplier


class CoreSamplerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads((HERE / "examples.json").read_text())["examples"][1]
        row = cls.fixture
        cls.G = Graph(row["vertices"], row["original_oriented_edges"])
        cls.reduction = SeriesReduction(cls.G, row["full_balanced_sources"])

    def test_gauss_moment_identities(self):
        # This tests the exact integration rule independently of all graph
        # determinants: integral_0^1 x^k dx = 1/(k+1).
        for loops in range(1, 7):
            prime = next(p for p in (2147483629, 2147483587, 2147483579, 2147483563)
                         if splitting_prime(p, loops))
            nodes = _gauss(loops, prime)
            self.assertEqual(len(nodes), (loops + 1) // 2)
            for degree in range(loops):
                with self.subTest(loops=loops, degree=degree):
                    moment = sum(weight * pow(node, degree, prime)
                                 for node, weight in nodes) % prime
                    self.assertEqual(moment, pow(degree + 1, -1, prime))

    def test_six_loop_stored_certificate_and_nonunit_gauge(self):
        row = self.fixture
        resistances = [F(x["numerator"], x["denominator"]) for x in row["retained_core_resistances"]]
        conductances = [1 / x for x in resistances]
        fibre = F(row["integrated_fibre_value"]["numerator"], row["integrated_fibre_value"]["denominator"])
        U = F(row["symanzik_U"]["numerator"], row["symanzik_U"]["denominator"])
        expected = fibre * U ** 8
        for chain, R, g in zip(self.reduction.chains, resistances, conductances):
            expected *= g ** 5 / R ** (chain.length - 1)
        degree = len(self.reduction.chains) * 5 - 35
        points = [conductances, [3 * x for x in conductances]]
        sampler = CoreSampler(self.reduction, threads=2)
        self.assertEqual(sampler.node_count, 27)
        for prime in (2147483629, 2147483587, 2147483579):
            values = sampler.samples(points, prime)
            self.assertEqual(list(values), [_mod(expected, prime), _mod(3 ** degree * expected, prime)])
        alternate = SeriesReduction(self.G, row["full_balanced_sources"], gauge=9)
        self.assertEqual(list(CoreSampler(alternate).samples(points, 2147483629)),
                         list(sampler.samples(points, 2147483629)))

    def test_one_two_and_three_node_rules(self):
        cases = [(Graph(2, [(0, 1), (1, 0)]), 1),
                 (theta(1, 2, 2)[0], 2),
                 (subdivide(4, complete_edges(4), [2, 1, 1, 1, 1, 1]), 2)]
        for G, H in cases:
            lam = rational_sources(G, H, 131 + G.loops)
            reduction = SeriesReduction(G, lam)
            conductances = [F(2 + i % 4, 1 + i % 2) for i in range(len(reduction.chains))]
            expected = rational_reciprocal(reduction, conductances)
            sampler = CoreSampler(reduction)
            prime = next(p for p in (2147483629, 2147483587, 2147483579, 2147483563)
                         if splitting_prime(p, G.loops))
            self.assertEqual(int(sampler.samples([conductances], prime)[0]), _mod(expected, prime))

    def test_edge_orientation_order_and_empty_fibre(self):
        G = subdivide(4, complete_edges(4), [2, 1, 1, 1, 1, 1])
        G = G.permute_edges([3, 1, 6, 0, 5, 2, 4]).reverse_edges([0, 3, 5])
        reduction = SeriesReduction(G, rational_sources(G, 2, 317))
        point = [i + 1 for i in range(len(reduction.chains))]
        expected = rational_reciprocal(reduction, point)
        prime = next(p for p in (2147483629, 2147483587, 2147483579) if splitting_prime(p, 3))
        self.assertEqual(int(CoreSampler(reduction).samples([point], prime)[0]), _mod(expected, prime))

        G = Graph(6, [(u, v) for u in range(3) for v in range(3, 6)])
        reduction = SeriesReduction(G, rational_sources(G, 2, 451))
        point = [1 + i % 3 for i in range(G.ne)]
        expected = rational_reciprocal(reduction, point)
        # No Gauss nodes are needed, so no quadratic-residue condition applies.
        prime = 101
        self.assertEqual(int(CoreSampler(reduction).samples([point], prime)[0]), _mod(expected, prime))

    def test_singular_and_invalid_data(self):
        sampler = CoreSampler(self.reduction)
        good = [1] * len(self.reduction.chains)
        bad = list(good)
        bad[-1] = F(-3, 2)  # K5 Kirchhoff polynomial = 75 + 50*g_last.
        with self.assertRaises(SingularSampleError) as caught:
            sampler.samples([good, bad], 2147483629)
        self.assertEqual(caught.exception.indices, (1,))
        for prime in (0, 9, 2147483549, 2 ** 31):
            with self.assertRaises(ValueError):
                sampler.samples([good], prime)
        with self.assertRaises(TypeError):
            sampler.samples([[float(x) for x in good]], 2147483629)
        with self.assertRaises(ValueError):
            sampler.samples([[0] + good[1:]], 2147483629)
        with self.assertRaises(ValueError):
            CoreSampler(SeriesReduction(theta(1, 3, 3)[0], [[1, 2, 3]] * 5))
        with self.assertRaises(ValueError):
            CoreSampler(SeriesReduction(fan(7), [[1, 2]] * (fan(7).nv - 1)))
        self.assertEqual(sampler.samples(np.empty((0, len(good)), dtype=np.int64), 2147483629).shape, (0,))

    def test_rational_denominator_and_zero_sources(self):
        prime = 2147483629
        lam = [list(row) for row in self.fixture["full_balanced_sources"][:-1]]
        lam[0][0] = F(1, prime)
        reduction = SeriesReduction(self.G, lam)
        with self.assertRaises(ArithmeticError):
            CoreSampler(reduction).samples([[1] * 10], prime)
        zero = SeriesReduction(self.G, [[0, 0] for _ in range(self.G.nv - 1)])
        self.assertEqual(list(CoreSampler(zero).samples([[1] * 10, [2] * 10], prime)), [0, 0])


if __name__ == "__main__":
    unittest.main(verbosity=2)
