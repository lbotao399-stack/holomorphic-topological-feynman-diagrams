"""Independent whole-graph oracle for the certified core reconstruction."""
from fractions import Fraction
from itertools import combinations
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(HERE.parent / "cores"))
from htk.series import SeriesReduction
import numer
from core_numerator import certified_core_numerator, integrate_core, primes_for
from run_examples import subdivide, sources


def whole_graph_oracle(graph, source_vectors, reduction):
    result = {}
    for exponents, coefficient in numer.bare_numerator(graph, source_vectors).items():
        powers = []
        for chain in reduction.chains:
            degree = sum(exponents[edge] for edge in chain.edges)
            coefficient *= Fraction(math.prod(math.factorial(exponents[edge]) for edge in chain.edges),
                                    math.factorial(degree + chain.length - 1))
            powers.append(degree + chain.length)
        powers = tuple(powers)
        result[powers] = result.get(powers, Fraction(0)) + coefficient
    return {key: value for key, value in result.items() if value}


def coefficient_dict(polynomial):
    return {tuple(map(int, powers)): Fraction(integer, polynomial.denominator)
            for powers, integer in zip(polynomial.powers, polynomial.integers) if integer}


class CertifiedCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph = subdivide(4, list(combinations(range(4), 2)), [2, 1, 1, 1, 1, 1])
        cls.sources = sources(cls.graph, 2, 317)
        cls.reduction = SeriesReduction(cls.graph, cls.sources)
        cls.expected = whole_graph_oracle(cls.graph, cls.sources, cls.reduction)
        cls.temporary = tempfile.TemporaryDirectory()
        cls.polynomial = certified_core_numerator(cls.reduction, checkpoint_dir=cls.temporary.name,
                                                   threads=2, verbose=False)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_every_coefficient_and_complete_period(self):
        self.assertEqual(coefficient_dict(self.polynomial), self.expected)
        self.assertTrue(self.polynomial.certificate["strict_uniqueness_inequality"])
        import per5
        expected = {}
        for powers, coefficient in self.expected.items():
            for label, value in per5.period(self.polynomial.edges, powers, self.polynomial.q).items():
                expected[label] = expected.get(label, Fraction(0)) + coefficient * value
        expected = {key: value for key, value in expected.items() if value}
        actual, _ = integrate_core(self.polynomial)
        self.assertEqual(actual, expected)
        self.assertTrue(actual)

    def test_rational_sources_and_different_gauge(self):
        rational = [[Fraction(value, i + 2) for i, value in enumerate(row)] for row in self.sources]
        reduction = SeriesReduction(self.graph, rational, gauge=5)
        polynomial = certified_core_numerator(reduction, threads=2, verbose=False)
        self.assertEqual(coefficient_dict(polynomial),
                         {key: value / 216 for key, value in self.expected.items()})

    def test_short_prime_supply_does_not_claim_certification(self):
        prime = next(primes_for(self.reduction))
        with self.assertRaisesRegex(ArithmeticError, "proved CRT bound"):
            certified_core_numerator(self.reduction, primes=[prime],
                                      checkpoint_dir=self.temporary.name, verbose=False)

    def test_checkpoint_replay_and_tampering(self):
        replay = certified_core_numerator(self.reduction, checkpoint_dir=self.temporary.name,
                                            verbose=False)
        self.assertEqual(coefficient_dict(replay), self.expected)
        checkpoint = next(Path(self.temporary.name).glob("*/coefficients_*.npy"))
        saved = np.load(checkpoint, allow_pickle=False)
        changed = saved.copy()
        prime = int(checkpoint.stem.rsplit("_", 1)[1])
        changed[0] = (int(changed[0]) + 1) % prime
        np.save(checkpoint, changed, allow_pickle=False)
        try:
            with self.assertRaisesRegex(ValueError, "checkpoint checksum"):
                certified_core_numerator(self.reduction, checkpoint_dir=self.temporary.name,
                                          verbose=False)
        finally:
            np.save(checkpoint, saved, allow_pickle=False)


if __name__ == "__main__":
    unittest.main(verbosity=2)
