"""Independent exact controls for bounded-state fan current/Taylor transfer.

The references are the established ratio integration, the previously
published ordered-simplex evaluator, and closed conservation identities.
Run: python code/tests/test_fan_moments.py
"""

from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import subprocess
import sys
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "code"))
sys.path.insert(0, str(REPO / "explore/iib-bcov"))

from htk import fan, sp_reduce
from htk.fan_moments import fan_current_moment, fan_endpoint_moment, fan_taylor
from scalar_fan_joint import ordered_uniform_moment


def endpoint_sources(n, first, last):
    """The grounded hub source is fixed by conservation in the htk API."""
    return [tuple(first)] + [(0, 0)] * (n - 1) + [tuple(last)]


def displacements_for(n):
    return tuple((F((-1) ** edge * (edge + 2), edge + 3),
                  F(2 * edge + 1, edge + 4)) for edge in range(2 * n + 1))


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), F(0))


class FanMomentChecks(unittest.TestCase):
    def test_optional_flint_import_boundary(self):
        program = """
import sys
from fractions import Fraction

class BlockFlint:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'flint' or fullname.startswith('flint.'):
            raise ImportError('python-flint is deliberately unavailable in this test')

sys.meta_path.insert(0, BlockFlint())
sys.path.insert(0, sys.argv[1])
import htk
result = htk.fan_taylor(1, (1, 0), (0, 1), ((0, 0),) * 3, 2)
assert result.coefficients == (Fraction(-1, 2), Fraction(0), Fraction(0))
try:
    htk.ForestMoments(htk.fan(1), [(1, 0), (0, 1)], (1, 2))
except ImportError as error:
    assert 'python-flint' in str(error)
else:
    raise AssertionError('forest moments must identify their unavailable optional dependency')
"""
        result = subprocess.run([sys.executable, "-c", program, str(REPO / "code")],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_all_original_edge_taylor_against_ratio_integration(self):
        source_pairs = (((2, -1), (-3, 4)),
                        ((F(2, 3), F(-1, 5)), (F(-3, 7), F(4, 9))))
        for n in range(1, 5):
            reduction = sp_reduce(fan(n))
            for first, last in source_pairs:
                sources = endpoint_sources(n, first, last)
                density = reduction.density(sources, check=0)
                displacements = displacements_for(n)
                reference = tuple(reduction.taylor(sources, displacements, 4,
                                                   dens=density, check=0))
                result = fan_taylor(n, first, last, displacements, 4)
                with self.subTest(loops=n, first=first, last=last):
                    self.assertEqual(result.coefficients, reference)
                    self.assertEqual(result.state_bound, comb(7, 3))
                    self.assertLessEqual(result.peak_states, result.state_bound)
                    self.assertTrue(all(isinstance(value, F) for value in result.coefficients))

    def test_mixed_raw_original_current_words(self):
        first, last = (F(2, 3), F(-1, 5)), (F(-3, 7), F(4, 9))
        for n in range(1, 5):
            reduction = sp_reduce(fan(n))
            sources = endpoint_sources(n, first, last)
            density = reduction.density(sources, check=0)
            words = ({}, {(0, 0): 2, (1, 1): 1, (2 * n, 0): 1},
                     {(2 * n - 2, 1): 1, (2 * n - 1, 0): 2},
                     {(0, 1): 1, (2 * n, 1): 2})
            for powers in words:
                with self.subTest(loops=n, powers=powers):
                    reference = reduction.integrate(sources, powers=powers,
                                                    dens=density, check=0)
                    self.assertEqual(fan_current_moment(n, first, last, powers), reference)

    def test_prior_ordered_uniform_identity(self):
        for n in range(1, 9):
            for shift in range(4):
                powers = tuple(((i + shift) % 3, (2 * i + shift) % 4)
                               for i in range(n + 1))
                with self.subTest(loops=n, powers=powers):
                    self.assertEqual(fan_endpoint_moment(powers),
                                     ordered_uniform_moment(powers))
            self.assertEqual(fan_endpoint_moment(((0, 0),) * (n + 1)), 1)

    def test_endpoint_exponents_beyond_one_byte(self):
        # Integrating the largest ordered uniform point gives this closed
        # factorial expression; no packed exponent key is used in the reference.
        for n in (1, 3, 7):
            for exponent in (255, 256, 257):
                powers = ((exponent, 0),) + ((0, 0),) * n
                expected = F(factorial(n + 1) * factorial(exponent + 1),
                             factorial(n + exponent + 1))
                self.assertEqual(fan_endpoint_moment(powers), expected)
        self.assertEqual(fan_endpoint_moment(((256, 0), (0, 0))), F(1, 129))
        self.assertEqual(fan_endpoint_moment(((0, 256), (0, 0))), F(1, 33153))
        self.assertEqual(fan_current_moment(1, (1, 0), (0, 1), {(1, 0): 256}),
                         -F(1, 258))

    def test_gradient_displacement_phase_identity(self):
        n, order = 5, 6
        G = fan(n)
        first, last = (F(3, 2), F(-5, 7)), (F(-2, 3), F(4, 5))
        potentials = tuple((F(vertex + 1, vertex + 2),
                            F((-1) ** vertex * (vertex + 2), vertex + 3))
                           for vertex in range(n + 2))
        gradients = tuple(tuple(potentials[left][j] - potentials[right][j]
                                for j in range(2)) for left, right in G.edges)
        hub = tuple(-first[j] - last[j] for j in range(2))
        constant = -(dot(first, potentials[0]) + dot(last, potentials[n])
                     + dot(hub, potentials[n + 1]))
        bare = fan_current_moment(n, first, last)
        pure = fan_taylor(n, first, last, gradients, order).coefficients
        self.assertEqual(pure, tuple(bare * constant ** degree / factorial(degree)
                                     for degree in range(order + 1)))
        displacements = displacements_for(n)
        shifted = tuple(tuple(displacements[edge][j] + gradients[edge][j] for j in range(2))
                        for edge in range(G.ne))
        original = fan_taylor(n, first, last, displacements, order).coefficients
        actual = fan_taylor(n, first, last, shifted, order).coefficients
        expected = tuple(sum((original[degree - j] * constant ** j / factorial(j)
                              for j in range(degree + 1)), F(0))
                         for degree in range(order + 1))
        self.assertEqual(actual, expected)

    def test_collinear_sources_and_zero_displacements(self):
        n, order = 9, 5
        result = fan_taylor(n, (1, 2), (3, 6), displacements_for(n), order)
        self.assertEqual(result.coefficients, (F(0),) * (order + 1))
        self.assertEqual(result.peak_states, 0)
        first, last = (1, 2), (-3, 4)
        result = fan_taylor(n, first, last, ((0, 0),) * (2 * n + 1), order)
        self.assertEqual(result.coefficients,
                         (fan_current_moment(n, first, last),) + (F(0),) * order)

    def test_budgets_and_exact_input_guards(self):
        with self.assertRaises(MemoryError):
            fan_taylor(1, (1, 0), (0, 1), ((0, 0),) * 3, 3, max_states=19)
        with self.assertRaises(MemoryError):
            fan_current_moment(2, (1, 1), (2, 3), {(0, 0): 1}, max_terms=1)
        for powers in (((-1, 0), (0, 0)), ((True, 0), (0, 0)),
                       ((F(1, 2), 0), (0, 0)), ((0, 0),)):
            with self.subTest(powers=powers), self.assertRaises(ValueError):
                fan_endpoint_moment(powers)
        for powers in ({(0, 0): -1}, {(0, 2): 1}, {(3, 0): 1}, {(0, 0): True}):
            with self.subTest(powers=powers), self.assertRaises(ValueError):
                fan_current_moment(1, (1, 0), (0, 1), powers)
        for kwargs in ({"loops": 0}, {"order": -1}, {"order": True}, {"max_states": 0}):
            arguments = dict(loops=1, first_source=(1, 0), last_source=(0, 1),
                             displacements=((0, 0),) * 3, order=1)
            arguments.update(kwargs)
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fan_taylor(**arguments)
        with self.assertRaises(ValueError):
            fan_taylor(2, (1, 0), (0, 1), ((0, 0),) * 4, 1)
        with self.assertRaises(ValueError):
            fan_taylor(1, (1, 0, 0), (0, 1), ((0, 0),) * 3, 1)
        with self.assertRaises(TypeError):
            fan_taylor(1, (1.0, 0), (0, 1), ((0, 0),) * 3, 1)
        with self.assertRaises(TypeError):
            fan_taylor(1, (1, 0), (0, 1), ((0, 0), (0.5, 0), (0, 0)), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
