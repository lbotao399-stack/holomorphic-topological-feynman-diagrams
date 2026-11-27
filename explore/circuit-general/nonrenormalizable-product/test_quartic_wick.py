"""Independent checks of the complete nonrenormalizable Wick calculation."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
import hashlib
import json
from math import factorial
from pathlib import Path
import unittest

import quartic_wick as calculation


def direct_differential_sum():
    """Apply every Delta to monomials; do not enumerate complete matchings."""
    result = defaultdict(int)
    possible_edges = tuple(edge for edge in combinations(range(5), 2) if edge != (3, 4))
    for physical in combinations(possible_edges, 6):
        edges = physical + ((3, 4),)
        if any(sum(a in vertices and b in vertices for a, b in edges) > 2 * len(vertices) - 3
               for size in range(2, 6) for vertices in combinations(range(5), size)):
            continue
        for first in range(1, 5):
            for second in range(1, 5):
                # This construction is independent of the matching enumerator.
                fields = tuple(("gamma", flavor, 0) for flavor in range(1, 5))
                fields += (("c", 0, 1), ("gamma", first, 1), ("beta", first, 1),
                           ("c", 0, 2), ("gamma", second, 2), ("beta", second, 2),
                           ("beta", 1, 3), ("beta", 2, 3), ("gamma5", 5, 3),
                           ("b", 0, 4), ("beta", 3, 4))
                states = {(fields, None): 1}
                for edge_index in range(5, -1, -1):
                    a, b = physical[edge_index]
                    following = defaultdict(int)
                    for (remaining, derivative_edge), coefficient in states.items():
                        for left in (field for field in remaining if field[2] == a):
                            for right in (field for field in remaining if field[2] == b):
                                species = (left[0], right[0])
                                if species in (("b", "c"), ("c", "b")):
                                    sign = 1
                                elif species in (("beta", "gamma"), ("gamma", "beta")) and left[1] == right[1]:
                                    sign = -1
                                else:
                                    continue
                                after = list(remaining)
                                for field in (right, left):
                                    if field[0] in ("beta", "c"):
                                        sign *= (-1) ** sum(old[0] in ("beta", "c")
                                                           for old in after[:after.index(field)])
                                    after.remove(field)
                                current = edge_index if ("beta", 2, 3) in (left, right) else derivative_edge
                                following[tuple(after), current] += sign * coefficient
                    states = {key: value for key, value in following.items() if value}
                for (remaining, derivative_edge), coefficient in states.items():
                    result[edges, (first, second), remaining, derivative_edge] += coefficient
    return {key: value for key, value in result.items() if value}


class QuarticWickTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads((Path(__file__).parent / "quartic_wick_certificate.json").read_text())

    def test_independent_full_delta_evaluation(self):
        valid, rejected = calculation.wick_matchings()
        self.assertEqual((len(valid), len(rejected)), (20, 28))
        expected = defaultdict(int)
        for record in valid:
            expected[record["ordered_edges"], record["matter_flavors"], record["surviving_fields"],
                     record["current_edge"]] += record["wick_sign"]
        self.assertEqual(direct_differential_sum(), {key: value for key, value in expected.items() if value})

    def test_all_interaction_placements_and_source_retention(self):
        groups, records, rejected = calculation.grouped_matchings()
        for vertex in (1, 2):
            self.assertEqual(calculation.grouped_matchings(vertex)[0], groups)
        self.assertEqual({pair for terms in groups.values() for pair in terms},
                         {(1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (2, 4)})
        for replay in self.certificate["independent_full_source_replays"]:
            for graph in replay["evaluation"]["graphs"]:
                rows = tuple(tuple(map(F, row)) for row in graph["full_balanced_sources"])
                self.assertEqual(tuple(sum(row[i] for row in rows) for i in (0, 1)), (0, 0))
                self.assertEqual(len(rows), 5)

    def test_complete_interpolation_lower_set(self):
        nodes = tuple((a, b) for a in range(5) for b in range(5 - a))
        for i, j in nodes:
            self.assertEqual(calculation.interpolate({(a, b): F(a ** i * b ** j) for a, b in nodes}),
                             {(i, j): F(1)})
        with self.assertRaises(ValueError):
            calculation.interpolate({node: F(0) for node in nodes[:-1]})

    def test_full_certificate_reconstruction(self):
        certificate = self.certificate
        self.assertTrue(certificate["zeta3_cancels_identically"])
        self.assertTrue(certificate["vanishes_identically_at_lambda_c_zero"])
        self.assertEqual(len(certificate["complete_covariant_charge_polynomial"]), 7)
        samples = certificate["source_polynomial_certificate"]["samples"]
        for sample in samples:
            self.assertEqual(set(sample["symbolic_charge_coefficients"]),
                             set(certificate["complete_covariant_charge_polynomial"]))
        for pair, channels in certificate["complete_covariant_charge_polynomial"].items():
            self.assertEqual(set(channels), {"1"})
            for sample in samples:
                vectors = sample["source_spinors_gamma4_c_gamma5"]
                value = calculation.evaluate_covariant(channels["1"], vectors)
                self.assertEqual(value, tuple(F(entry.get("1", 0))
                                              for entry in sample["symbolic_charge_coefficients"][pair]))
                self.assertFalse(any("z3" in entry for entry in sample["symbolic_charge_coefficients"][pair]))
        for relative, expected in certificate["implementation_sha256"].items():
            self.assertEqual(hashlib.sha256((calculation.ROOT / relative).read_bytes()).hexdigest(), expected)

    def test_fresh_gauge_invariant_input_replay(self):
        # Both inputs have charge zero, but no Q-closure is assumed.
        charges = (1, -1, 0, 0)
        self.assertEqual(sum(charges), 0)
        self.assertEqual(sum(charge ** 3 for charge in charges), 0)
        self.assertEqual(-charges[0] - charges[1], 0)
        self.assertEqual(-charges[2], 0)
        direct = calculation.evaluate(((1, 2), (3, -1), (-2, 4)))
        stored = self.certificate["independent_full_source_replays"][3]
        self.assertEqual(calculation.specialize(direct, charges),
                         tuple({label: F(coefficient) for label, coefficient in entry.items()}
                               for entry in stored["complete_raw_spinor"]))
        self.assertTrue(any(calculation.specialize(direct, charges)))

    def test_normalization_and_input_guards(self):
        # Eq. (eq-2-31), I/K at n=3, two matter factors, head derivative.
        self.assertEqual((-1) ** 3 * F(3, factorial(3)) * (-1) * (-2), -1)
        self.assertEqual((-1) * F(-1, 2) ** 7, F(1, 128))
        with self.assertRaises(TypeError):
            calculation._sources((0, 1, 3), ((1.0, 0), (0, 1), (1, 1)))
        with self.assertRaises(ValueError):
            calculation.wick_matchings(3)


if __name__ == "__main__":
    unittest.main()
