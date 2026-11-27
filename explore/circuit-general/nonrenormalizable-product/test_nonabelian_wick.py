"""Independent color, Wick, source, and normalization checks for SU(2)."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

import nonabelian_wick as calculation
import quartic_wick as base


def direct_delta_sum():
    """Apply the literal edge differentials to all possible graph monomials."""
    result = defaultdict(int)
    possible = tuple(edge for edge in combinations(range(5), 2) if edge != (3, 4))
    for physical in combinations(possible, 6):
        edges = physical + ((3, 4),)
        if any(sum(a in vertices and b in vertices for a, b in edges) > 2 * len(vertices) - 3
               for size in range(2, 6) for vertices in combinations(range(5), size)):
            continue
        for kinds in product(("M", "G"), repeat=2):
            for flavors in product(*(range(1, 5) if kind == "M" else (0,) for kind in kinds)):
                fields = [("gamma", flavor, 0, 0) for flavor in range(1, 5)]
                for vertex in (1, 2):
                    if kinds[vertex - 1] == "M":
                        flavor = flavors[vertex - 1]
                        fields += [("c", 0, vertex, 0), ("gamma", flavor, vertex, 0), ("beta", flavor, vertex, 0)]
                    else:
                        fields += [("b", 0, vertex, 0), ("c", 0, vertex, 0), ("c", 0, vertex, 1)]
                fields += [("beta", 1, 3, 0), ("beta", 2, 3, 0), ("gamma5", 5, 3, 0),
                           ("b", 0, 4, 0), ("beta", 3, 4, 0)]
                states = {(tuple(fields), None): 1}
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
                                current = edge_index if ("beta", 2, 3, 0) in (left, right) else derivative_edge
                                following[tuple(after), current] += coefficient * sign
                    states = {key: coefficient for key, coefficient in following.items() if coefficient}
                for (remaining, current), coefficient in states.items():
                    result[kinds, flavors, edges, remaining, current] += coefficient
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def independent_color_tensor(record):
    """Contract explicit integer epsilon/delta arrays, with no union-find."""
    fields = tuple(field for field in record["initial_fields"] if field[0] != "gamma5")
    index = {field: position for position, field in enumerate(fields)}
    identity = np.eye(3, dtype=np.int64)
    epsilon = np.zeros((3, 3, 3), dtype=np.int64)
    epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1
    epsilon[2, 1, 0] = epsilon[1, 0, 2] = epsilon[0, 2, 1] = -1
    args = []
    metric_pairs = [(('gamma', 1, 0, 0), ('gamma', 3, 0, 0)),
                    (('gamma', 2, 0, 0), ('gamma', 4, 0, 0)),
                    (('beta', 1, 3, 0), ('beta', 2, 3, 0)),
                    (('b', 0, 4, 0), ('beta', 3, 4, 0))] + list(record["pairs"])
    for left, right in metric_pairs:
        args += [identity, [index[left], index[right]]]
    prefactor = F(1)
    for vertex, kind in enumerate(record["cubic_types"], 1):
        if kind == "M":
            flavor = record["matter_flavors"][vertex - 1]
            slots = (('c', 0, vertex, 0), ('beta', flavor, vertex, 0), ('gamma', flavor, vertex, 0))
        else:
            slots = (('b', 0, vertex, 0), ('c', 0, vertex, 0), ('c', 0, vertex, 1))
            prefactor *= F(-1, 2)
        args += [epsilon, [index[field] for field in slots]]
    ghost = next(field for field in record["surviving_fields"] if field[0] == "c")
    gamma = next(field for field in record["surviving_fields"] if field[:2] == ("gamma", 4))
    # Every summand is 0,+1,-1 and there are at most 3**14 of them;
    # these fixed integer tensor contractions cannot overflow int64.
    tensor = np.einsum(*args, [index[ghost], index[gamma]], optimize="greedy")
    return tuple(tuple(prefactor * int(value) for value in row) for row in tensor)


class NonabelianWickTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = json.loads((Path(__file__).parent / "nonabelian_wick_certificate.json").read_text())

    def test_literal_delta_completeness(self):
        records, rejected = calculation.uncolored_matchings()
        self.assertEqual((len(records), len(rejected)), (28, 44))
        expected = defaultdict(int)
        for record in records:
            expected[record["cubic_types"], record["matter_flavors"], record["ordered_edges"],
                     record["surviving_fields"], record["current_edge"]] += record["wick_sign"]
        self.assertEqual(direct_delta_sum(), {key: value for key, value in expected.items() if value})
        double_gauge = [record for record in rejected if record["cubic_types"] == ("G", "G")]
        self.assertTrue(double_gauge)
        self.assertTrue(all(record["ordered_edges"].count((0, 3)) == 2 for record in double_gauge))

    def test_all_nine_color_entries_independently(self):
        for record in calculation.uncolored_matchings()[0]:
            direct = independent_color_tensor(record)
            self.assertEqual(direct, calculation.complete_color_tensor(record)["complete_output_color_matrix"])
            self.assertTrue(all(direct[a][b] == (direct[0][0] if a == b else 0)
                                for a in range(3) for b in range(3)))

    def test_representation_and_constant_gauge_invariance(self):
        # -i*t_a has entries epsilon(a,i,j) with row j, column i.
        generators = [np.array([[calculation.epsilon(a, i, j) for i in range(3)]
                                for j in range(3)], dtype=np.int64) for a in range(3)]
        for a in range(3):
            self.assertTrue(np.array_equal(generators[a].T, -generators[a]))
            for b in range(3):
                self.assertTrue(np.array_equal(generators[a] @ generators[b] - generators[b] @ generators[a],
                                               sum(calculation.epsilon(a, b, c) * generators[c] for c in range(3))))
                for c in range(3):
                    self.assertEqual(int(np.trace(generators[a] @ (generators[b] @ generators[c]
                                                                  + generators[c] @ generators[b]))), 0)

    def test_complete_source_polynomial_certificate(self):
        certificate = self.certificate
        self.assertTrue(certificate["zeta3_cancels_identically"])
        self.assertTrue(certificate["vanishes_identically_at_lambda_c_zero"])
        covariants = certificate["complete_covariant_spinors_by_vertex_type"]
        self.assertIn("z3", covariants["MM"])
        self.assertIn("z3", covariants["MG"])
        self.assertEqual(set(covariants["complete"]), {"1"})
        for sample in certificate["source_polynomial_certificate"]["samples"]:
            vectors = sample["source_spinors_gamma4_c_gamma5"]
            for kind in ("MM", "MG", "complete"):
                restored = [defaultdict(F), defaultdict(F)]
                for label, terms in covariants[kind].items():
                    for coordinate, value in enumerate(base.evaluate_covariant(terms, vectors)):
                        if value:
                            restored[coordinate][label] = value
                stored = tuple({label: F(value) for label, value in entry.items()}
                               for entry in sample["exact_spinors_by_vertex_type"][kind])
                self.assertEqual(tuple(map(dict, restored)), stored)
        for relative, expected in certificate["implementation_sha256"].items():
            self.assertEqual(hashlib.sha256((base.ROOT / relative).read_bytes()).hexdigest(), expected)

    def test_fresh_full_source_replay(self):
        vectors = ((2, 1), (-3, 2), (1, -1))
        direct = calculation.evaluate(vectors)
        for kind, channels in self.certificate["complete_covariant_spinors_by_vertex_type"].items():
            restored = [defaultdict(F), defaultdict(F)]
            for label, terms in channels.items():
                for coordinate, value in enumerate(base.evaluate_covariant(terms, vectors)):
                    if value:
                        restored[coordinate][label] = value
            self.assertEqual(tuple(map(dict, restored)), direct["exact_spinors_by_vertex_type"][kind])

    def test_compact_complete_coefficient_and_normalization(self):
        terms = self.certificate["complete_covariant_spinors_by_vertex_type"]["complete"]["1"]
        compact = {vector: F(coefficient) for vector, powers, coefficient in terms if powers == [3, 0, 0]}
        self.assertEqual(compact, {0: F(1, 90), 1: F(1, 144)})
        self.assertEqual((-1) ** 3 * F(1, 2) * (-2), 1)
        self.assertEqual(F(-1, 2) ** 7 / 720, F(-1, 92160))


if __name__ == "__main__":
    unittest.main()
