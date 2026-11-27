"""Independent open-index controls for the universal quartic tensor.

The first control contracts every original color index directly, using
arbitrary noncommuting matrices and an arbitrary non-invariant W tensor.
It does not use a Lie-algebra relation or an input color contraction.
The later controls use all three SU(2) invariant quartic pairings and all
729 entries of the six-index input/output tensor.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
import json
import unittest

import numpy as np

import quartic_tensor as universal


def direct_color_tensor(record, W, generators, structure):
    """Original Wick graph as a complete index contraction, without path tracing."""
    fields = tuple(field for field in record["initial_fields"] if field[0] != "gamma5")
    index = {field: i for i, field in enumerate(fields)}
    args = [W, [index[("gamma", flavor, 0, 0)] for flavor in range(1, 5)]]
    for left, right in record["pairs"]:
        dimension = structure.shape[0] if left[0] == "b" else W.shape[left[1] - 1]
        args += [np.eye(dimension, dtype=np.int64), [index[left], index[right]]]
    factor = Fraction(1)
    for vertex, kind in enumerate(record["cubic_types"], 1):
        if kind == "M":
            flavor = record["matter_flavors"][vertex - 1]
            tensor_fields = (("c", 0, vertex, 0), ("beta", flavor, vertex, 0),
                             ("gamma", flavor, vertex, 0))
            tensor = generators[flavor - 1].transpose(0, 2, 1)
        else:
            tensor_fields = (("b", 0, vertex, 0), ("c", 0, vertex, 0), ("c", 0, vertex, 1))
            tensor = structure
            factor *= Fraction(-1, 2)
        args += [tensor, [index[field] for field in tensor_fields]]
    ghost = next(field for field in record["surviving_fields"] if field[0] == "c")
    gamma = next(field for field in record["surviving_fields"] if field[:2] == ("gamma", 4))
    free = (("beta", 1, 3, 0), ("beta", 2, 3, 0), ("b", 0, 4, 0),
            ("beta", 3, 4, 0), ghost, gamma)
    result = np.einsum(*args, [index[field] for field in free], optimize="greedy")
    return {tuple(map(int, entry)): factor * int(result[entry])
            for entry in zip(*np.nonzero(result))}


def word_color_tensor(word, W, generators, structure):
    """Apply an ordered matrix word directly to each original W tensor index."""
    result = {}
    for external, ghost in product(range(structure.shape[0]), repeat=2):
        if word[0] == "MM":
            actions = []
            for flavor, letters in word[1]:
                matrix = np.eye(W.shape[flavor - 1], dtype=np.int64)
                for letter in letters:
                    matrix = matrix @ generators[flavor - 1][external if letter == "A" else ghost]
                actions.append((flavor, matrix))
        else:
            flavor = word[1]
            matrix = sum((int(structure[color, external, ghost]) * generators[flavor - 1][color]
                          for color in range(structure.shape[0])),
                         np.zeros((W.shape[flavor - 1], W.shape[flavor - 1]), dtype=np.int64))
            actions = [(flavor, matrix)]
        value = W
        for flavor, matrix in actions:
            value = np.moveaxis(np.tensordot(matrix, value, axes=(1, flavor - 1)), 0, flavor - 1)
        for entry in zip(*np.nonzero(value)):
            i1, i2, i3, i4 = map(int, entry)
            result[i1, i2, external, i3, ghost, i4] = Fraction(int(value[entry]))
    return result


def su2_data(pairing=None):
    structure = np.array([[[universal.wick.epsilon(c, a, b) for b in range(3)]
                           for a in range(3)] for c in range(3)], dtype=np.int64)
    generators = tuple(-structure.copy() for flavor in range(4))
    if pairing is None:
        W = np.array([((7 * i + 3) % 5) - 2 for i in range(81)], dtype=np.int64).reshape((3,) * 4)
    else:
        W = np.zeros((3,) * 4, dtype=np.int64)
        for indices in product(range(3), repeat=4):
            W[indices] = int(all(indices[a - 1] == indices[b - 1] for a, b in pairing))
    return W, generators, structure


def assemble_color_tensor(values, W, generators, structure):
    result = [defaultdict(lambda: defaultdict(Fraction)), defaultdict(lambda: defaultdict(Fraction))]
    for word, spinor in values.items():
        color = word_color_tensor(word, W, generators, structure)
        for coordinate in range(2):
            for label, coefficient in spinor[coordinate].items():
                for indices, tensor_entry in color.items():
                    result[coordinate][label][indices] += coefficient * tensor_entry
    return tuple({label: {indices: coefficient for indices, coefficient in entries.items() if coefficient}
                  for label, entries in channels.items() if any(entries.values())}
                 for channels in result)


class UniversalQuarticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geometry = universal.read_geometry()
        cls.groups, cls.records, cls.rejected = universal.grouped_matchings()
        cls.tensor = universal.covariant_tensor(cls.geometry, cls.groups)

    def test_every_open_index_and_every_matching_without_representation_relations(self):
        dimensions = (2, 3, 2, 3)
        W = np.array([((11 * i + 1) % 5) - 2 for i in range(36)], dtype=np.int64).reshape(dimensions)
        generators = tuple(np.array([((3 * i + flavor) % 5) - 2 for i in range(3 * dimension ** 2)],
                                    dtype=np.int64).reshape(3, dimension, dimension)
                           for flavor, dimension in enumerate(dimensions))
        structure = np.zeros((3, 3, 3), dtype=np.int64)
        for c in range(3):
            for a in range(3):
                for b in range(a + 1, 3):
                    structure[c, a, b] = ((c + 2 * a + b) % 5) - 2
                    structure[c, b, a] = -structure[c, a, b]
        for record in self.records:
            word, factor = universal.matching_action(record)
            direct = direct_color_tensor(record, W, generators, structure)
            traced = {indices: factor * coefficient
                      for indices, coefficient in word_color_tensor(word, W, generators, structure).items()
                      if factor * coefficient}
            self.assertEqual(direct, traced, (record["cubic_types"], record["matter_flavors"], record["pairs"]))

    def test_representation_relation_without_W_invariance(self):
        W, generators, structure = su2_data()
        for a, b in product(range(3), repeat=2):
            self.assertTrue(np.array_equal(generators[0][a] @ generators[0][b] - generators[0][b] @ generators[0][a],
                                           sum(structure[c, a, b] * generators[0][c] for c in range(3))))
        sample = self.geometry["independent_full_source_replays"][0]
        values = universal.evaluate_sample(sample, self.groups)
        actual = assemble_color_tensor(values, W, generators, structure)
        self.assertTrue(any(coordinate.get("1") for coordinate in actual))
        self.assertTrue(all(set(coordinate) <= {"1"} for coordinate in actual))
        # The exact arbitrary W used above is deliberately not invariant.
        ward = sum(np.moveaxis(np.tensordot(generators[flavor][0], W, axes=(1, flavor)), 0, flavor)
                   for flavor in range(4))
        self.assertTrue(np.any(ward))

    def test_all_three_invariant_quartic_tensors_with_all_indices_open(self):
        for pairing in (((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))):
            W, generators, structure = su2_data(pairing)
            for record in self.records:
                word, factor = universal.matching_action(record)
                self.assertEqual(direct_color_tensor(record, W, generators, structure),
                                 {indices: factor * coefficient for indices, coefficient
                                  in word_color_tensor(word, W, generators, structure).items() if factor * coefficient})
            for sample in self.geometry["source_polynomial_certificate"]["samples"]:
                values = universal.evaluate_sample(sample, self.groups)
                actual = assemble_color_tensor(values, W, generators, structure)
                self.assertTrue(all(set(coordinate) <= {"1"} for coordinate in actual),
                                (pairing, sample["source_spinors_gamma4_c_gamma5"]))

    def test_complete_old_SU2_input_contraction_and_full_source_replays(self):
        W, generators, structure = su2_data(((1, 3), (2, 4)))
        for sample in self.geometry["independent_full_source_replays"]:
            values = universal.evaluate_sample(sample, self.groups)
            opened = assemble_color_tensor(values, W, generators, structure)
            expected = sample["exact_spinors_by_vertex_type"]["complete"]
            for coordinate in range(2):
                contracted = defaultdict(lambda: defaultdict(Fraction))
                for label, tensor in opened[coordinate].items():
                    for (i1, i2, external, i3, ghost, gamma), coefficient in tensor.items():
                        if i1 == i2 and external == i3:
                            contracted[label][ghost, gamma] += coefficient
                for ghost, gamma in product(range(3), repeat=2):
                    self.assertEqual({label: entries[ghost, gamma] for label, entries in contracted.items()
                                      if entries[ghost, gamma]},
                                     {label: Fraction(coefficient) for label, coefficient in expected[coordinate].items()
                                      if ghost == gamma})

    def test_full_tensor_relation_is_exact_before_any_color_evaluation(self):
        reduced = universal.apply_representation_relation(self.tensor)
        self.assertEqual(len(self.tensor), 12)
        self.assertEqual(len(reduced), 10)
        self.assertTrue(all(set(channels) == {"1"} for channels in reduced.values()))
        for word, channels in self.tensor.items():
            for terms in channels.values():
                self.assertTrue(all(sum(powers) == 3 for vector, powers, coefficient in terms))

    def test_complete_certificate_replay(self):
        actual = universal.base._jsonable(universal.certify())
        expected = json.loads((universal.HERE / "quartic_tensor_certificate.json").read_text())
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
