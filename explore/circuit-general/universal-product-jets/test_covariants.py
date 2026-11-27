"""Independent complete-index controls for higher-spinor reconstruction."""
from fractions import Fraction
from itertools import permutations, product
import unittest

import covariants
import jet_tensor as jets


def direct_symmetric_entries(counts, vectors):
    """Definition by all distinct free-index permutations, without binary storage."""
    order = sum(counts)
    labels = tuple(label for label, count in enumerate(counts) for unused in range(count))
    orderings = sorted(set(permutations(labels)))
    result = []
    for second in range(order + 1):
        indices = (0,) * (order - second) + (1,) * second
        total = sum((product_value(ordering, indices, vectors) for ordering in orderings), Fraction(0))
        result.append(total / len(orderings))
    return tuple(result)


def product_value(ordering, indices, vectors):
    result = Fraction(1)
    for label, index in zip(ordering, indices):
        result *= vectors[label][index]
    return result


class CovariantTests(unittest.TestCase):
    def test_every_free_index_against_explicit_symmetrization(self):
        vectors = ((Fraction(2, 3), Fraction(-3, 2)), (Fraction(5, 4), Fraction(7, 3)),
                   (Fraction(-4, 5), Fraction(1, 6)))
        bracket = lambda left, right: left[0] * right[1] - left[1] * right[0]
        weight = bracket(vectors[0], vectors[1]) ** 3
        for order in range(5):
            for first in range(order + 1):
                for second in range(order + 1 - first):
                    counts = (first, second, order - first - second)
                    term = ((counts, (3, 0, 0), Fraction(7, 11)),)
                    expected = tuple(Fraction(7, 11) * weight * entry for entry in direct_symmetric_entries(counts, vectors))
                    self.assertEqual(covariants.evaluate(term, vectors, order), expected)

    def test_complete_canonical_basis_and_normalization_through_rank_six(self):
        for order in range(7):
            terms = []
            for spectator in range(order + 1):
                for second in range(order + 1 - spectator):
                    counts = order - spectator - second, second, spectator
                    for first in range(4):
                        for middle in range(4 - first):
                            powers = first, middle, 3 - first - middle
                            if spectator and first:
                                continue
                            coefficient = Fraction((-1) ** (spectator + second + middle) * (len(terms) + 1), len(terms) + 2)
                            terms.append((counts, powers, coefficient))
            self.assertEqual(len(terms), 2 * (order + 1) * (order + 5))
            polynomials = covariants.normalized_polynomials(terms, order)
            restored = covariants.restore(polynomials, order)
            self.assertEqual(tuple(sorted(terms)), restored)

    def test_immutable_first_current_certificate_as_independent_control(self):
        previous = jets.universal.covariant_tensor(jets.universal.read_geometry())
        for word, labels in previous.items():
            for label, old_terms in labels.items():
                terms = tuple((tuple(int(i == vector) for i in range(3)), powers, coefficient)
                              for vector, powers, coefficient in old_terms)
                restored = covariants.restore(covariants.normalized_polynomials(terms, 1), 1)
                self.assertEqual(tuple(sorted(terms)), restored)
                for vectors in (((1, 2), (3, -1), (-2, 4)), ((1, 2), (2, 4), (3, -1))):
                    self.assertEqual(covariants.evaluate(restored, vectors, 1),
                                     jets.base.evaluate_covariant(old_terms, vectors))

    def test_incomplete_non_covariant_and_inexact_inputs_fail_explicitly(self):
        with self.assertRaises(ValueError):
            covariants.restore([{}], 2)
        with self.assertRaises(ArithmeticError):
            covariants.restore([{(4, 0): Fraction(1)}, {}], 1)
        with self.assertRaises(ArithmeticError):
            covariants.restore([{(5, 0): Fraction(1)}, {}], 1)
        with self.assertRaises(TypeError):
            covariants.evaluate((), ((1.0, 0), (0, 1), (1, 1)), 1)
        with self.assertRaises(RuntimeError):
            covariants.restore([{(0, 0): 1, (1, 0): 1}], 0, max_terms=1)
        with self.assertRaises(RuntimeError):
            jets.calculate((100,), max_samples=1)


if __name__ == "__main__":
    unittest.main()
