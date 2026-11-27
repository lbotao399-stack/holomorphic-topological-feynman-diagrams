"""Complete direct-period tensors replayed by rational boundary integration.

The oracle certificates were obtained from the original six graphs.  The
new route is tested with the individual-period worker disabled, including
every source interpolation node and every free symmetric spinor entry.
"""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import covariants
import jet_tensor as jets
from boundary_faces import SourceTriangleTimeFaces, classify_faces
from rational_product import RationalJetProduct
from outer_faces import SourceTriangleOuterFaces
import htk.forest_periods as period_backend

HERE = Path(__file__).resolve().parent
GENERIC = ((1, 2), (3, -1), (-2, 4))


def _tuple(value):
    return tuple(_tuple(item) for item in value) if isinstance(value, list) else value


def _oracle(order):
    data = json.loads((HERE / ("jet_order_%d_certificate.json" % order)).read_text())
    tensor = {}
    for row in data["complete_reduced_tensor"]:
        channels = row["complete_covariant_tensor"]
        if set(channels) - {"1"}:
            raise ArithmeticError("the immutable complete tensor is not rational")
        tensor[_tuple(row["generator_word"])] = tuple(
            (tuple(counts), tuple(powers), Fraction(value)) for counts, powers, value in channels.get("1", ()))
    return data, tensor


class CompleteBoundaryProductTests(unittest.TestCase):
    def test_entire_rank_two_and_three_source_grids_without_period_lookup(self):
        oracles = {order: _oracle(order) for order in (2, 3)}
        source_orders = {}
        for order, (data, tensor) in oracles.items():
            self.assertEqual(len(tensor), 10)
            self.assertEqual(len(data["complete_source_samples"]), (order + 4) * (order + 5) // 2)
            for sample in data["complete_source_samples"] + data["independent_full_source_replays"]:
                sources = tuple(tuple(map(Fraction, row)) for row in sample["source_spinors_gamma4_c_gamma5"])
                source_orders.setdefault(sources, set()).add(order)
        with patch.object(period_backend, "_PeriodWorker", side_effect=AssertionError("an individual period was requested")):
            for sources, orders in sorted(source_orders.items()):
                with self.subTest(sources=sources), RationalJetProduct(sources, check=0) as evaluator:
                    for order in sorted(orders):
                        actual = evaluator.tensor({"beta2": order})
                        self.assertEqual(actual["complete_symmetric_entries_per_color_word"], order + 1)
                        self.assertEqual(set(actual["complete_open_color_tensor"]), set(oracles[order][1]))
                        for word, terms in oracles[order][1].items():
                            expected = covariants.evaluate(terms, sources, order)
                            entries = actual["complete_open_color_tensor"][word]
                            self.assertEqual(len(entries), order + 1)
                            self.assertEqual(tuple(entries[(0, second, 0, 0, 0)] for second in range(order + 1)), expected)
                    self.assertEqual(evaluator.last_stats["individual_period_lookups"], 0)

    def test_every_mixed_rank_four_entry_against_original_graphs(self):
        data = json.loads((HERE.parent / "forest-jet-identities/mixed_input_control.json").read_text())
        sources = tuple(tuple(map(Fraction, row)) for row in data["source_spinors_gamma4_c_gamma5"])
        expected = [{tuple(row["ordered_free_spinor_entries"]): Fraction(row["complete_rational_value"])
                     for row in pair["complete_rank_four_tensor"]} for pair in data["complete_rational_pair_tensors"]]
        with patch.object(period_backend, "_PeriodWorker", side_effect=AssertionError("an individual period was requested")):
            with RationalJetProduct(sources, check=1) as evaluator:
                for indices in product((0, 1), repeat=4):
                    powers = {(name, index): 1 for name, index in zip(("beta1", "beta2", "beta3", "b"), indices)}
                    actual = evaluator.pair_moment(evaluator._word(1, powers))
                    self.assertEqual(actual, tuple(pair[indices] for pair in expected))

    def test_free_spectator_is_a_physical_derivative_after_conversion(self):
        with patch.object(period_backend, "_PeriodWorker", side_effect=AssertionError("an individual period was requested")):
            with RationalJetProduct(GENERIC, check=0) as evaluator:
                base = evaluator.moment({("beta2", 0): 1, ("beta3", 1): 1})
                actual = evaluator.moment({("beta2", 0): 1, ("beta3", 1): 1,
                                           ("gamma5", 0): 2, ("gamma5", 1): 1})
                factor = Fraction(GENERIC[2][0] ** 2 * GENERIC[2][1])
                self.assertEqual(actual, {word: factor * value for word, value in base.items()})
                tensor = evaluator.tensor({"gamma5": 1})
                self.assertIn("after source conversion", tensor["spectator_derivatives"])
        # A raw physical spectator derivative is -2 lambda_gamma5,
        # not a propagator-head factor. Its own source conversion cancels
        # that -2, independently of the contracted derivative count.
        for contracted in range(9):
            for spectator in range(5):
                self.assertEqual(Fraction(1, 2) * (-2) ** contracted *
                                 Fraction(-1, 2) ** (6 + contracted), Fraction(1, 128))
                self.assertEqual((-2) ** spectator * Fraction(-1, 2) ** spectator, 1)

    def test_all_sixty_two_faces_have_exact_classification(self):
        rows = classify_faces()
        self.assertEqual(len(rows), 62)
        self.assertEqual(len({row["selected_edges"] for row in rows}), 62)
        active = [row for row in rows if row["kind"] == "polynomial_SP_chart"]
        self.assertEqual(sorted(len(row["selected_edges"]) for row in active), [1] * 6 + [3] * 4 + [5] * 6)
        for row in rows:
            if row["kind"] == "identically_zero":
                proof = row["exact_source_independent_rank_proof"]
                self.assertTrue(proof["relative_projective_direction_is_nontrivial"])
                self.assertEqual(set(proof["all_leading_weight_scale_exponents"]), {proof["common_weight_scale_exponent"]})
                self.assertEqual(len(proof["all_leading_weight_scale_exponents"]), len(row["all_leading_spanning_tree_indices"]))

    def test_strict_complete_output_cache_and_power_budgets(self):
        with RationalJetProduct(GENERIC, check=0, max_entries=9) as evaluator:
            self.assertEqual(len(evaluator.moment()), 10)
            self.assertEqual(evaluator._cache, {})
            with self.assertRaisesRegex(RuntimeError, "no partial tensor"):
                evaluator.tensor({})
        with RationalJetProduct(GENERIC, check=0, max_entries=15) as evaluator:
            evaluator.moment()
            evaluator.moment({("beta2", 0): 1})
            self.assertEqual(len(evaluator._cache), 1)
            with self.assertRaisesRegex(RuntimeError, "no partial tensor"):
                evaluator.tensor({"beta2": 1000000000})
        for cls in (SourceTriangleTimeFaces, SourceTriangleOuterFaces):
            evaluator = cls(GENERIC, "TWA")
            evaluator.max_terms = 16
            with self.assertRaisesRegex(RuntimeError, "max_terms"):
                evaluator.integrate({(1000000000, 0, 0, 0, 0, 0): Fraction(1)})

    def test_invalid_or_inexact_inputs_fail_before_partial_output(self):
        for value in (1.0, 1j, True):
            with self.assertRaises(TypeError):
                RationalJetProduct(((value, 2), (3, -1), (-2, 4)))
        with RationalJetProduct(GENERIC, check=0) as evaluator:
            for powers in ({("beta2", 0): -1}, {("beta2", 2): 1}, {("unknown", 0): 1}):
                with self.assertRaises(ValueError):
                    evaluator.moment(powers)
            with self.assertRaises(TypeError):
                evaluator.moment({("beta2", 0): True})
            with self.assertRaises(TypeError):
                evaluator.pair_moment({(0,) * 6: 1.0})
            with self.assertRaises(ValueError):
                evaluator.tensor({"unknown": 1})


if __name__ == "__main__":
    unittest.main()
