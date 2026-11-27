"""Exact CLI controls against elementary simplex integrals and full tensors."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from solve_product_forest_cubes import solve, write_complete


TRIANGLE = {"vertices": 3, "edges": [[0, 1], [0, 2], [1, 2]],
            "marked": [1, 2], "sources": [[1, 0], [0, 1], [-1, -1]]}
ZERO = {"vertices": 2, "edges": [[1, 0]], "marked": [0, 1],
        "sources": [[2, 3], [-2, -3]]}


def check_row(test, row, expected):
    partial, tail = F(row["partial_sum"]), F(row["absolute_tail_bound"])
    test.assertLessEqual(F(row["lower"]), expected)
    test.assertGreaterEqual(F(row["upper"]), expected)
    test.assertLessEqual(abs(expected - partial), tail)


class ProductForestCubeCliChecks(unittest.TestCase):
    def test_both_backends_enclose_independent_simplex_integral(self):
        document = dict(TRIANGLE, requests=[{"name": "mixed", "polynomial": [
            {"coefficient": 1},
            {"coefficient": 1, "powers": [[0, 0, 1]]},
            {"coefficient": -2, "powers": [[0, 0, 1], [1, 1, 1]]}]}])
        # On the original unit simplex the complete density is -1.
        # The three integrals are -1/2, -1/6 and -1/24, respectively.
        expected = -F(1, 2) - F(1, 6) + 2 * F(1, 24)
        for backend in ("midpoint", "geometric"):
            with self.subTest(backend=backend), patch(
                    "htk.forest_periods.CoreCurrentPeriods._source_numerator",
                    side_effect=AssertionError("a global source numerator was built")):
                result = solve(document, backend=backend, series_order=2)
            self.assertTrue(result["all_requests_complete"])
            self.assertEqual(result["full_balanced_sources"], [["1", "0"], ["0", "1"], ["-1", "-1"]])
            value = result["results"][0]["result"]
            self.assertEqual(value["type"], "complete_convergent_rational_series")
            self.assertEqual(len(value["original_current_polynomial"]), 3)
            self.assertTrue(value["all_original_current_terms_retained"])
            running = F(0)
            for row in value["certified_finite_terms"]:
                running += F(row["term"])
                self.assertEqual(running, F(row["partial_sum"]))
                check_row(self, row, expected)

    def test_whole_tensor_with_repeated_original_edge_and_exact_enclosure(self):
        document = dict(ZERO, requests=[{"name": "full_tensor", "tensor_edges": [0, 0],
                                       "indices": ["+dot_a", "+dot_b"]},
                                      {"name": "mixed", "polynomial": [
                                          {"coefficient": "5/7"},
                                          {"coefficient": "-2/3", "powers": [[0, 0, 2], [0, 1, 1]]}]}])
        for backend in ("midpoint", "geometric"):
            with self.subTest(backend=backend):
                result = solve(document, backend=backend, series_order=2, width="1/100")
                tensor = result["results"][0]["result"]
                self.assertEqual(tensor["free_lower_spinor_indices"], ["+dot_a", "+dot_b"])
                entries = {tuple(row["storage_indices"]): row["value"] for row in tensor["entries"]}
                self.assertEqual(set(entries), {(0, 0), (0, 1), (1, 0), (1, 1)})
                for coordinates, expected in {(0, 0): 4, (0, 1): 6, (1, 0): 6, (1, 1): 9}.items():
                    value = entries[coordinates]
                    self.assertEqual([F(row["term"]) for row in value["certified_finite_terms"]],
                                     [F(expected), F(0), F(0)])
                    enclosure = value["certified_enclosure"]
                    self.assertEqual(F(enclosure["lower"]), expected)
                    self.assertEqual(F(enclosure["upper"]), expected)
                mixed = result["results"][1]["result"]["certified_enclosure"]
                self.assertEqual(F(mixed["lower"]), F(61, 7))
                self.assertEqual(F(mixed["upper"]), F(61, 7))

    def test_specification_without_finite_evaluation_and_balanced_last_source(self):
        document = dict(TRIANGLE, sources=[["2/3", "-5/7"], ["11/13", "17/19"]],
                        requests=[{"name": "word", "powers": [[0, 0, 2], [0, 0, 3], [2, 1, 1]]}])
        for backend in ("midpoint", "geometric"):
            with self.subTest(backend=backend):
                result = solve(document, backend=backend)
                self.assertEqual(result["graph"], {key: TRIANGLE[key] for key in ("vertices", "edges", "marked")})
                self.assertEqual(result["full_balanced_sources"][-1], [str(-F(2, 3)-F(11, 13)), str(F(5, 7)-F(17, 19))])
                value = result["results"][0]["result"]
                self.assertNotIn("certified_finite_terms", value)
                self.assertNotIn("certified_enclosure", value)
                self.assertEqual(value["original_current_polynomial"],
                                 [{"exponents": [5, 0, 0, 0, 0, 1], "coefficient": "1"}])
                self.assertTrue(value["series_specification"])
                json.dumps(result)

    def test_original_edge_reversal_changes_every_odd_current(self):
        original = dict(ZERO, requests=[{"powers": [[0, 0, 1]]}])
        reversed_graph = dict(original, edges=[[0, 1]])
        for backend in ("midpoint", "geometric"):
            first = solve(original, backend=backend, series_order=0)
            second = solve(reversed_graph, backend=backend, series_order=0)
            value = lambda result: F(result["results"][0]["result"]["certified_finite_terms"][0]["partial_sum"])
            self.assertEqual(value(first), -2)
            self.assertEqual(value(second), 2)

    def test_invalid_or_partial_tensor_requests_rejected(self):
        for source in (True, 1.0):
            document = deepcopy(TRIANGLE)
            document["sources"][0][0] = source
            with self.assertRaises(TypeError):
                solve(document, backend="midpoint")
        document = dict(ZERO, requests=[{"tensor_edges": [0, 0]}],
                        budgets={"max_tensor_entries": 3})
        with self.assertRaisesRegex(ValueError, "complete tensor"):
            solve(document, backend="midpoint")
        for request in ({"powers": [], "polynomial": []}, {"powers": [[0, 2, 1]]},
                        {"tensor_edges": [0, 0], "indices": ["a", "a"]}):
            with self.assertRaises(ValueError):
                solve(dict(ZERO, requests=[request]), backend="midpoint")
        for width in ("0", "-1/3"):
            with self.assertRaises(ValueError):
                solve(ZERO, backend="midpoint", width=width)

    def test_atomic_cli_output_and_failure_preserve_existing_result(self):
        script = Path(__file__).resolve().parents[1] / "solve_product_forest_cubes.py"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            output = path / "result.json"
            prior = {"previous_complete_result": "7/13"}
            write_complete(output, prior)
            document = dict(ZERO, requests=[{"name": "bare", "powers": []},
                                           {"name": "invalid", "powers": [[9, 0, 1]]}])
            source = path / "input.json"
            source.write_text(json.dumps(document))
            failed = subprocess.run([sys.executable, str(script), str(source), "--output", str(output),
                                     "--backend", "midpoint"], text=True, capture_output=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertEqual(json.loads(output.read_text()), prior)
            document["requests"] = [{"name": "bare", "powers": []}]
            source.write_text(json.dumps(document))
            completed = subprocess.run([sys.executable, str(script), str(source), "--output", str(output),
                                        "--backend", "midpoint", "--series-order", "0"],
                                       text=True, capture_output=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads(output.read_text())
            self.assertTrue(result["all_requests_complete"])
            self.assertEqual(result["results"][0]["result"]["certified_finite_terms"][0]["partial_sum"], "1")
            self.assertEqual({item.name for item in path.iterdir()}, {"result.json", "input.json"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
