"""Independent lookup, graph invariance, convergence, and missing-data checks.

Run: python3 explore/circuit-general/period-dispatch/test_k5_bridge.py
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest

import k5_bridge as bridge


def subset_cycle_ranks():
    ranks = []
    for mask in range(1, (1 << 10) - 1):
        parent = list(range(5))
        def root(vertex):
            while parent[vertex] != vertex:
                vertex = parent[vertex]
            return vertex
        rank = 0
        for edge, (u, v) in enumerate(bridge.E5):
            if mask & (1 << edge):
                u, v = root(u), root(v)
                if u == v:
                    rank += 1
                else:
                    parent[u] = v
        ranks.append((tuple(e for e in range(10) if mask & (1 << e)), rank))
    return ranks


class K5BridgeChecks(unittest.TestCase):
    def test_independent_exact_values_across_levels(self):
        # References: existing k5/check_ref.py and k5/REPORT.md sections 4--5.
        examples = [
            (3, (1, 1, 2, 2, 2, 2, 2, 2, 2, 2), {"z3": F(5, 3), "z5": F(-5, 3)}),
            (3, (1, 2, 2, 2, 2, 2, 2, 1, 2, 2), {"z3": F(-2), "z5": F(5, 2)}),
            (4, (1, 1, 2, 2, 3, 3, 3, 3, 3, 3),
             {"1": F(-5, 12), "z3": F(-13, 12), "z5": F(5, 3)}),
            (5, (3,) * 10, {"1": F(-13, 15120), "z3": F(19, 720), "z5": F(-5, 168)}),
            (6, (1, 1, 3, 5, 5, 5, 3, 5, 5, 3),
             {"1": F(-719897, 24494400), "z3": F(-1003, 12600), "z5": F(38, 315)}),
            (7, (2, 2, 3, 4, 6, 3, 6, 6, 4, 6),
             {"1": F(-229357639, 2309472000), "z3": F(-2840093, 10692000), "z5": F(109, 270)}),
            (8, (5, 5, 6, 1, 6, 4, 5, 6, 5, 5),
             {"1": F(105943, 1569209241600), "z3": F(-17, 302702400)}),
        ]
        for q, powers, expected in examples:
            result = bridge.k5_period_with_provenance(bridge.E5, powers, q)
            self.assertEqual(result["value"], expected)
            self.assertEqual(len(result["table"]["sha256"]), 64)
            result["value"].clear()
            self.assertEqual(bridge.k5_period(bridge.E5, powers, q), expected)

    def test_all_vertex_permutations_and_edge_orientations(self):
        powers = (1, 1, 2, 2, 3, 3, 3, 3, 3, 3)
        expected = {"1": F(-5, 12), "z3": F(-13, 12), "z5": F(5, 3)}
        rng = random.Random(906)
        for permutation in permutations(range(5)):
            order = list(range(10))
            rng.shuffle(order)
            edges, reordered = [], []
            for e in order:
                u, v = bridge.E5[e]
                edge = (11 + 7 * permutation[u], 11 + 7 * permutation[v])
                edges.append(edge if rng.randrange(2) else edge[::-1])
                reordered.append(powers[e])
            self.assertEqual(bridge.k5_period(edges, reordered, 4), expected)

    def test_short_convergence_criterion_against_every_subgraph(self):
        ranks = subset_cycle_ranks()
        counts = {}
        for q in (3, 4):
            checked = 0
            canonical = set()
            for powers in product(range(1, q), repeat=10):
                if sum(powers) != 6 * q:
                    continue
                exact = all(sum(powers[e] for e in selected) > q * rank
                            for selected, rank in ranks)
                self.assertEqual(bridge.convergence_failure(powers, q) is None, exact)
                if exact:
                    self.assertTrue(bridge.k5_period(bridge.E5, powers, q))
                    canonical.add(bridge.canonical_powers(powers))
                checked += 1
            counts[q] = (checked, len(canonical))
        self.assertEqual(counts[3][1], 2)
        self.assertEqual(counts[4][1], 30)
        print("exhaustive K5 convergence checks:", counts, flush=True)

    def test_signature_key_against_full_permutation_oracle(self):
        correspondence = {}
        inverse = {}
        for q in (3, 4):
            for powers in product(range(1, q), repeat=10):
                if sum(powers) != 6 * q or bridge.convergence_failure(powers, q) is not None:
                    continue
                key = bridge.signature_key(powers)
                full = bridge.canonical_powers(powers)
                self.assertEqual(correspondence.setdefault((q, key), full), full)
                self.assertEqual(inverse.setdefault((q, full), key), key)
        rows = json.loads((bridge.TABLE_DIRECTORY / "K5_q8.json").read_text())["values"]
        rng = random.Random(277)
        for powers, _, _ in rng.sample(rows, 300):
            powers = tuple(powers)
            key = bridge.signature_key(powers)
            full = bridge.canonical_powers(powers)
            for _ in range(4):
                permutation = rng.choice(bridge._PERMUTATIONS)
                changed = tuple(powers[e] for e in permutation)
                self.assertEqual(bridge.signature_key(changed), key)
                self.assertEqual(bridge.canonical_powers(changed), full)

    def test_public_period_dispatch_without_table_builds(self):
        from unittest.mock import patch
        folder = Path(__file__).resolve().parents[1] / "fiveloop"
        sys.path.insert(0, str(folder))
        import per5
        expected = {"1": F(-13, 15120), "z3": F(19, 720), "z5": F(-5, 168)}
        with patch.object(per5.msys.MultiSystem, "build", side_effect=AssertionError("unexpected table build")):
            self.assertEqual(per5.period(bridge.E5, [3] * 10, 5), expected)
            self.assertEqual(per5.lower(list(bridge.E5), [3] * 10, 5), expected)
            constant, unknown = per5.reducer(5).classify(list(bridge.E5), [3] * 10, 5)
            self.assertEqual(constant, expected)
            self.assertFalse(unknown)
            with self.assertRaises(per5.UnknownCore):
                per5.period(bridge.E5, [5] * 6 + [6] * 4, 9)

    def test_bulk_orbit_contraction(self):
        lookup = bridge.K5Lookup(bridge.E5, 4)
        powers = (1, 1, 2, 2, 3, 3, 3, 3, 3, 3)
        changed = tuple(powers[e] for e in bridge._PERMUTATIONS[41])
        self.assertEqual(lookup.key(powers), lookup.key(changed))
        grouped = {}
        for coefficient, vector in ((17, powers), (-9, changed)):
            key = lookup.key(vector)
            grouped[key] = grouped.get(key, 0) + coefficient
        result = {}
        for key, coefficient in grouped.items():
            for label, value in lookup.values[key].items():
                result[label] = result.get(label, F(0)) + coefficient * value
        self.assertEqual(result, {label: 8 * value for label, value in lookup.period(powers).items()})
        with self.assertRaises(TypeError):
            lookup.values[lookup.key(powers)]["1"] = F(7)
        with self.assertRaises(bridge.K5NotConvergent):
            lookup.key([1] * 10)

    def test_missing_and_divergent_are_not_zero(self):
        with self.assertRaises(bridge.K5NotConvergent):
            bridge.k5_period(bridge.E5, (1,) * 10, 3)
        with self.assertRaises(bridge.K5NotConvergent):
            bridge.k5_period(bridge.E5, (0, 2) + (2,) * 8, 3)
        with self.assertRaises(bridge.K5Unavailable):
            bridge.k5_period(bridge.E5, (5,) * 6 + (6,) * 4, 9)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(bridge.K5Unavailable):
                bridge.k5_period(bridge.E5, (1, 1) + (2,) * 8, 3, table_directory=directory)
            path = Path(directory) / "K5_q3.json"
            path.write_text(json.dumps(dict(core="K5", edges=bridge.E5, q=3, values=[])))
            with self.assertRaises(bridge.K5Unavailable):
                bridge.k5_period(bridge.E5, (1, 1) + (2,) * 8, 3, table_directory=directory)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "K5_q3.json"
            row = [[1, 1] + [2] * 8, "T1", {}]
            path.write_text(json.dumps(dict(core="K5", edges=bridge.E5, q=3, values=[row])))
            with self.assertRaises(bridge.K5TableError):
                bridge.k5_period(bridge.E5, row[0], 3, table_directory=directory)
        self.assertFalse(bridge.is_k5_graph(bridge.E5[:-1] + (bridge.E5[0],)))
        with self.assertRaises(ValueError):
            bridge.k5_period(bridge.E5, (1.0, 1) + (2,) * 8, 3)

    def test_tracked_inventory(self):
        records = bridge.inventory()
        self.assertEqual({row["q"]: row["orbits"] for row in records},
                         {3: 2, 4: 30, 5: 312, 6: 2107, 7: 10512, 8: 41428})
        self.assertEqual(sum(row["orbits"] for row in records), 54391)


if __name__ == "__main__":
    unittest.main()
