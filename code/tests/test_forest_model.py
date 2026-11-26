"""Exact marked-forest recognition, canonical graph maps, and source checks.

Run: python3 code/tests/test_forest_model.py
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import networkx as nx
from htk import brute, catalog
from htk.forest_model import (
    CubicForest, NonCubicProductForest, NotProductForest, UnsupportedForestSources,
    forest_counts, forest_graph, forests, recognize_product_forest, trees,
)
from htk.graph import Graph


class ForestModelChecks(unittest.TestCase):
    def test_complete_small_enumeration_and_count_recurrence(self):
        tree_counts, counts = forest_counts(9)
        self.assertEqual(tree_counts, (0, 1, 2, 5, 12, 32, 86, 243, 696, 2045))
        self.assertEqual(counts, (1, 1, 3, 8, 23, 65, 194, 576, 1753, 5369))
        for n in range(7):
            objects = list(forests(n))
            self.assertEqual(len(objects), counts[n])
            self.assertEqual(len(set(objects)), counts[n])
            self.assertEqual(len(trees(n)), tree_counts[n])
            networks = []
            for forest in objects:
                G, roots = forest_graph(forest)
                model = CubicForest(G, marked=(n, n + 1))
                self.assertEqual(model.forest, forest)
                self.assertEqual(model.canonical_roots, roots)
                self.assertEqual(model.canonical_graph.edges, G.edges)
                self.assertEqual(model.kernel_sign, 1)
                self.assertTrue(catalog.tight_ok(G.to_networkx(), 2))
                if n <= 4:
                    graph = nx.Graph(G.to_networkx())
                    nx.set_node_attributes(graph, {v: "A" if v == n else "B" if v == n + 1 else "internal"
                                                   for v in graph}, "mark")
                    networks.append(graph)
            for first, second in combinations(networks, 2):
                self.assertFalse(nx.is_isomorphic(first, second,
                    node_match=lambda a, b: a["mark"] == b["mark"]))
        # This count uses only the recurrence, not an exponential census.
        self.assertGreater(forest_counts(100)[1][-1], 10**40)

    def test_broader_non_cubic_forest_is_distinct(self):
        # One degree-four interaction, two degree-two surviving fields.
        G = Graph(7, [(0, 1), (0, 2), (0, 3), (0, 4),
                      (1, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
        self.assertTrue(catalog.tight_ok(G.to_networkx(), 2))
        structure = recognize_product_forest(G, (5, 6))
        self.assertFalse(structure.cubic)
        self.assertEqual(structure.degree_two_vertices, ((1, 2),))
        self.assertEqual(structure.higher_valence_vertices, ((0,),))
        with self.assertRaises(NonCubicProductForest) as failure:
            CubicForest(G, marked=(5, 6))
        self.assertEqual(failure.exception.structure, structure)

    def test_invalid_topologies_and_exact_subgraph_bound(self):
        cyclic = Graph(5, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 3), (3, 4)])
        # The four-vertex set {0,1,A,B} has six edges, exceeding its bound five.
        dense_subtree = Graph(6, [(0, 1), (1, 2), (2, 3), (0, 4), (0, 5),
                                  (1, 4), (1, 5), (3, 4), (4, 5)])
        self.assertFalse(catalog.tight_ok(dense_subtree.to_networkx(), 2))
        for G, marked in ((cyclic, (3, 4)), (dense_subtree, (4, 5)),
                          (Graph(3, [(0, 1), (1, 2)]), (1, 2)),
                          (Graph(3, [(0, 1), (1, 2), (0, 1)]), (1, 2))):
            with self.assertRaises(NotProductForest):
                CubicForest(G, marked=marked)
        G, _ = forest_graph((("A", "B"),))
        for marked in ((1, 1), (1, 3), (True, 2), (1,), (1, 2, 3)):
            with self.assertRaises((TypeError, ValueError)):
                CubicForest(G, marked=marked)
        for invalid in (True, 1.0, -1):
            with self.assertRaises((TypeError, ValueError)):
                trees(invalid)

    def test_linear_subtree_criterion_against_every_induced_subset(self):
        checked = 0
        for n in range(1, 5):
            interior_trees = [nx.empty_graph(1)] if n == 1 else [
                nx.from_prufer_sequence(sequence) for sequence in product(range(n), repeat=n - 2)]
            for interior in interior_trees:
                for masks in product(range(4), repeat=n):
                    if sum(mask.bit_count() for mask in masks) != n + 1:
                        continue
                    edges = list(interior.edges()) + [(n, n + 1)]
                    edges += [(v, n + boundary) for v, mask in enumerate(masks)
                              for boundary in range(2) if mask & (1 << boundary)]
                    G = Graph(n + 2, edges)
                    if any(G.degree(v) < 2 for v in range(n)) or sum(G.degree(v) > 3 for v in range(n)) > 1:
                        continue
                    exact = catalog.tight_ok(G.to_networkx(), 2)
                    try:
                        recognize_product_forest(G, (n, n + 1))
                        recognized = True
                    except NotProductForest:
                        recognized = False
                    self.assertEqual(recognized, exact)
                    checked += 1
        self.assertGreater(checked, 200)

    def test_arbitrary_graph_orders_against_exact_original_density(self):
        rng = random.Random(3621)
        checks = 0
        for n in range(1, 5):
            for forest in forests(n):
                G, roots = forest_graph(forest)
                sources = [[F(0), F(0)] for _ in range(G.nv)]
                sources[n] = [F(3), F(-2)]
                for i, vertex in enumerate(roots):
                    sources[vertex] = [F(i + 1, i + 2), F(2 * i - 3, i + 1)]
                sources[n + 1] = [-sum(row[i] for row in sources[:-1]) for i in range(2)]
                labels = list(range(G.nv))
                order = list(range(G.ne))
                rng.shuffle(labels)
                rng.shuffle(order)
                changed = G.relabel(labels).permute_edges(order).reverse_edges(
                    [edge for edge in range(G.ne) if rng.randrange(2)])
                moved = [None] * G.nv
                for old, new in enumerate(labels):
                    moved[new] = sources[old]
                model = CubicForest(changed, marked=(labels[n], labels[n + 1]))
                canonical_sources = model.normalize_sources(moved)
                times = [F(edge + 2, edge + 1) for edge in range(G.ne)]
                canonical_times = [times[edge] for edge in model.canonical_to_original_edges]
                original = brute.density(changed, times, moved).rational()
                canonical = brute.density(model.canonical_graph, canonical_times, canonical_sources).rational()
                self.assertEqual(original, model.kernel_sign * canonical)
                currents = brute.currents(changed, times, moved)
                canonical_currents = brute.currents(model.canonical_graph, canonical_times, canonical_sources)
                for edge in range(G.ne):
                    for component in range(2):
                        self.assertEqual(currents[edge][component], model.original_edge_signs[edge]
                                         * canonical_currents[model.original_to_canonical_edges[edge]][component])
                powers = {(0, 0): 3, (G.ne - 1, 1): 2}
                mapped, sign = model.map_powers(powers)
                self.assertEqual(original * brute.monomial(currents, powers),
                                 sign * canonical * brute.monomial(canonical_currents, mapped))
                checks += 1
        self.assertEqual(checks, 35)

    def test_source_conservation_formats_and_unsupported_sources(self):
        G, roots = forest_graph((("A", "B"), ("A", ("A", "B"))))
        model = CubicForest(G, marked=(3, 4))
        full = [[1, 2], [-3, 4], [0, 0], [5, -1], [-3, -5]]
        expected = tuple(tuple(map(F, row)) for row in full)
        for value in (full, full[:-1], dict(enumerate(full)), dict(enumerate(full[:-1]))):
            self.assertEqual(model.normalize_sources(value), expected)
        self.assertEqual(model.normalize_sources(None), ((F(0), F(0)),) * G.nv)
        unsupported = [list(row) for row in full]
        unsupported[2][0] = 1
        unsupported[4][0] -= 1
        with self.assertRaises(UnsupportedForestSources):
            model.normalize_sources(unsupported)
        with self.assertRaises(ValueError):
            model.normalize_sources(full[:-1] + [[0, 0]])
        with self.assertRaises(ValueError):
            model.normalize_sources([[1]] * G.nv)
        with self.assertRaises(TypeError):
            model.normalize_sources([[0.0, 0.0]] * G.nv)
        with self.assertRaises(ValueError):
            model.normalize_sources({-1: (0, 0)})
        for powers in ({(-1, 0): 1}, {(0, 2): 1}, {(0, 0): -1}, {(0, 0): 1.5}):
            with self.assertRaises((TypeError, ValueError)):
                model.map_powers(powers)


if __name__ == "__main__":
    unittest.main()
