"""Independent exact controls for complete product-forest leaf elimination."""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk import brute
from htk.forest_kron import ForestKronPlan, KronBudgetExceeded
from htk.forest_model import NotProductForest, forest_graph
from htk.graph import Graph, fan
from htk.linalg import det, perm_sign
from htk.treeflow import tree_flow


FORK = ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
HIGHER = ((0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
          (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7))


def all_sources(graph):
    rows = [(F(2 * vertex + 1, vertex + 1), F((-1) ** vertex * (vertex + 3), vertex + 2))
            for vertex in range(graph.nv - 1)]
    rows.append(tuple(-sum(row[a] for row in rows) for a in range(2)))
    return tuple(rows)


def times(graph):
    return tuple(F(edge + 2, edge + 1) for edge in range(graph.ne))


def tree_average_derivatives(graph, resistances):
    """Differentiate the original tree sum, independently of any Schur solve."""
    size, edges = graph.nv - 1, graph.ne
    units = [[int(i == j) for j in range(size)] for i in range(size)]
    numerator = [[F(0)] * size for unused in range(edges)]
    derivative = [[[F(0)] * edges for unused in range(size)] for edge in range(edges)]
    U, dU = F(0), [F(0)] * edges
    for tree in combinations(range(edges), size):
        try:
            flow = tree_flow(graph, units, tree)
        except ValueError:
            continue
        complement = tuple(edge for edge in range(edges) if edge not in tree)
        weight = prod(resistances[edge] for edge in complement)
        U += weight
        for edge in complement:
            dU[edge] += weight / resistances[edge]
        for edge in range(edges):
            for source in range(size):
                numerator[edge][source] += weight * flow[edge][source]
                for f in complement:
                    derivative[edge][source][f] += weight * flow[edge][source] / resistances[f]
    response = tuple(tuple(value / U for value in row) for row in numerator)
    derivatives = tuple(tuple(tuple((derivative[edge][source][f] * U - numerator[edge][source] * dU[f]) / U ** 2
                                   for f in range(edges)) for source in range(size)) for edge in range(edges))
    return U, response, derivatives


class ForestKronTests(unittest.TestCase):
    def test_complete_sources_on_branched_multicomponent_and_higher_valence_graphs(self):
        multiple, unused = forest_graph((("A", "B"), (("A", "B"), ("A", "B"))))
        graphs = ((fan(2), (2, 3)), (Graph(5, FORK), (3, 4)),
                  (multiple, (multiple.nv - 2, multiple.nv - 1)), (Graph(8, HIGHER), (6, 7)))
        for graph, marked in graphs:
            with self.subTest(graph=graph):
                plan = ForestKronPlan(graph, marked)
                result = plan.response(times(graph))
                sources = all_sources(graph)
                self.assertEqual(len(plan.steps), graph.nv - 2)
                self.assertEqual(result.current_response, tuple(map(tuple, brute.response(graph, result.resistances))))
                self.assertEqual(result.currents(sources), tuple(map(tuple, brute.currents(graph, result.resistances, sources))))
                self.assertEqual(result.potentials(sources), tuple(map(tuple, brute.potentials(graph, result.resistances, sources))))
                self.assertEqual(result.laplacian_determinant, det(brute.laplacian(graph, result.resistances)))
                self.assertEqual(result.U, brute.symanzik_U(graph, result.resistances))
                self.assertEqual(result.U, prod(result.resistances) * prod(result.pivots) * result.retained_conductance)
                self.assertEqual(result.density(sources), brute.density(graph, result.resistances, sources).rational())
                self.assertEqual(result.stats["expanded_Kirchhoff_terms"], 0)
        self.assertGreater(ForestKronPlan(Graph(5, FORK), (3, 4)).stats["non_Laman_intermediate_stages"], 0)

    def test_every_response_derivative_against_original_spanning_tree_flow(self):
        graph = Graph(5, FORK)
        result = ForestKronPlan(graph, (3, 4)).response(times(graph))
        U, response, derivatives = tree_average_derivatives(graph, result.resistances)
        self.assertEqual(result.U, U)
        self.assertEqual(result.current_response, response)
        self.assertEqual(result.source_jacobian(), derivatives)
        sources = all_sources(graph)
        jacobian = result.current_jacobian(sources)
        for edge in range(graph.ne):
            for a in range(2):
                for f in range(graph.ne):
                    expected = sum(derivatives[edge][v][f] * sources[v][a] for v in range(graph.nv - 1))
                    self.assertEqual(jacobian[edge][a][f], expected)
                self.assertEqual(sum(result.resistances[f] * jacobian[edge][a][f] for f in range(graph.ne)), 0)

    def test_all_source_transfers_conserve_full_rows_and_all_root_choices_agree(self):
        graph = Graph(8, HIGHER)
        resistances, original = times(graph), all_sources(graph)
        reference = None
        for root in range(6):
            result = ForestKronPlan(graph, (6, 7), roots=(root,)).response(resistances)
            source_rows = {vertex: row for vertex, row in enumerate(original)}
            removed = result.transferred_sources(original)
            for transfer, expected in zip(result.transfers, removed["removed_sources"]):
                source = source_rows.pop(transfer.vertex)
                self.assertEqual(source, expected)
                self.assertEqual(sum(transfer.weights), 1)
                for neighbor, weight in zip(transfer.neighbors, transfer.weights):
                    source_rows[neighbor] = tuple(source_rows[neighbor][a] + weight * source[a] for a in range(2))
                self.assertEqual(tuple(sum(row[a] for row in source_rows.values()) for a in range(2)), (0, 0))
            self.assertEqual(tuple(source_rows[vertex] for vertex in (6, 7)), removed["retained_sources"])
            value = (result.current_response, result.U, result.density(original))
            if reference is None:
                reference = value
            self.assertEqual(value, reference)

    def test_original_labels_orientations_gauges_and_scale(self):
        graph, marked = Graph(5, FORK), (3, 4)
        resistances, sources = times(graph), all_sources(graph)
        original = ForestKronPlan(graph, marked).response(resistances)
        density = original.density(sources)
        for gauge in range(graph.ne):
            self.assertEqual(original.density(sources, gauge_edge=gauge), density)
        scaled = ForestKronPlan(graph, marked).response(tuple(3 * value for value in resistances))
        self.assertEqual(scaled.current_response, original.current_response)
        self.assertEqual(scaled.U, 3 ** graph.loops * original.U)
        self.assertEqual(scaled.density(sources), density / 3 ** graph.ne)
        rng = random.Random(3001)
        order, labels = list(range(graph.ne)), list(range(graph.nv))
        rng.shuffle(order)
        rng.shuffle(labels)
        reversed_edges = {0, 2, 5}
        changed = graph.reverse_edges(reversed_edges).permute_edges(order).relabel(labels)
        new_sources = [None] * graph.nv
        for vertex, source in enumerate(sources):
            new_sources[labels[vertex]] = source
        result = ForestKronPlan(changed, tuple(labels[v] for v in marked)).response(tuple(resistances[edge] for edge in order))
        expected = tuple(tuple((-1 if edge in reversed_edges else 1) * value for value in original.currents(sources)[edge])
                         for edge in order)
        self.assertEqual(result.currents(new_sources), expected)
        self.assertEqual(result.U, original.U)
        self.assertEqual(result.density(new_sources), perm_sign(order) * density)
        self.assertEqual(result.current_response, tuple(map(tuple, brute.response(changed, result.resistances))))

    def test_high_valence_and_zero_loop_without_global_polynomial(self):
        # A superpotential centre of degree13, with eleven fully sourced
        # degree-two children, plus the two ordered marked vertices.
        leaves, A, B = 11, 12, 13
        edges = [(0, A), (0, B), (A, B)]
        for vertex in range(1, leaves + 1):
            edges.extend(((0, vertex), (vertex, A if vertex % 2 else B)))
        graph = Graph(14, edges)
        plan = ForestKronPlan(graph, (A, B))
        result = plan.response(times(graph))
        self.assertEqual(graph.degree(0), 13)
        self.assertEqual(len(result.transfers), 12)
        self.assertEqual(result.current_response, tuple(map(tuple, brute.response(graph, result.resistances))))
        self.assertEqual(result.U, brute.symanzik_U(graph, result.resistances))
        empty = ForestKronPlan(Graph(2, [(1, 0)]), (0, 1)).response((F(5, 3),))
        self.assertEqual(empty.U, 1)
        self.assertEqual(empty.currents(((2, 7), (-2, -7))), ((-2, -7),))
        self.assertEqual(empty.density(((2, 7), (-2, -7))), F(3, 5))

    def test_source_formats_exact_inputs_and_explicit_budgets(self):
        graph = Graph(8, HIGHER)
        plan = ForestKronPlan(graph, (6, 7))
        result, sources = plan.response(times(graph)), all_sources(graph)
        expected = result.currents(sources)
        self.assertEqual(result.currents(sources[:-1]), expected)
        self.assertEqual(result.currents(dict(enumerate(sources))), expected)
        self.assertEqual(result.currents(dict(enumerate(sources[:-1]))), expected)
        self.assertEqual(result.currents(None), ((0, 0),) * graph.ne)
        with self.assertRaises(KronBudgetExceeded):
            ForestKronPlan(graph, (6, 7), max_entries=10)
        bounded = ForestKronPlan(graph, (6, 7), max_entries=plan.response_entry_bound).response(times(graph))
        with self.assertRaises(KronBudgetExceeded):
            bounded.source_jacobian()
        self.assertEqual(bounded.currents(sources), expected)
        for value in (0, -1, F(0)):
            with self.assertRaises(ValueError):
                plan.response((value,) + times(graph)[1:])
        for value in (1.0, True, 1j):
            with self.assertRaises(TypeError):
                plan.response((value,) + times(graph)[1:])
            changed = [list(row) for row in sources]
            changed[0][0] = value
            with self.assertRaises(TypeError):
                result.currents(changed)
        with self.assertRaises(ValueError):
            plan.response(times(graph)[:-1])
        with self.assertRaises(ValueError):
            result.currents(sources[:-2])
        with self.assertRaises(ValueError):
            result.currents([(1, 1)] * graph.nv)
        with self.assertRaises(ValueError):
            ForestKronPlan(graph, (6, 7), roots=(7,))
        with self.assertRaises(ValueError):
            result.density(sources, gauge_edge=graph.ne)
        with self.assertRaises(NotProductForest):
            ForestKronPlan(Graph(5, [(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 4), (3, 4)]), (3, 4))
        with self.assertRaises(NotProductForest):
            ForestKronPlan(Graph(3, [(0, 1), (1, 0), (1, 2)]), (1, 2))


if __name__ == "__main__":
    unittest.main()
