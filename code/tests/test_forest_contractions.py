"""Exact local polynomial controls, including a sourced higher-valence parent.

Run: python3 code/tests/test_forest_contractions.py
"""
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from htk import Graph
from htk.forest_contractions import (ForestContractions, ContractionBudgetExceeded,
                                    normalize_forest_sources)
from htk.forest_model import forest_graph
from htk.forest_moments import ForestMoments
from htk.linalg import perm_sign
from htk.spred import plan_generic, sp_reduce


PAIRED_GRAPH = Graph(7, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6),
                       (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)])
PAIRED_SOURCES = ((1, 2), (3, -1), (-2, 4), (0, 0), (0, 0), (5, 3), (-7, -8))
# The passive centre is the original grounded vertex. The planner must
# relabel vertices without altering edge order, orientations, or sources.
SINGLE_GRAPH = Graph(6, [(0, 1), (0, 2), (0, 3), (1, 3), (2, 4), (3, 4),
                       (5, 0), (5, 4), (5, 3)])
SINGLE_SOURCES = ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8), (0, 0))


def evaluate(polynomial, values):
    result = F(0)
    for exponents, coefficient in polynomial.items():
        term = coefficient
        for value, power in zip(values, exponents):
            term *= value ** power
        result += term
    return result


def integrate_remainder(result):
    reduction = sp_reduce(result.graph, D=2, plan=plan_generic(result.graph))
    density = reduction.density(result.sources, T=0, method="generic", check=0)

    def insertion(currents):
        value = 0
        for exponents, coefficient in result.polynomial.items():
            term = coefficient
            for position, power in enumerate(exponents):
                if power:
                    term *= currents[position // 2][position % 2] ** power
            value += term
        return value

    return reduction.integrate(result.sources, T=0, insertion=insertion,
                               method="generic", check=0, dens=density)


class ForestContractionChecks(unittest.TestCase):
    def test_paired_polynomial_and_complete_sources(self):
        engine = ForestContractions(PAIRED_GRAPH, PAIRED_SOURCES, (5, 6), check=2)
        bare = engine.contract()
        self.assertEqual(bare.polynomial, {})
        self.assertEqual(bare.graph.loops, 3)
        self.assertEqual(bare.sources, tuple(tuple(map(F, row)) for i, row in enumerate(PAIRED_SOURCES)
                                             if i not in (3, 4)))
        self.assertEqual(bare.original_to_reduced_vertices, (0, 1, 2, None, None, 3, 4))
        self.assertEqual(bare.steps[0]["original_centres"], (3, 4))
        self.assertEqual(bare.steps[0]["kind"], "paired")
        self.assertEqual(bare.marked, (3, 4))
        first = engine.contract({(6, 0): 1})
        mixed = engine.contract({(6, 0): 1, (10, 1): 2})
        self.assertIs(first.graph, bare.graph)
        self.assertEqual(first.degree, 5)
        self.assertEqual(mixed.degree, 7)
        self.assertEqual(len(first.polynomial), 4)
        self.assertTrue(all(type(power) is int for exponent in first.polynomial for power in exponent))
        for offset in (1, 4):
            values = tuple(F((i + offset) ** 2 + 3, i + offset + 1)
                           for i in range(2 * first.graph.ne))
            left, right, marked = values[:2], values[2:4], values[4:6]
            factor = ((left[0] * marked[1] - left[1] * marked[0])
                      * (right[0] * marked[1] - right[1] * marked[0]))
            self.assertNotEqual(factor, 0)
            self.assertEqual(evaluate(first.polynomial, values), -left[0] * factor / 72)
            self.assertEqual(evaluate(mixed.polynomial, values), -left[0] * marked[1] ** 2 * factor / 240)
        self.assertEqual(first.stats["global_ratio_coordinates"], 0)
        self.assertEqual(first.stats["largest_local_ratio_count"], 4)

    def test_single_with_sourced_higher_parent_and_original_root_centre(self):
        engine = ForestContractions(SINGLE_GRAPH, SINGLE_SOURCES, (3, 4), check=2)
        result = engine.contract()
        self.assertEqual(result.graph.loops, 3)
        self.assertEqual(result.steps[0]["kind"], "single")
        self.assertEqual(result.steps[0]["original_centres"], (5,))
        self.assertEqual(result.original_to_reduced_vertices, (0, 1, 2, 3, 4, None))
        self.assertEqual(result.sources, tuple(tuple(map(F, row)) for row in SINGLE_SOURCES[:-1]))
        values = tuple(F(i * i + 1, i + 2) for i in range(2 * result.graph.ne))
        bracket = values[2] * values[5] - values[3] * values[4]
        self.assertNotEqual(bracket, 0)
        self.assertEqual(evaluate(result.polynomial, values), result.steps[0]["kernel_sign"] * bracket / 4)

    def test_complete_cubic_controls_with_edge_and_vertex_changes(self):
        rng = random.Random(913)
        cases = (("A", ("A", "B")), (("A", "B"), ("A", "B")),
                 ("B", (("A", "B"), ("A", "B"))))
        for tree in cases:
            G, roots = forest_graph((tree,))
            sources = [[0, 0] for _ in range(G.nv)]
            sources[roots[0]], sources[-2], sources[-1] = [1, 2], [5, 3], [-6, -5]
            marked = (G.nv - 2, G.nv - 1)
            inserted_edge = 2 if G.loops == 4 else 0
            word = {(inserted_edge, 0): 1, (G.ne - 1, 1): 2}
            reference = ForestMoments(G, sources, marked).moment(word)
            self.assertNotEqual(reference, 0)
            for changed in (False, True):
                if changed:
                    order, vertices = list(range(G.ne)), list(range(G.nv))
                    rng.shuffle(order)
                    rng.shuffle(vertices)
                    reversals = set(range(0, G.ne, 2))
                    changed_G = G.reverse_edges(reversals).permute_edges(order).relabel(vertices)
                    changed_sources = [None] * G.nv
                    for vertex, row in enumerate(sources):
                        changed_sources[vertices[vertex]] = row
                    changed_word = {(order.index(edge), component): power
                                    for (edge, component), power in word.items()}
                    sign = perm_sign(order) * (-1) ** sum(power for (edge, _), power in word.items()
                                                         if edge in reversals)
                    changed_marked = tuple(vertices[v] for v in marked)
                else:
                    changed_G, changed_sources, changed_word, changed_marked = G, sources, word, marked
                    sign = 1
                result = ForestContractions(changed_G, changed_sources, changed_marked, check=1).contract(changed_word)
                self.assertEqual(integrate_remainder(result), sign * reference)

    def test_sourced_stars_remain_and_sources_formats_agree(self):
        G, _ = forest_graph(((("A", "B"), ("A", "B")),))
        full = ((1, 2), (3, -1), (0, 0), (5, 3), (-9, -4))
        for sources in (full, full[:-1], dict(enumerate(full))):
            engine = ForestContractions(G, sources, (3, 4), check=2)
            result = engine.contract({(0, 0): 2, (6, 1): 1})
            self.assertFalse(engine.steps)
            self.assertEqual(result.sources, tuple(tuple(map(F, row)) for row in full))
            self.assertEqual(result.graph.edges, G.edges)
            self.assertEqual(len(result.polynomial), 1)
            self.assertEqual(result.degree, 3)

    def test_nineteen_loop_higher_valence_comb_exact_boundary_formula(self):
        # Invert the paired move repeatedly on one passive child. The
        # original four-valent parent and all five nonzero source rows stay.
        for pairs in (1, 3, 8):
            edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (2, 4), (3, 4)]
            sources = [[1, 2], [3, -1], [-2, 4], [5, 3], [-7, -8]]
            parent = 0
            for _ in range(pairs):
                edges.remove((parent, 3))
                edges.remove((parent, 4))
                first, second = len(sources), len(sources) + 1
                sources += [[0, 0], [0, 0]]
                edges += [(first, parent), (first, 3), (first, 4),
                          (second, parent), (second, 3), (second, 4)]
                parent = first
            G = Graph(len(sources), edges)
            word = {(edges.index((parent, 3)), 0): 1}
            result = ForestContractions(G, sources, (3, 4), check=0).contract(word)
            coefficient = F(2 * (-1) ** pairs, factorial(pairs) * factorial(pairs + 1) ** 2
                            * factorial(pairs + 2) * factorial(2 * pairs + 1))
            expected = {}
            for first in range(pairs + 1):
                for second in range(pairs + 1):
                    powers = (pairs - first + 1, first, pairs - second, second,
                              first + second, 2 * pairs - first - second) + (0,) * 8
                    expected[powers] = (coefficient * (-1) ** (first + second)
                                        * comb(pairs, first) * comb(pairs, second))
            self.assertEqual(result.polynomial, expected)
            self.assertEqual(result.graph.loops, 3)
            self.assertEqual(result.sources, tuple(tuple(map(F, row)) for row in sources[:5]))
            self.assertEqual(result.stats["step_term_counts"], [j * j for j in range(2, pairs + 2)])
            self.assertLessEqual(result.stats["largest_polynomial_terms"], 15360)

    def test_polynomial_linearity_validation_and_budget(self):
        engine = ForestContractions(PAIRED_GRAPH, PAIRED_SOURCES, (5, 6), check=0)
        left = [0] * (2 * PAIRED_GRAPH.ne)
        right = list(left)
        left[12], right[15] = 1, 1
        polynomial = {tuple(left): F(2, 3), tuple(right): F(-7, 5)}
        result = engine.contract(polynomial=polynomial)
        terms = {}
        for factor, word in ((F(2, 3), {(6, 0): 1}), (F(-7, 5), {(7, 1): 1})):
            for exponents, value in engine.contract(word).polynomial.items():
                terms[exponents] = terms.get(exponents, F(0)) + factor * value
        self.assertEqual(result.polynomial, {key: value for key, value in terms.items() if value})
        for powers in ({(-1, 0): 1}, {(0, 2): 1}, {(0, 0): -1}, {(0, 0): 1.5}):
            with self.assertRaises((TypeError, ValueError)):
                engine.contract(powers)
        with self.assertRaises(ValueError):
            engine.contract(polynomial={tuple(left): 1, (0,) * len(left): 2})
        with self.assertRaises(TypeError):
            engine.contract(polynomial={tuple(left): 0.1})
        with self.assertRaises(ValueError):
            engine.contract({}, polynomial={tuple(left): 1})
        with self.assertRaises(ContractionBudgetExceeded):
            ForestContractions(PAIRED_GRAPH, PAIRED_SOURCES, (5, 6), check=0,
                               max_terms=1).contract({(6, 0): 1})
        for sources in ([[0.0, 0]] * 7, [[0]] * 7, [[0, 0]] * 6 + [[1, 0]]):
            with self.assertRaises((TypeError, ValueError)):
                normalize_forest_sources(PAIRED_GRAPH, sources)


if __name__ == "__main__":
    unittest.main()
