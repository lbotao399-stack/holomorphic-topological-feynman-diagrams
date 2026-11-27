"""Independent source-complete local/core integration controls.

Run: python explore/circuit-general/sourced-forests/tests_exact.py
"""

from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

import networkx as nx
from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx, fmpq_series

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "code"))
from htk import Graph, ProductForestMoments, brute
from htk.forest_contractions import ForestContractions
from htk.forest_model import recognize_product_forest
from htk.forest_periods import CoreCurrentPeriods, ExactForestSeries
from htk.treeflow import tree_flow
from independent_duality import EDGES as HARD_EDGES, SOURCES as HARD_SOURCES, EXPECTED, replay

PAIR_EDGES = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (2, 6),
              (3, 5), (3, 6), (4, 5), (4, 6), (5, 6))
PAIR_SOURCES = ((1, 2), (3, -1), (-2, 4), (0, 0), (0, 0), (5, 3), (-7, -8))
PAIR_CORE = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (2, 4), (3, 4))
SIX_EDGES = ((0, 1), (1, 3), (0, 2), (2, 4), (0, 5), (0, 6),
             (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7))
SIX_SOURCES = ((0, 0), (1, 2), (3, -1), (0, 0), (0, 0), (0, 0), (5, 3), (-9, -4))


def bracket(left, right):
    return left[0] * right[1] - left[1] * right[0]


def word(G, powers=None):
    exponent = [0] * (2 * G.ne)
    for (edge, component), power in (powers or {}).items():
        exponent[2 * edge + component] = power
    return {tuple(exponent): F(1)}


def paired_forward(times):
    """Independent original-resistance to ratios/core map, with old AB=1."""
    def star(a, b, c):
        total = a * b + a * c + b * c
        return total / c, total / b, total / a
    left = star(times[2], times[6], times[7])
    right = star(times[3], times[8], times[9])
    ua = left[0] * right[0] / (left[0] + right[0])
    ub = left[1] * right[1] / (left[1] + right[1])
    ab = 1 / (1 / times[10] + 1 / left[2] + 1 / right[2])
    ratios = [ua / left[0], ub / left[1], ab / left[2], ab / right[2]]
    return ratios + [value / ab for value in (times[0], times[1], ua, ub, times[4], times[5])]


def determinant_by_subsets(matrix, context):
    """Polynomial determinant without division or the production determinant code."""
    states = {0: context.constant(1)}
    for row in matrix:
        following = {}
        for mask, coefficient in states.items():
            for column, value in enumerate(row):
                if value and not mask & (1 << column):
                    target = mask | (1 << column)
                    signed = coefficient * value * (-1) ** ((mask >> (column + 1)).bit_count())
                    following[target] = following.get(target, context.constant(0)) + signed
        states = following
    return states.get((1 << len(matrix)) - 1, context.constant(0))


class SourcedForestIntegration(unittest.TestCase):
    def test_paired_five_loop_nonzero_currents_and_conservation(self):
        G = Graph(7, PAIR_EDGES)
        contraction = ForestContractions(G, PAIR_SOURCES, (5, 6))
        self.assertEqual([step["kind"] for step in contraction.steps], ["paired"])
        self.assertEqual(contraction.reduced_sources,
                         tuple(tuple(map(F, PAIR_SOURCES[vertex])) for vertex in (0, 1, 2, 5, 6)))
        backend = CoreCurrentPeriods(contraction.reduced_graph, contraction.reduced_sources,
                                     contraction.reduced_marked)
        self.addCleanup(backend.close)
        self.assertEqual(backend.integrate(contraction.contract().polynomial), 0)
        up = (F(106215469, 1814400), F(164932973, 2721600))
        pA = (F(82364363, 2721600), F(294587393, 5443200))
        pB = (F(153917681, 5443200), F(35278553, 5443200))
        expected = {2: up, 3: tuple(-value for value in up),
                    6: pA, 7: pB, 8: tuple(-value for value in pA),
                    9: tuple(-value for value in pB)}
        currents = []
        for edge in range(G.ne):
            row = []
            for component in range(2):
                result = contraction.contract({(edge, component): 1})
                value = backend.integrate(result.polynomial)
                self.assertEqual(value, expected.get(edge, (F(0), F(0)))[component])
                row.append(value)
            currents.append(row)
        for vertex in range(G.nv):
            for component in range(2):
                self.assertEqual(sum((int(left == vertex) - int(right == vertex))
                                     * currents[edge][component]
                                     for edge, (left, right) in enumerate(G.edges)), 0)
        for powers, expected in (({(6, 0): 1, (10, 1): 2}, -F(2542321651, 47628000)),
                                 ({(6, 0): 2, (0, 1): 1}, F(79343557121, 6858432000))):
            self.assertEqual(backend.integrate(contraction.contract(powers).polynomial), expected)

    def test_independent_full_pair_current_and_jacobian_identity(self):
        G, reduced = Graph(7, PAIR_EDGES), Graph(5, PAIR_CORE)
        sources = tuple(PAIR_SOURCES[vertex] for vertex in (0, 1, 2, 5, 6))
        for offset in (0, 1, 3, 7):
            times = [F(1 + (i + offset) ** 2) for i in range(10)] + [F(1)]
            coordinates = paired_forward(times)
            r1, r2, r3, r4 = coordinates[:4]
            core_times = coordinates[4:] + [F(1)]
            actual = brute.currents(G, times, PAIR_SOURCES)
            core_currents = brute.currents(reduced, core_times, sources)
            left, right, marked = core_currents[2], core_currents[3], core_currents[6]
            predicted = [core_currents[0], core_currents[1],
                         tuple(r1 * left[j] + r2 * right[j] for j in range(2)),
                         tuple((1 - r1) * left[j] + (1 - r2) * right[j] for j in range(2)),
                         core_currents[4], core_currents[5],
                         tuple(r1 * left[j] - r3 * marked[j] for j in range(2)),
                         tuple(r2 * right[j] + r3 * marked[j] for j in range(2)),
                         tuple((1 - r1) * left[j] - r4 * marked[j] for j in range(2)),
                         tuple((1 - r2) * right[j] + r4 * marked[j] for j in range(2)),
                         tuple((1 - r3 - r4) * marked[j] for j in range(2))]
            self.assertEqual(tuple(map(tuple, actual)), tuple(map(tuple, predicted)))
            jacobian = []
            for variable in range(10):
                dual = [fmpq_series([fmpq(str(value)), int(edge == variable)], prec=2)
                        for edge, value in enumerate(times)]
                jacobian.append([value[1] for value in paired_forward(dual)])
            jacobian = F(str(fmpq_mat(jacobian).det()))
            self.assertGreater(jacobian, 0)
            factor = ((r1 - r2) * (1 - r3 - r4)
                      * bracket(left, marked) * bracket(right, marked))
            original_density = brute.density(G, times, PAIR_SOURCES).rational()
            self.assertNotEqual(original_density, 0)
            self.assertEqual(original_density,
                             -factor * brute.density(reduced, core_times, sources).rational() * jacobian)

    def test_bordered_current_bracket_polynomial_identity(self):
        fixtures = (
            (Graph(5, ((0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))),
             ((1, 2), (3, -1), (-2, 4), (5, 3), (-7, -8))),
            (Graph(7, HARD_EDGES),
             ((1, 2), (3, -1), (-2, 4), (1, -3), (2, 1), (5, 3), (-10, -6))),
        )
        for G, sources in fixtures:
            context = fmpq_mpoly_ctx.get(tuple("t%d" % edge for edge in range(G.ne)))
            times = context.gens()
            cycles = G.cycle_matrix()
            flow = tree_flow(G, sources)
            matrix = [[sum((cycles[i][edge] * cycles[j][edge] * times[edge]
                            for edge in range(G.ne)), context.constant(0))
                       for j in range(G.loops)] for i in range(G.loops)]
            rhs = [[sum((cycles[i][edge] * int(flow[edge][component]) * times[edge]
                         for edge in range(G.ne)), context.constant(0))
                    for component in range(2)] for i in range(G.loops)]
            U = determinant_by_subsets(matrix, context)
            for first, second in ((0, G.ne - 1), (G.ne // 2, G.ne - 2), (1, 3)):
                numerators = []
                for edge in (first, second):
                    pair = []
                    for component in range(2):
                        augmented = [row + [rhs[i][component]] for i, row in enumerate(matrix)]
                        augmented.append([context.constant(cycles[i][edge]) for i in range(G.loops)]
                                         + [context.constant(int(flow[edge][component]))])
                        pair.append(determinant_by_subsets(augmented, context))
                    numerators.append(pair)
                augmented = [row + rhs[i] for i, row in enumerate(matrix)]
                for edge in (first, second):
                    augmented.append([context.constant(cycles[i][edge]) for i in range(G.loops)]
                                     + [context.constant(int(value)) for value in flow[edge]])
                bordered = determinant_by_subsets(augmented, context)
                self.assertNotEqual(bordered, 0)
                self.assertEqual(bracket(*numerators), U * bordered)
                self.assertTrue(all(sum(exponent) == G.loops for exponent, _ in bordered.terms()))

    def test_irreducible_physical_five_loop_complete_period(self):
        G = Graph(7, HARD_EDGES)
        structure = recognize_product_forest(G, (5, 6))
        self.assertFalse(structure.cubic)
        self.assertEqual(structure.degree_two_vertices, ((1, 2),))
        self.assertEqual(structure.higher_valence_vertices, ((0,),))
        contraction = ForestContractions(G, HARD_SOURCES, (5, 6))
        self.assertEqual(contraction.steps, ())
        backend = CoreCurrentPeriods(G, HARD_SOURCES, (5, 6))
        self.addCleanup(backend.close)
        self.assertEqual(dict(backend.integrate(word(G))), EXPECTED)
        independent, details = replay()
        self.assertEqual(independent, EXPECTED)
        self.assertEqual(details["spanning_tree_monomials"], 75)
        self.assertEqual(details["minimum_convergence_margin"], 1)
        self.assertEqual(details["core_period_orbits"], 2474)

    def test_six_loop_nonplanar_physical_forest_retains_hard_core(self):
        G = Graph(8, SIX_EDGES)
        self.assertFalse(nx.check_planarity(nx.Graph(G.edges))[0])
        recognize_product_forest(G, (6, 7))
        contraction = ForestContractions(G, SIX_SOURCES, (6, 7), check=4)
        self.assertEqual([step["kind"] for step in contraction.steps], ["single"])
        self.assertEqual(contraction.steps[0]["original_centres"], (5,))
        self.assertEqual(contraction.reduced_sources, tuple(tuple(map(F, row)) for row in HARD_SOURCES))
        self.assertEqual({frozenset(edge) for edge in contraction.reduced_graph.edges},
                         {frozenset(edge) for edge in HARD_EDGES})
        result = contraction.contract()
        self.assertEqual(result.degree, 2)
        self.assertEqual(len(result.polynomial), 2)
        self.assertEqual({abs(value) for value in result.polynomial.values()}, {F(1, 4)})

    def test_actual_unavailable_core_returns_complete_series(self):
        # The non-root interaction source forbids the passive pairing.
        # Suppressing the one series vertex leaves K3,3 plus one edge at q=7.
        # This exercises the real read-only table lookup, without a mock.
        G = Graph(7, ((0, 1), (1, 2), (0, 3), (0, 4), (2, 5), (2, 6),
                      (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)))
        sources = ((0, 0), (1, 0), (0, 0), (0, 1), (0, 0), (0, 0), (-1, -1))
        with ProductForestMoments(G, sources, (5, 6)) as engine:
            self.assertEqual(engine.method, "sourced_core_periods")
            value = engine.moment()
            self.assertIsInstance(value, ExactForestSeries)
            self.assertEqual([(q, len(polynomial)) for q, polynomial in value.components],
                             [(7, 7368)])
            self.assertEqual(engine.last_stats["bare_numerator_terms"], 10577)
            self.assertEqual(engine.last_stats["finite_period_unavailable"],
                             "K33+e at q=7 not tabulated")
            self.assertTrue(engine.last_stats["all_original_sources_retained"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
