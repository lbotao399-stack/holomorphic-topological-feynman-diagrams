"""Independent exact controls for the local cubic-forest moment recursion.

Run from any directory:
    python explore/circuit-general/forest-moments/tests_exact.py

The reference evaluator expands the original global ratio chart.  The new
engine contracts only local three-current messages.  No reference values
are obtained from the new engine itself.
"""

from fractions import Fraction
from math import factorial
from pathlib import Path
import sys
import unittest

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "code"))
sys.path.insert(0, str(REPO / "explore/nonSP/verify/full_ward"))

from htk import Graph
from htk.forest_moments import ForestBudgetExceeded, ForestMoments
from forest_jets import chart, integrate_factors
from forest_orders import forests, graph


def source_assignment(n, roots, insertion=(2, -3)):
    """Every root direction is distinct, and none is the A-source direction."""
    source = [(0, 0)] * (n + 2)
    source[n] = tuple(insertion)
    for index, root in enumerate(roots, 1):
        source[root] = (index, index * index + 1)
    source[n + 1] = tuple(-sum(row[j] for row in source) for j in range(2))
    return source


def old_chart_moment(form, edges, powers):
    factors, currents, simplices = form
    inserted = list(factors)
    for (edge_index, component), power in powers.items():
        if power:
            inserted.append(currents[edges[edge_index]][component] ** power)
    return Fraction(str(integrate_factors(inserted, simplices)))


def old_reference(edges, A, B, sources, powers):
    return old_chart_moment(chart(edges, A, B, sources, factored=True), edges, powers)


def fan_tree(n, terminal="A"):
    result = "(AB)"
    for _ in range(n - 1):
        result = "(" + terminal + result + ")"
    return result


def fork_family(k):
    A, B = 3 * k, 3 * k + 1
    edges, roots, internal, physical = [], [], [], []
    for block in range(k):
        root, left, right = 3 * block, 3 * block + 1, 3 * block + 2
        roots.append(root)
        internal.append((len(edges), 0))
        physical.append((len(edges) + 1, 0))
        edges.extend(((left, root), (left, A), (left, B),
                      (right, root), (right, A), (right, B)))
    edges.append((A, B))
    sources = [(0, 0)] * (3 * k) + [(1, 0), (-1, -k)]
    for root in roots:
        sources[root] = (0, 1)
    return Graph(3 * k + 2, edges), sources, (A, B), internal, physical


class ExactForestMoments(unittest.TestCase):
    def test_complete_one_and_two_loop_bare_controls(self):
        def bracket(left, right):
            return left[0] * right[1] - left[1] * right[0]

        for root, second, insertion in (
                ((1, 2), (2, 5), (2, -3)),
                ((Fraction(2, 7), Fraction(-3, 5)), (Fraction(1, 2), Fraction(3, 11)),
                 (Fraction(7, 3), Fraction(1, 11)))):
            for forest in (("(AB)",), ("(AB)", "(AB)"),
                           ("(A(AB))",), ("(B(AB))",)):
                edges, A, B, roots, _ = graph(forest)
                sources = [(0, 0)] * (A + 2)
                sources[A], sources[roots[0]] = insertion, root
                if len(roots) == 2:
                    sources[roots[1]] = second
                sources[B] = tuple(-sum(row[j] for row in sources) for j in range(2))
                first_pair = bracket(root, insertion)
                if A == 1:
                    expected = -Fraction(first_pair, 2)
                elif len(roots) == 2:
                    second_pair = bracket(second, insertion)
                    roots_pair = bracket(root, second)
                    expected = (-Fraction(first_pair * second_pair, 6)
                                + Fraction((first_pair - second_pair) * roots_pair, 12)
                                + Fraction(roots_pair ** 2, 24))
                else:
                    expected = Fraction(first_pair ** 2, 24)
                    if forest[0][1] == "B":
                        expected = -expected
                with self.subTest(forest=forest, root=root):
                    self.assertEqual(ForestMoments(Graph(A + 2, edges), sources, (A, B)).moment(),
                                     expected)

    def test_all_topologies_through_five_loops(self):
        """All 100 topologies; every linear edge word and two mixed higher words."""
        total, comparisons = 0, 0
        for n in range(1, 6):
            order_count = 0
            for forest in forests(n):
                edges, A, B, roots, _ = graph(forest)
                sources = source_assignment(n, roots)
                form = chart(edges, A, B, sources, factored=True)
                engine = ForestMoments(Graph(n + 2, edges), sources, (A, B))
                words = [{}]
                words.extend({(edge_index, component): 1}
                             for edge_index in range(len(edges))
                             for component in range(2))
                words.extend((
                    {(0, 0): 2, (len(edges) // 2, 1): 1},
                    {(len(edges) - 2, 1): 2, (len(edges) - 1, 0): 2, (0, 1): 1},
                ))
                values = []
                for powers in words:
                    with self.subTest(loops=n, forest=forest, powers=powers):
                        reference = old_chart_moment(form, edges, powers)
                        value = engine.moment(powers)
                        self.assertEqual(value, reference)
                        self.assertEqual(engine.last_stats["global_ratio_coordinates"], 0)
                        self.assertEqual(engine.last_stats["current_interface_coordinates"], 6)
                        self.assertLessEqual(engine.last_stats["largest_local_ratio_count"], 4)
                        values.append(value)
                        comparisons += 1
                # All first-current moments obey the original vertex equations.
                bare = values[0]
                for vertex, source in enumerate(sources):
                    for component in range(2):
                        divergence = sum(
                            ((int(left == vertex) - int(right == vertex))
                             * values[1 + 2 * edge_index + component]
                             for edge_index, (left, right) in enumerate(edges)),
                            Fraction(0),
                        )
                        self.assertEqual(divergence, source[component] * bare)
                order_count += 1
                total += 1
            self.assertEqual(order_count, (1, 3, 8, 23, 65)[n - 1])
            print("exact original-edge word replay: loops", n,
                  "forests", order_count, flush=True)
        self.assertEqual(total, 100)
        self.assertEqual(comparisons, 2292)

    def test_source_translation_polynomials(self):
        for forest in forests(3):
            edges, A, B, roots, _ = graph(forest)
            source = source_assignment(3, roots)
            engine = ForestMoments(Graph(5, edges), source, (A, B))
            powers = {(0, 0): 1, (len(edges) - 1, 1): 2}
            polynomial = engine.source_polynomial(powers)
            self.assertTrue(all(sum(exponent) <= 6 for exponent in polynomial))
            for insertion in ((0, 0), (Fraction(3, 5), Fraction(-7, 4)), (-1, 4)):
                with self.subTest(forest=forest, insertion=insertion):
                    changed = source_assignment(3, roots, insertion)
                    evaluated = sum((coefficient * insertion[0] ** exponent[0]
                                     * insertion[1] ** exponent[1]
                                     for exponent, coefficient in polynomial.items()), Fraction(0))
                    self.assertEqual(evaluated, old_reference(edges, A, B, changed, powers))

    def test_oriented_words_and_original_edge_order(self):
        edges, A, B, roots, _ = graph(("(AB)", "(A(AB))"))
        sources = source_assignment(3, roots)
        G = Graph(5, edges)
        powers = {(0, 0): 1, (3, 1): 2, (6, 0): 1}
        reference = old_reference(edges, A, B, sources, powers)
        self.assertNotEqual(reference, 0)
        for reversed_edges in ((), (0,), (0, 6), (2, 3, 6)):
            changed = G.reverse_edges(reversed_edges)
            order = tuple(reversed(range(G.ne)))
            changed = changed.permute_edges(order)
            remapped = {(order.index(edge_index), component): power
                        for (edge_index, component), power in powers.items()}
            sign = (-1) ** (G.ne * (G.ne - 1) // 2
                            + sum(power for (edge_index, _), power in powers.items()
                                  if edge_index in reversed_edges))
            with self.subTest(reversed_edges=reversed_edges):
                self.assertEqual(ForestMoments(changed, sources, (A, B)).moment(remapped),
                                 sign * reference)

    def test_complete_fork_first_spinors(self):
        edges, A, B, roots, _ = graph(("((AB)(AB))",))
        for root, insertion in (((0, 1), (1, 0)), ((2, -1), (-3, 4)),
                                ((Fraction(2, 7), Fraction(-3, 5)),
                                 (Fraction(7, 3), Fraction(1, 11)))):
            sources = [(0, 0)] * 5
            sources[0], sources[A] = root, insertion
            sources[B] = tuple(-root[j] - insertion[j] for j in range(2))
            engine = ForestMoments(Graph(5, edges), sources, (A, B))
            bracket = root[0] * insertion[1] - root[1] * insertion[0]
            self.assertEqual(engine.moment(), 0)
            for endpoints in ((1, A), (1, B), (0, 1)):
                for component in range(2):
                    # Fraction division prevents even the integer fixtures
                    # from passing through binary floating-point arithmetic.
                    numerator = ({(1, A): 2 * root[component] - insertion[component],
                                  (1, B): -insertion[component] - 3 * root[component],
                                  (0, 1): -2 * insertion[component] - root[component]}
                                 [endpoints])
                    wanted = Fraction(bracket ** 3 * numerator, 8640)
                    self.assertEqual(engine.moment({(edges.index(endpoints), component): 1}), wanted)
                    twin = tuple(2 if vertex == 1 else vertex for vertex in endpoints)
                    self.assertEqual(engine.moment({(edges.index(twin), component): 1}), -wanted)
            self.assertEqual(engine.moment({(edges.index((A, B)), 0): 4}), 0)

    def test_known_fan_formula_at_higher_loops(self):
        for n in (6, 8, 12):
            edges, A, B, roots, _ = graph((fan_tree(n),))
            source = [(0, 0)] * (n + 2)
            source[roots[0]], source[A], source[B] = (0, 1), (1, 0), (-1, -1)
            engine = ForestMoments(Graph(n + 2, edges), source, (A, B))
            expected = Fraction((-1) ** ((n - 1) * (n - 2) // 2),
                                factorial(n) ** 3 * (n + 1))
            self.assertEqual(engine.moment(), expected)

    def test_known_repeated_fork_formula_at_higher_loops(self):
        for k in range(1, 5):
            G, sources, marked, internal, physical = fork_family(k)
            engine = ForestMoments(G, sources, marked)
            n = 3 * k
            phase = (-1) ** (n * (n + 1) // 2)
            self.assertEqual(engine.moment(dict.fromkeys(physical[:-1], 1)), 0)
            self.assertEqual(engine.moment(dict.fromkeys(physical, 1)),
                             Fraction(phase, 72 ** k * factorial(4 * k + 1)))
            self.assertEqual(engine.moment(dict.fromkeys(internal, 1)),
                             Fraction(phase * (-1) ** k, 36 ** k * factorial(4 * k + 1)))

    def test_asymmetric_branched_trees(self):
        alternating = "(AB)"
        for terminal in ("B", "A", "B"):
            alternating = "(" + terminal + alternating + ")"
        examples = (("(" + fan_tree(3) + fan_tree(2, "B") + ")",),
                    ("(" + fan_tree(4) + alternating + ")",))
        for forest in examples:
            edges, A, B, roots, _ = graph(forest)
            sources = source_assignment(A, roots)
            engine = ForestMoments(Graph(A + 2, edges), sources, (A, B))
            powers = {(0, 0): 1, (len(edges) // 2, 1): 1, (len(edges) - 2, 0): 1}
            with self.subTest(forest=forest):
                self.assertEqual(engine.moment(powers),
                                 old_reference(edges, A, B, sources, powers))

    def test_exact_budgets_and_input_validation(self):
        edges, A, B, roots, _ = graph(("((AB)(A(AB)))",))
        sources = source_assignment(A, roots)
        G = Graph(A + 2, edges)
        for budget in ({"max_terms": 1}, {"max_states": 1}):
            with self.assertRaises(ForestBudgetExceeded):
                ForestMoments(G, sources, (A, B), **budget).moment({(0, 0): 2})
        engine = ForestMoments(G, sources, (A, B))
        for powers in ({(0, 0): -1}, {(True, 0): 1}, {(0, 2): 1},
                       {(G.ne, 0): 1}, {(0, 0): Fraction(1, 2)}):
            with self.subTest(powers=powers), self.assertRaises(ValueError):
                engine.moment(powers)


if __name__ == "__main__":
    unittest.main(verbosity=2)
