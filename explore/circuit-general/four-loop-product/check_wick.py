"""Independent literal contraction and full-source reconstruction controls."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import factorial
from pathlib import Path
import json
import unittest

import quartic_b_squared as calculation


def literal_differential_sum():
    """Enumerate graphs first and apply Delta, without complete matching search."""
    result = defaultdict(int)
    all_edges = tuple(edge for edge in combinations(range(6), 2) if edge != (4, 5))
    graph_count = 0
    for physical in combinations(all_edges, 8):
        edges = physical + ((4, 5),)
        degrees = tuple(sum(vertex in edge for edge in edges) for vertex in range(6))
        if degrees[4:] != (3, 4) or degrees[0] not in (3, 4) or any(degree not in (2, 3) for degree in degrees[1:4]):
            continue
        if any(sum(a in vertices and b in vertices for a, b in edges) > 2 * size - 3
               for size in range(2, 7) for vertices in combinations(range(6), size)):
            continue
        graph_count += 1
        for flavors in product(range(1, 5), repeat=3):
            # Independent construction of the original even/odd monomial.
            original = tuple(('gamma', flavor, 0, 0) for flavor in range(1, 5))
            for vertex in (1, 2, 3):
                flavor = flavors[vertex - 1]
                original += (('c', 0, vertex, 0), ('gamma', flavor, vertex, 0), ('beta', flavor, vertex, 0))
            original += (('beta', 1, 4, 0), ('beta', 2, 4, 0), ('gamma5', 5, 4, 0),
                         ('b', 0, 5, 0), ('b', 0, 5, 1), ('beta', 3, 5, 0))
            states = {(original, None): 1}
            for edge_index in reversed(range(8)):
                left_vertex, right_vertex = physical[edge_index]
                following = defaultdict(int)
                for (remaining, current_edge), coefficient in states.items():
                    for left in (field for field in remaining if field[2] == left_vertex):
                        for right in (field for field in remaining if field[2] == right_vertex):
                            species = (left[0], right[0])
                            if species in (('b', 'c'), ('c', 'b')):
                                sign = 1
                            elif species in (('beta', 'gamma'), ('gamma', 'beta')) and left[1] == right[1]:
                                sign = -1
                            else:
                                continue
                            after = list(remaining)
                            for field in (right, left):
                                if field[0] in ('beta', 'c'):
                                    sign *= (-1) ** sum(old[0] in ('beta', 'c') for old in after[:after.index(field)])
                                after.remove(field)
                            inserted = edge_index if ('beta', 2, 4, 0) in (left, right) else current_edge
                            following[tuple(after), inserted] += sign * coefficient
                states = {key: value for key, value in following.items() if value}
                if not states:
                    break
            for (remaining, current_edge), coefficient in states.items():
                result[edges, flavors, remaining, current_edge] += coefficient
    return {key: value for key, value in result.items() if value}, graph_count


class FourLoopChecks(unittest.TestCase):
    def test_every_complete_matching_against_literal_delta_operators(self):
        records, rejected = calculation.wick_matchings()
        self.assertEqual((len(records), len(rejected)), (324, 732))
        expected = defaultdict(int)
        for record in records:
            expected[record['ordered_edges'], record['matter_flavors'], record['surviving_fields'],
                     record['current_edge']] += record['wick_sign']
        actual, graph_count = literal_differential_sum()
        self.assertEqual(actual, {key: value for key, value in expected.items() if value})
        self.assertGreater(graph_count, 0)
        # Two original b occurrences give the derivative multiplicity of b**2.
        self.assertTrue(all(abs(coefficient) == 2 for coefficient in actual.values()))

    def test_all_W_placements_and_source_rows(self):
        groups = calculation.grouped_matchings()[0]
        for slot in (1, 2, 3):
            self.assertEqual(calculation.grouped_matchings(slot)[0], groups)
        self.assertEqual(len(groups), 14)
        for edges, locations in groups:
            vectors = ((1, 2), (3, -1), (-2, 4))
            rows = calculation.source_rows(locations, vectors)
            self.assertEqual(tuple(map(sum, zip(*rows))), (0, 0))
            for vertex in range(5):
                self.assertEqual(rows[vertex], tuple(sum(vectors[j][a] for j in range(3) if locations[j] == vertex)
                                                     for a in range(2)))

    def test_complete_covariant_reconstruction(self):
        data = json.loads(calculation.CERTIFICATE.read_text())
        calculation.verify_saved(data)
        self.assertEqual(len(data['complete_covariant_charge_polynomial']), 16)
        self.assertTrue(any(label == 'z3' for sample in data['samples'] for graph in sample['graphs']
                            for spinor in graph['full_first_current_moments'].values()
                            for entry in spinor for label in entry))
        for _, channels in data['complete_covariant_charge_polynomial']:
            self.assertEqual(set(channels), {'1'})
        generic_graphs = data['independent_full_source_replays'][0]['graphs']
        self.assertEqual(sum(any('z3' in component for spinor in graph['full_first_current_moments'].values()
                                 for component in spinor) for graph in generic_graphs), 7)
        # Exact covariant identities, with no numerical source specialization.
        diagonal = defaultdict(Fraction)
        spectator_zero = defaultdict(Fraction)
        for vector, powers, value in data['complete_physical_covariant_spinor']['1']:
            if powers[2] == 0:
                diagonal[0 if vector == 2 else vector] += Fraction(value)
            if vector != 2 and powers[1] == powers[2] == 0:
                spectator_zero[vector] += Fraction(value)
        self.assertEqual(dict(diagonal), {0: Fraction(1, 960), 1: Fraction(-1, 3200)})
        self.assertFalse(any(spectator_zero.values()))

    def test_normalization_and_exact_source_guards(self):
        # Factor i is explicit; all remaining quantities are exact rationals.
        raw_without_i = (-1) ** 4 * Fraction(4, factorial(4)) * (-1) ** 10 * (-2)
        self.assertEqual(raw_without_i, Fraction(-1, 3))
        self.assertEqual(raw_without_i * Fraction(-1, 2) ** 9, Fraction(1, 1536))
        self.assertEqual(sum((1, -1, 0, 0)), 0)
        self.assertEqual(sum(value ** 3 for value in (1, -1, 0, 0)), 0)
        with self.assertRaises(TypeError):
            calculation.source_rows((0, 1, 4), ((1.0, 0), (0, 1), (1, 1)))
        with self.assertRaises(ValueError):
            calculation.wick_matchings(4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
