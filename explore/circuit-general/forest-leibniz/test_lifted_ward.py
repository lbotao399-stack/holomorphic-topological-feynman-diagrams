"""Independent exact controls for the higher-loop forest Ward lift.

Run: python3 explore/circuit-general/forest-leibniz/test_lifted_ward.py
"""
from fractions import Fraction as F
import unittest

from lifted_ward import (LiftedWard, lifted_fixture, face_census, inner_boundary_spinor,
                         boundary_spinor, reconstruct_last, _fraction,
                         GENERIC_SOURCES, PHYSICAL_SOURCES)
from ward_three import diamond_first_spinor, boundary_first_spinor
from htk.graph import Graph
from htk.forest_moments import ForestMoments


class LiftedWardChecks(unittest.TestCase):
    def test_independent_cubic_boundary_moments_through_eighteen_loops(self):
        # This tree-message engine does not use ForestContractions,
        # CoreCurrentPeriods, the factorial comb transfer, or our beta sum.
        P, Q, A = (F(2, 3), F(-3)), (F(-1), F(4, 5)), (F(3), F(1, 7))
        for pairs in (0, 1, 2, 3, 5, 8):
            original, _, _, word = lifted_fixture(pairs)
            vertices = tuple(v for v in range(original.nv) if v not in (0, 5))
            inside = tuple(e for e, (a, b) in enumerate(original.edges) if a in vertices and b in vertices)
            graph = Graph(len(vertices), [(vertices.index(original.edges[e][0]),
                                          vertices.index(original.edges[e][1])) for e in inside])
            rows = [(F(0), F(0))] * graph.nv
            for vertex, row in ((1, P), (2, Q), (3, A),
                                (4, tuple(-P[a] - Q[a] - A[a] for a in range(2)))):
                rows[vertices.index(vertex)] = row
            engine = ForestMoments(graph, rows, (vertices.index(3), vertices.index(4)), max_terms=1000000)
            edge = inside.index(next(iter(word))[0])
            actual = tuple(engine.moment({(edge, a): 1}) for a in range(2))
            expected = tuple(_fraction(value) for value in inner_boundary_spinor(pairs, P, Q, A))
            self.assertEqual(actual, expected)

    def test_base_formula_and_complete_spinor_covariance(self):
        P, Q, A = (F(2, 3), -3), (-1, F(4, 5)), (3, F(1, 7))
        self.assertEqual(tuple(_fraction(v) for v in inner_boundary_spinor(0, P, Q, A)),
                         diamond_first_spinor(P, Q, A))
        for sources in (GENERIC_SOURCES, PHYSICAL_SOURCES):
            self.assertEqual(boundary_spinor(0, sources), boundary_first_spinor(sources))
            transform = ((1, 2), (3, 4))
            changed = tuple(tuple(sum(transform[a][b] * row[b] for b in range(2))
                                  for a in range(2)) for row in sources)
            for pairs in (1, 2, 8):
                value = boundary_spinor(pairs, sources)
                expected = tuple((-2) ** (3 + 2 * pairs) * sum(transform[a][b] * value[b] for b in range(2))
                                 for a in range(2))
                self.assertEqual(boundary_spinor(pairs, changed), expected)

    def test_complete_five_and_seven_loop_zero_face_witnesses(self):
        for pairs, expected_count in ((1, 9), (2, 21)):
            census = face_census(pairs)
            self.assertEqual(len(census["proper_faces"]), expected_count)
            nonzero = [row for row in census["tree_faces"] if row["zero_witness"] is None]
            self.assertEqual(tuple(row["position"] for row in nonzero), (0, 1, 4))
            self.assertEqual(tuple(row["sign"] for row in nonzero), (1, -1, 1))
            self.assertTrue(any(row["zero_witness"] is not None
                                and len(row["zero_witness"]["vertices"]) > 2
                                for row in census["tree_faces"]))
            for row in census["tree_faces"]:
                witness = row["zero_witness"]
                if witness is None:
                    continue
                vertices = set(witness["vertices"])
                count = sum(set(edge) <= vertices for edge in row["quotient_edges"])
                self.assertEqual(count, witness["edge_count"])
                self.assertGreater(count, 2 * len(vertices) - 3)

    def test_actual_five_loop_spinor_and_period_reconstruction(self):
        expected = (F(-1979627, 4536000), F(165839, 113400))
        with LiftedWard(1, PHYSICAL_SOURCES) as ward:
            for component in range(2):
                result = ward.verify(component)
                self.assertEqual(result["rational_boundary"], expected[component])
                first, middle, last = (term["coefficients"] for term in result["terms"])
                self.assertNotEqual(first["z3"], 0)
                self.assertEqual(first["z3"], -last["z3"])
                self.assertNotIn("z3", middle)
                self.assertEqual(ward.reconstruct_last(component, first, middle), last)
        # The helper retains every supplied period label, not merely zeta(3).
        self.assertEqual(reconstruct_last({"1": F(2, 3), "z3": F(-7, 11), "z5": F(5, 17)},
                                          F(4, 9), F(2, 5)),
                         {"1": F(8, 45), "z3": F(7, 11), "z5": F(-5, 17)})
        with self.assertRaises(ValueError):
            reconstruct_last(F(1), F(2), F(3), signs=(1.0, -1, 1))


if __name__ == "__main__":
    unittest.main()
