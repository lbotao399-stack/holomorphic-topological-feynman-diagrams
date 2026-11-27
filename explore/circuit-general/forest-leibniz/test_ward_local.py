"""Independent local oracles for the finite-word forest Ward calculation.

Run: python3 explore/circuit-general/forest-leibniz/test_ward_local.py
"""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "code"))

from flint import fmpq
from htk import brute
from htk.graph import Graph
from htk.linalg import solve
from htk.spred import plan_generic, sp_reduce
from htk.treeflow import tree_flow
from ward_three import (EDGES, INNER_EDGES, GENERIC_SOURCES, PHYSICAL_SOURCES,
                        WardThree, _diamond_polynomial, boundary_value,
                        boundary_first_spinor, diamond_first_spinor)


def zero_resistance_currents(graph, sources, times):
    """Independent cycle solve, valid with one time exactly zero.

    This uses no vertex contraction or lifted-current formula. The cycle
    period matrix remains invertible at these contraction faces, so its
    exact value is the original graph's actual zero-resistance limit.
    """
    cycle = graph.cycle_matrix()
    particular = tree_flow(graph, sources)
    period = [[sum(cycle[i][edge] * times[edge] * cycle[j][edge]
                   for edge in range(graph.ne)) for j in range(graph.loops)]
              for i in range(graph.loops)]
    offset = [[sum(cycle[i][edge] * times[edge] * particular[edge][component]
                   for edge in range(graph.ne)) for component in range(2)]
              for i in range(graph.loops)]
    currents = solve(period, offset)
    return tuple(tuple(particular[edge][component]
                       - sum(cycle[i][edge] * currents[i][component] for i in range(graph.loops))
                       for component in range(2)) for edge in range(graph.ne))


class WardLocalChecks(unittest.TestCase):
    def test_complete_diamond_spider_chart_against_original_sp_density(self):
        words = ({}, {(0, 0): 1}, {(0, 1): 1}, {(0, 0): 1, (4, 1): 1},
                 {(2, 1): 2}, {(1, 0): 1, (2, 1): 1, (4, 0): 1})
        for P, Q, A in (((1, 2), (3, -1), (5, 3)),
                        ((F(2, 3), -3), (-1, F(4, 5)), (3, F(1, 7)))):
            sources = (P, Q, A, tuple(-sum(vector[i] for vector in (P, Q, A)) for i in range(2)))
            graph = Graph(4, INNER_EDGES)
            reduction = sp_reduce(graph, D=2, plan=plan_generic(graph))
            density = reduction.density(sources, T=0, method="generic", check=2)
            arguments = tuple(fmpq(F(value).numerator, F(value).denominator)
                              for row in (P, Q, A) for value in row)
            for word in words:
                polynomial = _diamond_polynomial(tuple((edge, component, power)
                                                        for (edge, component), power in sorted(word.items())),
                                                 1000000)
                expected = reduction.integrate(sources, T=0, powers=word, method="generic", check=0, dens=density)
                self.assertEqual(F(str(polynomial(*arguments))), expected)
            spinor = tuple(reduction.integrate(sources, T=0, powers={(0, component): 1},
                                                method="generic", check=0, dens=density)
                            for component in range(2))
            self.assertEqual(diamond_first_spinor(P, Q, A), spinor)

    def test_all_tree_face_current_limits_include_the_collapsed_edge(self):
        graph = Graph(6, EDGES)
        for sources in (GENERIC_SOURCES, PHYSICAL_SOURCES):
            with WardThree(sources, check=0) as ward:
                for quotient in ward.quotients:
                    position = quotient["position"]
                    times = [F(edge + 2, edge + 1) for edge in range(len(EDGES))]
                    times[position] = F(0)
                    expected = zero_resistance_currents(graph, sources, times)
                    remaining = times[:position] + times[position + 1:]
                    reduced = brute.currents(Graph(5, quotient["edges"]), remaining, quotient["sources"])
                    values = tuple(fmpq(value.numerator, value.denominator) for row in reduced for value in row)
                    lifted = tuple(tuple(F(str(polynomial(*values))) for polynomial in row)
                                   for row in quotient["lifted"])
                    self.assertEqual(lifted, expected)
                    self.assertTrue(brute.check_kirchhoff(graph, times, sources, expected))

    def test_universal_boundary_spinor_and_covariance(self):
        for sources in (GENERIC_SOURCES, PHYSICAL_SOURCES):
            value = boundary_first_spinor(sources)
            self.assertEqual(value, tuple(boundary_value(sources, {(2, component): 1})
                                          for component in range(2)))
            # The scalar density has determinant weight three; the inserted
            # current additionally transforms as the full original spinor.
            transform = ((1, 2), (3, 4))
            determinant = -2
            changed = tuple(tuple(sum(transform[i][j] * vector[j] for j in range(2))
                                  for i in range(2)) for vector in sources)
            expected = tuple(determinant ** 3 * sum(transform[i][j] * value[j] for j in range(2))
                             for i in range(2))
            self.assertEqual(boundary_first_spinor(changed), expected)


if __name__ == "__main__":
    unittest.main()
