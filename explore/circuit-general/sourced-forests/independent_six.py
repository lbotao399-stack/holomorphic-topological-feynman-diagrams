"""Independent bordered-determinant and face-dual replay of the six-loop forest.

This does not call ForestContractions, CoreCurrentPeriods, per5, or exact
division by the Kirchhoff polynomial. The local sign and surviving bracket
are specified explicitly, and the numerator comes directly from its bordered
determinant before a separate face-dual prism-table contraction.
"""

from fractions import Fraction as F
from flint import fmpq, fmpq_mpoly_ctx

from independent_duality import (EDGES, SOURCES, Graph, exact_original_numerator,
                                replay, tree_flow)

EXPECTED_SIX = {"1": -F(2553462100217, 369515520000),
                "z3": -F(107288658589, 1710720000),
                "z5": F(7666393, 95040)}


def determinant(matrix, context):
    states = {0: context.constant(1)}
    for row in matrix:
        following = {}
        for mask, coefficient in states.items():
            for column, value in enumerate(row):
                if value and not mask & (1 << column):
                    target = mask | (1 << column)
                    term = coefficient * value * (-1) ** ((mask >> (column + 1)).bit_count())
                    following[target] = following.get(target, context.constant(0)) + term
        states = following
    return states.get((1 << len(matrix)) - 1, context.constant(0))


def inserted_numerator():
    graph = Graph(7, EDGES)
    context = fmpq_mpoly_ctx.get(tuple("t%d" % edge for edge in range(graph.ne)))
    times = context.gens()
    cycles = graph.cycle_matrix()
    flow = tree_flow(graph, SOURCES)
    matrix = [[sum((cycles[i][edge] * cycles[j][edge] * times[edge]
                    for edge in range(graph.ne)), context.constant(0))
               for j in range(graph.loops)] for i in range(graph.loops)]
    rhs = [[sum((cycles[i][edge] * int(flow[edge][component]) * times[edge]
                 for edge in range(graph.ne)), context.constant(0))
            for component in range(2)] for i in range(graph.loops)]
    augmented = [row + rhs[i] for i, row in enumerate(matrix)]
    for edge in (4, 10):  # retained edges (0,5) and (5,6), in EDGES order
        augmented.append([context.constant(cycles[i][edge]) for i in range(graph.loops)]
                         + [context.constant(int(value)) for value in flow[edge]])
    bracket_numerator = determinant(augmented, context)
    assert bracket_numerator and all(sum(exponent) == 5 for exponent, _ in bracket_numerator.terms())

    # The six-loop single move gives -<y_05,y_56>/4 in EDGES order:
    # its local permutation is odd; the retained 65 arrow is reversed;
    # restoring EDGES from (06,05,65,01,13,02,24,35,36,45,46) is odd.
    # Thus three minus signs multiply. No production contraction is called.
    bare = context.from_dict({exponent: fmpq(coefficient)
                              for exponent, coefficient in exact_original_numerator().items()})
    inserted = -bare * bracket_numerator / 4
    return {tuple(map(int, exponent)): F(str(coefficient))
            for exponent, coefficient in inserted.terms()}


def replay_six():
    numerator = inserted_numerator()
    value, details = replay(numerator, q=8, expected=EXPECTED_SIX)
    assert details["numerator_terms"] == 287184
    assert details["core_period_orbits"] == 7784
    return value, details


if __name__ == "__main__":
    value, details = replay_six()
    print(details)
    print({label: str(coefficient) for label, coefficient in value.items()})
