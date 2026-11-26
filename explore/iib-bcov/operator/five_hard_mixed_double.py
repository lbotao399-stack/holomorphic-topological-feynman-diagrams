"""Five hard marks and three identical soft jumps on the genus-four two-double core."""

from fractions import Fraction as F
from itertools import combinations
from math import factorial
from pathlib import Path
import json

from quartic_seven_plane_double import graph_data, image, local_moment, polynomials, shapes, symmetries
from quartic_seven_plane_double_three_two import MARKS, integrate_polynomials
from fourfold_mixed_soft_ward import integrated_soft_placements
from point_graph import solve_graph


HARD = ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0),
        (0, 0, 0, 1), (-1, -1, -1, -1))
SOFT = (0, 1, 0, 1)
# The integral change (x,y,z,w) -> (y,x,z,w-y) puts SOFT on the first axis.
SLOPES = tuple(y for x, y, z, w in HARD)
NORMAL = tuple((x, z, w-y) for x, y, z, w in HARD)
CASES = (((0, 0, 4, 7, 7), (0, 2)), ((0, 1, 4, 7, 8), (0, 2)),
         ((0, 0, 0, 7, 7), (0, 1)))
SHAPE_TOTALS = (
    F(-26748313, 24841789378560000000), F(-4844761, 7452536813568000000),
    F(-6127859, 35488270540800000000), F(-10773253, 53232405811200000000),
    F(-4468973, 7452536813568000000), F(-19820357, 55894026101760000000),
    F(-2322697, 39924304358400000000), F(-1277, 24535100620800000),
    F(-247, 7360530186240000), F(-15331567, 74525368135680000000),
    F(-632561, 5081275100160000000))


def hard_order(pair):
    first = iter((1, 3))
    other = iter((0, 2, 4))
    return tuple(next(first if i in pair else other) for i in range(5))


def four_voltage(distance, pair):
    first = [6+i for i in pair]
    other = [6+i for i in range(5) if i not in pair]
    return (sum(distance(i, j) for i in first for j in other) - 3*distance(*first)
            - sum(distance(i, j) for i, j in combinations(other, 2)))


def verify_metric(data, ds, voltage, pair, times):
    graph, _, u, adj, _, _, _ = data
    order = hard_order(pair)
    rays, _ = solve_graph(graph.edges, [(0,)*4]*6+[HARD[i] for i in order], times, 0)
    bare = integrated_soft_placements(rays, graph.C, times, SOFT)/factorial(3)
    U = int(u(*times))
    values = [int(d(*times)) for d in ds]
    square = sum(values[i]*int(adj[i][j](*times))*values[j] for i in range(4) for j in range(4))
    assert bare == F(square, U**8)

    # Inject a unit current at each hard vertex and absorb it at vertex zero.
    # Columns of this grounded Green function give all pair resistances.
    sources = [[int(i == 6+j)-int(i == 0) for j in range(5)] for i in range(11)]
    _, potential = solve_graph(graph.edges, sources, times, 0)
    def resistance(i, j):
        return (potential[i][i-6]+potential[j][j-6]
                - potential[i][j-6]-potential[j][i-6])
    K = four_voltage(resistance, pair)/4
    assert int(voltage(*times)) == 4*U*K
    assert bare*K/U**4 == F(square*int(voltage(*times)), 4*U**13)


def chamber_polynomials(data, pair):
    _, ctx, _, _, dist, _, _ = data
    raw = four_voltage(lambda i, j: dist[i][j], pair)
    assert all(c.denominator == 1 for _, c in raw.terms())
    voltage = ctx.from_dict({m: int(c) for m, c in raw.terms()})
    ds, _ = polynomials(data, hard_order(pair), normal=NORMAL, slopes=SLOPES,
                        polarization=(0, 0, 0))
    for times in (tuple(range(1, 15)), tuple(2+(i*i+3*i)%17 for i in range(14))):
        verify_metric(data, ds, voltage, pair, times)
    return ds, voltage


def chamber(data, pair):
    return integrate_polynomials(data, *chamber_polynomials(data, pair))/4


def color_orbit(marks, pair):
    config = tuple(tuple(int(i in pair) for i, edge in enumerate(marks) if edge == k)
                   for k in range(9))
    return {image(config, action) for action in symmetries()}


def complete_shapes():
    # The three hard current rows on a bubble live in the time-orthogonal
    # subspace of its segments. Hence each bubble needs at least two hard marks.
    return [marks for marks, _ in shapes()]+list(MARKS)


def complete_shape(marks):
    data = graph_data(marks)
    unseen = set(combinations(range(5), 2))
    total = F()
    count = 0
    while unseen:
        pair = min(unseen)
        orbit = color_orbit(marks, pair)
        for config in orbit:
            projected = tuple(k for k, path in enumerate(config) for _ in path)
            if projected == marks:
                unseen.discard(tuple(i for i, color in enumerate(sum(config, ())) if color))
        value = chamber(data, pair)
        weight = F(12*len(orbit), len(symmetries()))
        count += weight
        total += weight*value
        print('color', marks, pair, weight, value, flush=True)
    shape_size = len({tuple(sorted(action[k][0] for k in marks)) for action in symmetries()})
    assert count == F(120*shape_size, len(symmetries()))
    local_moment.cache_clear()
    return total


def recorded_shapes(topology, marks):
    from quartic_scalar_response import interval

    data = json.loads(Path(__file__).with_name('five_hard_mixed_periods.json').read_text())[topology]
    rows = data['shape_orbits']
    assert [tuple(r['marks']) for r in rows] == list(marks)
    values = [tuple(map(F, r['channels'])) for r in rows]
    assert all(interval(value)[1] < 0 for value in values)
    assert tuple(map(sum, zip(*values))) == tuple(map(F, data['total']))
    return values


if __name__ == '__main__':
    from sys import argv

    if '--all' in argv:
        total = F()
        for marks, expected in zip(complete_shapes(), SHAPE_TOTALS, strict=True):
            value = complete_shape(marks)
            assert value == expected
            total += value
            print('shape', marks, value, flush=True)
        assert total == F(-22550149, 6387888697344000000)
        print('two-double total', total)
    else:
        for case in map(int, argv[1:] or [0]):
            marks, pair = CASES[case]
            value = chamber(graph_data(marks), pair)
            multiplicity = F(2*6*len(color_orbit(marks, pair)), len(symmetries()))
            print(marks, pair, 'chamber', value, 'multiplicity', multiplicity,
                  'orbit', multiplicity*value, flush=True)
