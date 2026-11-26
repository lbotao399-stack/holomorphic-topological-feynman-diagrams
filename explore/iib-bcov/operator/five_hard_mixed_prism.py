"""Five-hard mixed Ward periods on the five prism placements with two unmarked stars."""

from fractions import Fraction as F
from itertools import combinations
from math import prod

from flint import fmpq, fmpq_mat
from five_hard_mixed_double import NORMAL, SLOPES, chamber_polynomials, four_voltage, hard_order, recorded_shapes
from fourfold_bulk import symmetries
from quartic_seven_plane_prism import CORE, MARKS, cut_weights, cuts, graph_data, project_currents
from quartic_seven_plane_prism_paired import MARKS as PAIRED, project as paired_project


def projected(marks, pairs):
    data = graph_data(marks)
    if marks in MARKS:
        def terms():
            for pair in pairs:
                ds, four = chamber_polynomials(data, pair)
                yield ds, 3*four
        return project_currents(marks, data, terms())
    else:
        # Twelve labelings per color and the voltage divisor four give factor three.
        edges = list(combinations(range(5), 2))
        table = fmpq_mat([[3*four_voltage(lambda x, y: int(tuple(sorted((x-6, y-6))) == edge), pair)
                           for edge in edges] for pair in pairs])
        return paired_project(marks, [hard_order(pair) for pair in pairs], NORMAL,
                              slopes=SLOPES, table=table)


def shape_channels(marks):
    integral = projected(marks, list(combinations(range(5), 2)))
    actions = symmetries(CORE)
    images = {tuple(sorted(action[k] for k in marks)) for action in actions}
    return tuple(F(len(images), len(actions))*x for x in cuts(cut_weights(integral)))


def angular_check(marks, pair, times):
    from quartic_seven_plane_lift_positions import CTX, data as electric, sampled
    from quartic_seven_plane_prism_check import chamber, dualize, kirchhoff, pair_graph

    u, adj, cycles, intercept, slope, green = electric(CORE, times, marks)
    order = hard_order(pair)
    d = sampled([NORMAL[i] for i in order], [SLOPES[i] for i in order], marks,
                cycles, intercept, slope)
    square = sum((adj[i][j]*d[i]*d[j] for i in range(4) for j in range(4)), CTX.constant(0))
    s = CTX.gens()

    def kernel(i, j):
        a, b = CORE[marks[i]]
        c, e = CORE[marks[j]]
        value = (int(u*green[a][c])*(1-s[i])*(1-s[j]) + int(u*green[a][e])*(1-s[i])*s[j]
                 + int(u*green[b][c])*s[i]*(1-s[j]) + int(u*green[b][e])*s[i]*s[j])
        if marks[i] == marks[j]:
            value += u*times[marks[i]]*(s[min(i, j)]-s[i]*s[j])
        return value

    def distance(i, j):
        i, j = i-6, j-6
        return kernel(i, i)+kernel(j, j)-2*kernel(i, j)

    numerator = 3*four_voltage(distance, pair)*square
    integrated = sum(int(c)*chamber(tuple(map(int, powers)), marks) for powers, c in numerator.terms())
    measured = integrated*prod(times[i] for i in marks)/(6**10*u**21*sum(times)**2)
    dual, jac = dualize(times)
    dual_u = kirchhoff(pair_graph(), dual)
    target = projected(marks, [pair])(*(fmpq(str(x)) for x in dual))
    assert measured and measured*jac == F(str(target))/dual_u**13


if __name__ == '__main__':
    from sys import argv

    shapes = (*MARKS, *PAIRED)
    recorded = recorded_shapes('two_unmarked_stars_on_prism', shapes)
    if '--check-angular' in argv:
        for i in (0, 3):
            angular_check(shapes[i], (0, 2), (2, 3, 5, 7, 11, 13, 17, 19, 23))
        print('Two independent ordered-position and star Jacobian checks.')
        raise SystemExit
    indices = range(len(shapes)) if '--all' in argv else map(int, argv[1:] or [0])
    total = [F(), F()]
    for i in indices:
        value = shape_channels(shapes[i])
        assert value == recorded[i]
        total = [x+y for x, y in zip(total, value)]
        print('prism', i, shapes[i], *value, flush=True)
    print('selected shapes', *total, flush=True)
