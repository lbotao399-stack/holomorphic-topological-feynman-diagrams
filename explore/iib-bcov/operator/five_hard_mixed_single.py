"""Five-hard mixed Ward periods of the single-double four-loop core."""

from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement

from flint import fmpq_mpoly_ctx
from five_hard_mixed_double import chamber_polynomials, recorded_shapes
from quartic_seven_plane_current_basis import expand, integer_poly
from quartic_seven_plane_double import local_moment
from quartic_seven_plane_single import cuts, graph_data
from quartic_seven_plane_single_projection import project_currents, shape_orbits


def shape_weights(marks, *, fold=True):
    data = graph_data(marks)
    ctx = fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
    bases = [[] for _ in range(4)]
    voltages = {}
    for pair in combinations(range(5), 2):
        ds, voltage = chamber_polynomials(data, pair)
        rows = []
        for i, d in enumerate(ds):
            rows.extend((i, j, c) for j, c in expand(bases[i], ctx.from_dict(d.to_dict())))
        v = ctx.from_dict(voltage.to_dict())
        for (i, j, c), (k, h, b) in combinations_with_replacement(rows, 2):
            key = (i, j, k, h)
            # Twelve labelings per color placement, and voltage = 4 U K.
            voltages[key] = voltages.get(key, 0) + 3*(1+((i, j) != (k, h)))*c*b*v
    prepared = {(i, j): integer_poly(p, data[1])
                for i, row in enumerate(bases) for j, (_, p) in enumerate(row)}
    return project_currents(data, prepared, voltages, fold=fold)


if __name__ == '__main__':
    from sys import argv

    shapes = list(shape_orbits())
    recorded = recorded_shapes('single_double', [m for m, _ in shapes])
    if '--star-check' in argv:
        from quartic_seven_plane_single_star import star_integral
        for index in (0, 4):
            marks, size = shapes[index]
            assert (F(size, 8)*star_integral(shape_weights(marks, fold=False)), 0) == recorded[index]
        print('Two independent star-triangle integrations.')
        raise SystemExit
    selected = list(range(len(shapes))) if '--all' in argv else list(map(int, argv[1:] or [0]))
    total = [F(), F()]
    for index in selected:
        marks, size = shapes[index]
        channels = tuple(F(size, 8)*x for x in cuts(shape_weights(marks)))
        assert channels == recorded[index]
        print('shape', index, marks, *channels, flush=True)
        total = [x+y for x, y in zip(total, channels)]
        graph_data.cache_clear()
        local_moment.cache_clear()
    print('selected shape total', *total, flush=True)
