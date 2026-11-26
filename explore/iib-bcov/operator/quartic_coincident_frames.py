"""Exact seven-planar interpolation on conserved normal waves with coincident Gale weights."""
from fractions import Fraction
from itertools import combinations
from math import prod

from flint import fmpq_mat, fmpq_mpoly_ctx

from planar_normal_tensors import GENERATORS, generators, powers, tensor_row
from quartic_seven_plane_space import (
    as_matrix, evaluation_row, gale, normal_frame, plane_rows, standard_lifts,
)


FRAMES = (
    (0, 1, 2, 0, 1), (0, 1, 3, 0, 1), (0, 1, 4, 0, 1),
    (0, 1, 2, -1, 0), (0, 1, 3, -1, 0),
)
CONSTANT = ((0, (0, 1)), (1, (0, 1)), (2, (0, 1)),
            (3, (0, 1)), (3, (0, 2)), (4, (0, 2)))
LINEAR = CONSTANT[:4]
QUADRATIC = (CONSTANT[0], CONSTANT[1], CONSTANT[3])
PLANES = tuple(((1, 0),) * 3 + ((0, 1),) * 3 + (last,)
               for last in ((1, 0), (1, 1), (2, 2)))


def measurements():
    """Arguments of evaluation_row, in reconstruction order."""
    def node(frame, pair, plane, sign):
        z = FRAMES[frame]
        return plane, normal_frame(z), standard_lifts(plane, z, (sign, 0)), pair

    nodes = [node(frame, pair, plane, 0) for frame, pair in CONSTANT for plane in PLANES]
    for frame, pair in LINEAR:
        nodes.append(node(frame, pair, PLANES[1], 1))
        if (frame, pair) in QUADRATIC:
            nodes.append(node(frame, pair, PLANES[1], -1))
    return nodes


def restriction(degree):
    """Coefficient matrix on the full hyperplane z[0]=z[4] for test pair (0,1)."""
    ctx = fmpq_mpoly_ctx.get(('x', 'y'))
    x, y = ctx.gens()
    z = (ctx.constant(0), ctx.constant(1), x, y, ctx.constant(0))
    centered = tuple(t - sum(z) / 5 for t in z)
    elementary = [sum(prod(c) for c in combinations(centered, k)) for k in range(2, 6)]
    columns = [(v * prod(e**a for e, a in zip(elementary, exponent))).to_dict()
               for (r, s), v in zip(GENERATORS, generators(centered, (0, 1)))
               for exponent in powers(degree - 1 - r - s)]
    monomials = sorted(set().union(*columns))
    return fmpq_mat([[column.get(key, 0) for column in columns] for key in monomials])


def reconstruct(values):
    """Recover the 25 coefficients in the original evaluation_row basis."""
    if len(values) != 25:
        raise ValueError("Expected 25 seven-planar measurements")
    return as_matrix([evaluation_row(*node) for node in measurements()]).solve(
        as_matrix([[v] for v in values]))


def check():
    for z in FRAMES:
        q = normal_frame(z)
        assert all(sum(component) == 0 for component in zip(*q))
        assert gale(q) == tuple(v-z[-1] for v in z)
        assert all(v.denominator == 1 for p in q for v in p)
    for m, keys, expected in ((7, QUADRATIC, '-1728/125'),
                              (8, LINEAR, '1410048/3125'),
                              (9, CONSTANT, '216710774784/1953125')):
        rows = [tensor_row(FRAMES[frame], pair, m) for frame, pair in keys]
        assert str(as_matrix(rows).det()) == expected
    nodes = measurements()
    matrix = as_matrix([evaluation_row(*node) for node in nodes])
    assert matrix.rank() == 25
    # Independently evaluate the four-generator basis on both sets of nodes.
    def generator_row(node):
        plane, q, shifts, pair = node
        z = gale(q)
        mean = Fraction(sum(z), 5)
        lift = tuple(sum((v-mean)*b for v, b in zip(z, column)) for column in zip(*shifts))
        return [p * t for k, planar in enumerate(plane_rows(plane, lift))
                for p in planar for t in tensor_row(z, pair, 9-k)]

    from quartic_seven_plane_space import interpolation
    old_nodes, old_inverse = interpolation()
    change = matrix.solve(as_matrix([generator_row(node) for node in nodes]))
    assert change.det() and old_inverse * as_matrix([generator_row(node) for node in old_nodes]) == change
    check_column = as_matrix([[i**3 - 3*i + 1] for i in range(25)])
    assert reconstruct((matrix * check_column).entries()) == check_column

    for degree in range(4, 16):
        table = restriction(degree)
        assert table.ncols()-table.rank() == int(degree == 15)
    # The first lost direction is nonzero on a simple (all distinct) normal frame.
    table = restriction(15)
    outside = tensor_row((0, 1, 2, 3, 5), (0, 1), 15)
    assert as_matrix(table.tolist() + [outside]).rank() == 24
    print("Sparse interpolation rank 25. Coincidence restriction injective through degree 14.")
    print("Degree 15: rank 23 of 24, restored by the five-distinct-weight test.")


if __name__ == '__main__':
    check()
