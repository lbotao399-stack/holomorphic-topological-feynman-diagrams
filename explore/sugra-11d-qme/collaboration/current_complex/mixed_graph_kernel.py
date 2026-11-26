"""Nine-term determinant of a theta graph on the positive-jet mixed-current probe."""

from functools import cache
from itertools import combinations, permutations
from math import factorial, prod

from block_homotopy import transvectant
from mixed_planar import graph_test
from math import comb
from fractions import Fraction


def pure(degree, p):
    x, y = p
    return {
        (j, degree - j): comb(degree, j) * x**j * y ** (degree - j)
        for j in range(degree + 1)
        if x**j * y ** (degree - j)
    }


def invariants(forms):
    a, b, c, d = forms
    tr = transvectant

    def scalar(g):
        return g.get((0, 0), 0)

    return (
        scalar(tr(a, b, 3)) * scalar(tr(c, d, 4)),
        scalar(tr(tr(a, b, 1), tr(c, d, 2), 4)),
    )


def verify_planar():
    target = [0, 0]
    for word, c in graph_test().items():
        forms = [{(i, h + 2 - i): 1} for _, h, i in word]
        for j, v in enumerate(invariants(forms)):
            target[j] += c * v
    assert target == [-1036800, -8709120]
    values = [
        invariants([pure(n, p) for n, p in zip((3, 3, 4, 4), momenta)])
        for momenta in PLANAR
    ]
    assert values == [(20736, 124416), (331776, -995328)]
    assert [
        Fraction(-190, 3) * values[0][j] + Fraction(5, 6) * values[1][j]
        for j in range(2)
    ] == target


PLANAR = (((1, 0), (0, 1), (1, 0), (0, 1)), ((1, 0), (0, 1), (1, 1), (1, -1)))
SHAPES = ((1, 6, 6), (2, 5, 6), (3, 4, 6), (3, 5, 5), (4, 4, 5))


def det(columns):
    n = len(columns)
    rows = [list(x) for x in zip(*columns)]
    sign = prev = 1
    for i in range(n - 1):
        if not rows[i][i]:
            j = next((j for j in range(i + 1, n) if rows[j][i]), None)
            if j is None:
                return 0
            rows[i], rows[j] = rows[j], rows[i]
            sign = -sign
        a = rows[i][i]
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                rows[j][k] = (a * rows[j][k] - rows[j][i] * rows[i][k]) // prev
        prev = a
    return sign * rows[-1][-1]


def lift(p, path):
    z = (0,) * len(p)
    return (z if path == 2 else p) + (z if path == 0 else p)


@cache
def small(directions, paths):
    return det(tuple(lift(p, c) for p, c in zip(directions, paths)))


def directions(faxis):
    vectors = ((1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1))
    function = tuple((2 if faxis == 3 else 1) * int(i == faxis % 3) for i in range(3))
    return vectors, function


def coefficients(planar, vertical, f, paths):
    # Vertex labels in the test order H1,H2,H3,mark,Y1,...,Y5,gamma.
    h0 = tuple(paths[:3])
    mark = paths[3]
    ys = tuple(paths[4:9])
    gamma = paths[9]
    on_mark = small(planar, h0 + (mark,))
    vertical_mark = small((f, *vertical), (mark, *ys))
    rx = tuple(
        -vertical_mark * small((*planar[:3], p), h0 + (gamma,)) for p in planar[:3]
    )
    ry = tuple(-on_mark * small((*vertical, p), ys + (gamma,)) for p in vertical)
    last = -on_mark * small((*vertical, f), ys + (gamma,)) - vertical_mark * small(
        planar, h0 + (gamma,)
    )
    return (*rx, *ry, last)


def paths(shape):
    positions = tuple(range(10))
    for first in combinations(positions, shape[0] - 1):
        rest = tuple(i for i in positions if i not in first)
        for second in combinations(rest, shape[1] - 1):
            yield tuple(0 if i in first else 1 if i in second else 2 for i in positions)


def direct(planar, vertical, f, path, x, y, a, b):
    hz = (0, 0)
    vz = (0, 0, 0)
    cols = [tuple(s * m for m in p) + vz for s, p in zip(x, planar[:3])]
    cols += [tuple(a * m for m in planar[3]) + tuple(b * m for m in f)]
    cols += [hz + tuple(s * m for m in p) for s, p in zip(y, vertical)]
    cols.append(tuple(-sum(v) for v in zip(*cols)))
    return det(tuple(lift(p, c) for p, c in zip(cols, path)))


def check(shape, configs):
    counts = [0] * len(configs)
    union = 0
    x = (1, 2, -1)
    y = (1, 2, 3, -2, -1)
    a = -2
    b = 3
    for i, path in enumerate(paths(shape)):
        any_nonzero = False
        for j, (planar, vertical, f) in enumerate(configs):
            c = coefficients(planar, vertical, f, path)
            if any(c):
                counts[j] += 1
                any_nonzero = True
            if i % 71 == 0:
                predicted = (
                    prod(x)
                    * prod(y)
                    * (
                        b * sum(u * v for u, v in zip(x, c[:3]))
                        + a * sum(u * v for u, v in zip(y, c[3:8]))
                        + a * b * c[-1]
                    )
                )
                assert predicted == direct(planar, vertical, f, path, x, y, a, b)
        union += any_nonzero
    nper = prod(factorial(k - 1) for k in shape)
    total = factorial(10) // nper
    return total, union, union * nper, counts


if __name__ == "__main__":
    import sys

    verify_planar()
    table = sys.argv[1] if len(sys.argv) > 1 else "mixed_graph_probe.tsv"
    with open(table) as stream:
        axes = {int(row.split()[2]) for row in stream if row.strip()}
    configs = [(p, *directions(v)) for p in PLANAR for v in sorted(axes)]
    print("momentum configurations", len(configs))
    for shape in SHAPES:
        total, union, ordered, _ = check(shape, configs)
        print(
            shape,
            "path partitions",
            total,
            "surviving union",
            union,
            "ordered",
            ordered,
        )
