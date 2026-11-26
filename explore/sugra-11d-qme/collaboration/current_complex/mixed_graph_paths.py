"""Exact subset recursion for the ordered theta path with arbitrary current moments."""

from itertools import permutations
from math import comb

D = 12
N = 5


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def matrix(f):
    return tuple(tuple(f(i, j) for j in range(N)) for i in range(N))


ZERO = matrix(lambda i, j: 0)
IDENTITY = matrix(lambda i, j: int(i == j))


def product(a, b):
    cols = tuple(zip(*b))
    return tuple(tuple(dot(row, col) for col in cols) for row in a)


def accumulate(a, b, c=1):
    return matrix(lambda i, j: a[i][j] + c * b[i][j])


def epsilon(indices):
    if len(set(indices)) < N:
        return 0
    return (-1) ** sum(x > y for i, x in enumerate(indices) for y in indices[i + 1 :])


def transfer(p, u, enter, forward, gamma):
    if gamma == "source":
        return matrix(lambda i, j: u[i] * p[j] - p[i] * u[j])
    if gamma == "sink":
        return matrix(
            lambda i, j: sum(
                epsilon((i, j, k, l, s)) * u[k] * enter[l] * p[s]
                for k in range(N)
                for l in range(N)
                for s in range(N)
            )
        )
    current = enter if forward else add(enter, p)
    trace = dot(current, u)
    if forward:
        return matrix(lambda i, j: u[i] * current[j] - (trace if i == j else 0))
    return matrix(lambda i, j: -u[j] * current[i] + (trace if i == j else 0))


def subsets(momenta, polar, initial, direction, kind, gamma_position=None):
    # Label 0 carries the one-form. Other labels carry divergence-free vectors.
    n = len(momenta)
    endpoint = 1 << n
    prefix = [(0,) * N] * endpoint
    path = {0: [IDENTITY] + [ZERO] * D}
    for s in range(1, endpoint):
        bit = s & -s
        prefix[s] = add(prefix[s ^ bit], momenta[bit.bit_length() - 1])
        row = [ZERO] * (D + 1)
        for j, (p, u) in enumerate(zip(momenta, polar)):
            if not s & (1 << j):
                continue
            r = s ^ (1 << j)
            if gamma_position is not None and (
                (j == 0) != (r.bit_count() == gamma_position)
            ):
                continue
            forward = (kind == "sink") != bool(r & 1)
            v = transfer(
                p, u, add(initial, prefix[r]), forward, kind if j == 0 else None
            )
            for d in range(D + 1):
                row[d] = accumulate(row[d], product(path[r][d], v))
        q = dot(direction, prefix[s])
        for d in range(1, D + 1):
            row[d] = accumulate(row[d], row[d - 1], q)
        path[s] = row
    return path[endpoint - 1]


def ordered(momenta, polar, initial, direction, kind, gamma_position=None):
    # Integrate the exponential generating series directly over ordered tail times.
    # Scaling its coefficient of t^(k+d) by (k+d)! keeps every entry integral.
    n = len(momenta)
    results = [ZERO] * (D + 1)
    for order in permutations(range(n)):
        if gamma_position is not None and order[gamma_position] != 0:
            continue
        coefficients = [1] + [0] * D
        for k, j in enumerate(reversed(order)):
            a = dot(direction, momenta[j])
            coefficients = [
                sum(comb(k + d, r) * a**r * coefficients[d - r] for r in range(d + 1))
                for d in range(D + 1)
            ]
        mat = IDENTITY
        p = initial
        forward = kind == "sink"
        for j in order:
            mat = product(
                mat,
                transfer(momenta[j], polar[j], p, forward, kind if j == 0 else None),
            )
            p = add(p, momenta[j])
            if j == 0:
                forward = not forward
        for d, c in enumerate(coefficients):
            results[d] = accumulate(results[d], mat, c)
    return results


if __name__ == "__main__":
    # The first label is the external one-form. The other polarizations are transverse.
    momenta = [
        (-2, 1, -3, 1, -2),
        (1, 0, 0, 0, 0),
        (1, 1, 1, 0, 0),
        (0, 0, 0, 1, 0),
        (0, 0, 0, 1, 1),
    ]
    polar = [
        (0, 0, 1, 0, 0),
        (0, 1, 0, 0, 0),
        (-1, 1, 0, 0, 0),
        (0, 0, 0, 0, 1),
        (0, 0, 0, 1, -1),
    ]
    assert all(not dot(p, u) for p, u in zip(momenta[1:], polar[1:]))
    initial = (1, -2, 3, 4, -1)
    direction = (2, -1, 3, 1, -2)
    for kind in ["source", "sink"]:
        for gamma_position in [None, 2]:
            args = (momenta, polar, initial, direction, kind, gamma_position)
            a = subsets(*args)
            b = ordered(*args)
            assert a == b
            assert any(x for mat in a for row in mat for x in row)
            print(kind, gamma_position, "matched tensor coordinates", len(a) * N * N)
