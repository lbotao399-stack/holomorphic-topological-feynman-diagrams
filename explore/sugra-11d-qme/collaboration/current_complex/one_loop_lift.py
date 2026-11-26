"""The one-loop two-one-form lifting map on conserved divergence-free wave symbols."""

from itertools import permutations, product
from math import prod
from random import Random

import numpy as np

from one_loop_vector_counterterms import PRIME, cycles, det, pivots, waves


PERMS = np.asarray([x + g for x in permutations(range(4)) for g in permutations((4, 5))])
SIGNS = np.asarray([(-1) ** (4 - len(cycles(p[:4]))) for p in PERMS])
CHART = [tuple(j for j in range(6) if j != i) for i in range(6)]
MONOMIALS = np.array(list(product(range(4), repeat=6))).T


def source_dimension():
    trace = 0
    for p, sign in zip(PERMS, SIGNS):
        sizes = cycles(p)
        trace += sign * prod(sum(k for k in sizes if r % k == 0) - 2 for r in sizes)
    assert trace == 48 * 83


def solve_frame(p, eta):
    # The six conserved momenta give homogeneous coordinates modulo (1,...,1).
    # Replace one row of the momentum determinant for integral coordinates.
    base = p[:5]
    coords = np.zeros((2, 6), dtype=np.int64)
    for i in range(5):
        for j in range(2):
            repl = base.copy()
            repl[i] = eta[j]
            coords[j, i] = det(repl)
    return coords


def two_forms(p, u, monomials):
    coordinates = solve_frame(p, u[4:])
    pairs = u[:4] @ p.T % PRIME
    # A chart change under a label permutation affects a determinant coordinate
    # by its sign. There are two of them, so their common sign squares to one.
    output = np.ones((48, monomials.shape[1]), dtype=np.int64)
    for i in range(6):
        js = PERMS[:, CHART[i]]
        if i < 4:
            values = pairs[PERMS[:, i, None], js[:, :4]]
        else:
            row = PERMS[:, i, None] - 4
            values = coordinates[row, js[:, :4]] - coordinates[row, js[:, 4:]]
        output = output * values[:, monomials[i]] % PRIME
    return SIGNS @ output % PRIME


def primitive(p, u, monomials):
    perms = np.asarray([x + g for x in permutations(range(3)) for g in permutations((3, 4))])
    signed = np.asarray([(-1) ** (3 - len(cycles(x[:3]))) for x in perms])
    volume = prod(det(np.vstack((p[:4], eta))) for eta in u[3:]) % PRIME
    pairings = u[:3] @ p.T % PRIME
    values = np.ones((12, monomials.shape[1]), dtype=np.int64)
    for i in range(3):
        targets = [j for j in range(5) if j != i][:3]
        row = pairings[perms[:, i, None], perms[:, targets]]
        values = values * row[:, monomials[i]] % PRIME
    return volume * (signed @ values) % PRIME


def outgoing(p, u, monomials, cochain=two_forms):
    out = np.zeros(monomials.shape[1], dtype=np.int64)
    n = len(p) - 2
    for i in range(n):
        for j in range(i + 1, n):
            rest = [k for k in range(n + 2) if k not in (i, j)]
            q = np.vstack((p[i] + p[j], p[rest])) % PRIME
            bracket = (u[i] @ p[j]) * u[j] - (u[j] @ p[i]) * u[i]
            v = np.vstack((bracket, u[rest])) % PRIME
            out += (-1) ** (i + j) * cochain(q, v, monomials)
        for j in (n, n + 1):
            rest = [k for k in range(n) if k != i]
            other = 2 * n + 1 - j
            q = np.vstack((p[rest], p[i] + p[j], p[other])) % PRIME
            # Lie derivative of a one-form, modulo exact one-forms.
            action = (u[i] @ p[j]) * u[j] + (u[i] @ u[j]) * p[i]
            v = np.vstack((u[rest], action, u[other])) % PRIME
            out += (-1) ** (i + 1) * cochain(q, v, monomials)
    return out % PRIME


def scalar_source(p, u, trace_maps):
    four = np.array((p[5], u[5], p[6], u[6]))
    internal = np.array([(-1) ** i * det(np.delete(four, i, axis=1)) for i in range(5)])
    q = np.vstack((p[:5], p[5] + p[6])) % PRIME
    v = np.vstack((u[:5], internal))
    pairings = v @ q.T % PRIME
    out = []
    for maps in trace_maps:
        values = np.ones(len(maps), dtype=np.int64)
        for i in range(6):
            values = values * pairings[i, maps[:, i]] % PRIME
        out.append(int(values.sum()) * det(q[:5]) % PRIME)
    return out


def inputs(n, rng):
    p, u = waves(n, rng)
    # Last two variables are independent one-forms, with no transversality equation.
    u[-2:] = [[rng.randrange(PRIME) for _ in range(5)] for _ in range(2)]
    return p, u


def verify():
    source_dimension()
    rng = Random(21110305)
    trace_types = ((2, 2, 2), (3, 3), (2, 4))
    traces = [np.array([p for p in permutations(range(6)) if cycles(p) == t]) for t in trace_types]
    samples = [two_forms(*inputs(6, rng), MONOMIALS) for _ in range(89)]
    columns = pivots(samples)
    assert len(columns) == 83
    monomials = MONOMIALS[:, columns]
    rows = []
    for _ in range(92):
        p, u = inputs(7, rng)
        rows.append(np.concatenate((outgoing(p, u, monomials), scalar_source(p, u, traces))))
    basis = pivots(rows)
    lower = np.array(list(product(range(3), repeat=3))).T
    boundaries = [outgoing(*inputs(6, rng), lower, primitive) for _ in range(9)]
    exact = len(pivots(boundaries))
    assert (exact, len([i for i in basis if i < 83]), basis[-3:]) == (3, 80, [83, 84, 85])
    def coboundary(p, u, choices):
        return outgoing(p, u, choices, primitive)

    assert not np.any(outgoing(*inputs(7, rng), lower, coboundary))
    p, u = inputs(6, rng)
    original = two_forms(p, u, monomials)
    u[4:] += np.array([[13], [21]]) * p[4:]
    assert np.array_equal(two_forms(p, u, monomials), original)
    print(PRIME, "two-form cochains, incoming rank, outgoing rank, rank with scalar sources",
          len(columns), exact, len(basis) - 3, len(basis))
    return rows, columns


if __name__ == "__main__":
    verify()
