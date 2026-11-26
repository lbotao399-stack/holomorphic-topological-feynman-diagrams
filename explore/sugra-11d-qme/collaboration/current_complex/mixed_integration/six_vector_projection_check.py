"""Check Taylor extraction against the original two-scale integrated Wick polynomial."""

import argparse
from pathlib import Path
from subprocess import run

from ordered_detector import geometry
from ordered_shapes import determinant
from six_vector_nodes import PRIME, SHAPES, fields, fixture, payload, records


def raw(i, shape, order, choice, a):
    (p, u), _ = fields(0, choice)
    p[2] = (a, *p[2][1:])
    p[9] = tuple(-sum(c) for c in zip(*p[:9]))
    d = determinant(p, shape, order)
    coefficient = -60480 * (-1) ** (9 - order.index(9))
    return geometry(i, shape, order, (p, u), d * coefficient, fixture()[choice])


def verify(solver):
    jobs = []
    for shape, _, _ in SHAPES:
        order, det = next(records(shape))
        for choice in (0, 1):
            i = len(jobs)
            jobs.append(payload(i, shape, order, det, choice))
            jobs.append(raw(i + 1, shape, order, choice, 1))
            jobs.append(raw(i + 2, shape, order, choice, 2))
    out = run([solver], input="".join(jobs), capture_output=True, text=True, check=True)
    lines = [tuple(map(int, row.split())) for row in out.stdout.splitlines()]
    assert [i for i, x in lines] == list(range(len(jobs)))
    for i in range(0, len(lines), 3):
        v, one, two = [x for _, x in lines[i : i + 3]]
        assert (4 * one - two - 2 * v) % PRIME == 0
    print(
        "integrated Taylor identities",
        len(jobs) // 3,
        "original Wick evaluations",
        len(jobs) * 2 // 3,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ordered_solver", type=Path)
    verify(p.parse_args().ordered_solver)
