"""Independent integrated test of a linear insertion across the five theta shapes."""

import argparse
from fractions import Fraction
from itertools import permutations
from random import Random
import subprocess

from detector_node import single_fields
from ordered_detector import geometry
from ordered_shapes import determinant, expectation, numerator
from theta_denominators import SHAPES

DENOMINATORS = {s: k for s, k, _ in SHAPES}


def fixtures():
    rng = Random(867)
    profiles = ((0, 4, 5), (1, 3, 5), (1, 4, 4), (2, 2, 5), (2, 3, 4), (3, 3, 3))
    for i in range(120):
        bits, eta, br = rng.randrange(32), rng.randrange(3), rng.randrange(2)
        (momentum, polar), _ = single_fields(br, bits, 0, eta, (1, 2, 3, -1), 63)
        lin = rng.randrange(9)
        momentum = list(momentum)
        q = momentum[lin]
        momentum[lin] = (0,) * 5
        momentum[9] = tuple(a + b for a, b in zip(momentum[9], q))
        columns = momentum.copy()
        columns[lin] = q
        seq = [j for j in range(10) if j != lin]
        rng.shuffle(seq)
        counts = rng.choice(profiles)
        cut = counts[0] + counts[1]
        yield (
            i,
            lin,
            (seq[: counts[0]], seq[counts[0] : cut], seq[cut:]),
            (momentum, polar),
            columns,
        )


def contraction(groups, data, columns, solver):
    shape = tuple(len(g) + 1 for g in groups)
    order = tuple(i for path in groups for i in path)
    det = determinant(columns, shape, order) * (-1) ** (9 - order.index(9))
    if not det:
        return 0
    den = DENOMINATORS[tuple(sorted(shape))]
    if solver is None:
        return Fraction(det, den) * expectation(
            numerator(data, shape, order, True), shape, order, data[0]
        )
    solver.stdin.write(geometry(0, shape, order, data, det, [(1, data[1])]))
    solver.stdin.flush()
    row = solver.stdout.readline().split()
    if len(row) != 2 or row[0] != "0":
        raise RuntimeError(row)
    return int(row[1]) * pow(den, -1, 1009) % 1009


def check(solver):
    checked = nonzero = 0
    shapes = set()
    for i, lin, groups, data, cols in fixtures():
        if solver is None and i not in (1, 9, 12):
            continue
        sums = []
        for h, g in enumerate(groups):
            subtotal = 0
            for j in range(len(g) + 1):
                gg = [v.copy() for v in groups]
                gg[h].insert(j, lin)
                val = contraction(gg, data, cols, solver)
                if val:
                    shapes.add(tuple(sorted(len(path) + 1 for path in gg)))
                if solver is not None and val and checked < 60:
                    for perm in permutations(range(3)):
                        assert (
                            contraction([gg[t] for t in perm], data, cols, solver)
                            == val
                        )
                        checked += 1
                subtotal += val
            sums.append(subtotal if solver is None else subtotal % 1009)
        assert sum(sums) == 0 if solver is None else sum(sums) % 1009 == 0
        nonzero += int(any(sums))
        if solver is None:
            print(
                "rational fixture",
                i,
                "path subtotals",
                *[str(v) for v in sums],
                flush=True,
            )
    if solver is not None:
        assert nonzero == 70 and checked == 60 and shapes == set(DENOMINATORS)
        print(
            "120 insertions checked, 70 nontrivial path cancellations, 60 permutation checks"
        )
    else:
        assert nonzero == 3


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("solver", nargs="?")
    args = parser.parse_args()
    if args.solver:
        with subprocess.Popen(
            [args.solver], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True
        ) as solver:
            try:
                check(solver)
            finally:
                solver.stdin.close()
                if solver.wait():
                    raise RuntimeError("solver exit")
    else:
        check(None)
