"""Certify the integer six-vector detector on the second one-mark graph test."""

import argparse
from collections import defaultdict
from functools import reduce
from math import gcd
from pathlib import Path
from sys import stdin

from mixed_probe_jets import matrix, pivots


def verify(lines, weights, balanced=False, linear=False):
    data = defaultdict(list)
    for line in lines:
        tag, *row = map(int, line.split())
        data[tag].append(row)
    if balanced:
        assert len(data[0]) == 45 and len(data[1]) == 34
        target, *others = data[1]
        columns = pivots(r[:150] for r in others)
        assert len(columns) == 33
        square = matrix([[r[j] for r in others] for j in columns])
        rhs = matrix([[-target[j]] for j in columns])
        num, den = square.solve_den(rhs)
        factors = [int(den), *map(int, num.to_Matrix())]
        common = reduce(gcd, factors)
        factors = [x // common for x in factors]
        total = [sum(v * w for v, w in zip(factors, col)) for col in zip(*data[1])]
        assert not any(total[:150])
        assert len(pivots(data[0])) == 28 and len(pivots(data[0] + [total[150:]])) == 29
        print(
            "balanced integer test",
            len(factors),
            "boundary columns",
            150,
            "augmented rank",
            29,
        )
        return
    assert [len(data[i]) for i in range(3)] == [45, 64, 64]
    outgoing = data[0]
    evaluations = data[1] + data[2]
    assert len(pivots(outgoing)) == 28
    assert len(pivots(row[:150] for row in evaluations)) == 76

    total = [0] * 255
    seen = set()
    for line in weights.read_text().splitlines():
        weight, choice, mask = map(int, line.split())
        assert (choice, mask) not in seen and choice in (0, 1) and 0 <= mask < 64
        seen.add((choice, mask))
        for i, value in enumerate(data[choice + 1][mask]):
            total[i] += weight * value
    assert len(seen) == 77 and not any(total[:150])
    assert len(pivots(outgoing + [total[150:]])) == 29
    if linear:
        (grades,) = data[3]
        assert [grades.count(s) for s in (1, 2, 3)] == [70, 28, 7]
        keep = [i for i, s in enumerate(grades) if s == 1]
        restricted = [[row[i] for i in keep] for row in outgoing]
        detector = [total[150 + i] for i in keep]
        columns = pivots(restricted + [detector])
        row_ids = pivots(zip(*[[row[i] for i in columns] for row in restricted]))
        assert len(columns) == len(row_ids) + 1
        square = matrix(
            [[restricted[j][i] for i in columns] for j in row_ids]
            + [[detector[i] for i in columns]]
        )
        num, den = square.solve_den(matrix([[0] for _ in row_ids] + [[1]]))
        factors = list(map(int, num.to_Matrix()))
        assert all(
            sum(row[i] * v for i, v in zip(columns, factors)) == 0 for row in restricted
        )
        assert sum(detector[i] * v for i, v in zip(columns, factors)) == den and den
        print(
            "closed scalar class supported on first derivative",
            len(keep),
            "closure rank",
            len(row_ids),
        )
    print(
        "integer boundary pairings",
        len(total[:150]),
        "nonzero",
        sum(x != 0 for x in total[:150]),
    )
    print("six-to-seven rank and augmented rank modulo 251", 28, 29)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--weights",
        type=Path,
        default=Path(__file__).with_name("mixed_six_weights.tsv"),
    )
    parser.add_argument("--balanced", action="store_true")
    parser.add_argument("--linear", action="store_true")
    args = parser.parse_args()
    verify(stdin, args.weights, args.balanced, args.linear)
