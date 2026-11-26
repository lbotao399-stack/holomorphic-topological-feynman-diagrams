"""Full symmetry-weighted rational 166 kernel on three or four linear vectors."""

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from math import factorial, lcm, prod
from pathlib import Path
import os
import subprocess
from tempfile import TemporaryDirectory

from detector_node import single_fields
from exact_projection import compile_solver, primes
from projective_shape166 import partition_orbits
from single_projection import input_geometry, residual_det

MASK = 35
ROWS = ((1, 0, 1),)


def assignments(scales):
    (mom, pol), columns = single_fields(0, 0, 0, 1, scales, MASK)
    groups = defaultdict(list)
    for order in sum(partition_orbits().values(), []):
        det = residual_det(order, 0, columns)
        if det:
            groups[tuple(v for v in order[:5] if any(mom[v]))].append((order, det))
    return groups, mom, pol


def evaluate(solver, prime, groups, scales):
    total = 0
    with subprocess.Popen(
        [solver],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
    ) as proc:
        try:
            for members in groups:
                multiplicity = sum(d for _, d in members)
                if not multiplicity:
                    continue
                order, det = members[0]
                proc.stdin.write(
                    input_geometry(0, order, 0, scales, MASK, 120, prime, ROWS)
                )
                proc.stdin.flush()
                i, value = map(int, proc.stdout.readline().split())
                if i:
                    raise ValueError("unexpected response")
                total += multiplicity * value * pow(det, -1, prime)
        finally:
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError("native detector stopped")
    return total % prime


def bound(groups, mom, pol):
    # Bound the unprojected rational normalized sum with the pure columns removed.
    path = lcm(
        *(
            factorial(5 + i) * factorial(5 + j) // factorial(5) ** 2
            for i in range(8)
            for j in range(8)
            if i + j <= 8
        )
    )
    den = path * prod(range(13, 27, 2))
    coordinates = [sum(abs(p[i]) for p in mom) for i in range(5)]
    largest = max(coordinates)
    ordinary = 1
    for p, u in zip(mom[:9], pol[:9]):
        trace = sum(abs(x) * y for x, y in zip(u, coordinates))
        rank = max(
            max(map(abs, u)) * sum(map(abs, p)), max(map(abs, p)) * sum(map(abs, u))
        )
        ordinary *= trace + rank
    pg = mom[9]
    endpoints = 4 * largest**3 * (
        sum(map(abs, pg)) + max(map(abs, pg))
    ) + 60 * largest**3 * max(map(abs, pg))
    total_det = sum(abs(d) for row in groups for _, d in row)
    return den, den * total_det * factorial(5) ** 2 * 5**5 * endpoints * ordinary


def main(four):
    scales = (1, 1, 0 if four else 2, 1)
    groups, mom, pol = assignments(scales)
    rows = list(groups.values())
    den, limit = bound(rows, mom, pol)
    number, modulus = 0, 1
    with TemporaryDirectory(prefix="mixed-three-linear-") as directory:
        solver = Path(directory) / "native"
        for p in primes():
            compile_solver(solver, p)
            res = evaluate(solver, p, rows, scales)
            number += modulus * ((den * res - number) * pow(modulus, -1, p) % p)
            modulus *= p
            if modulus > 2 * limit:
                break
        if 2 * number > modulus:
            number -= modulus
        if abs(number) > limit:
            raise ArithmeticError("numerator bound")
        answer = Fraction(number, den)
        compile_solver(solver, 1009)
        # Every placement is evaluated independently at the held-out prime.
        check = evaluate(solver, 1009, [[pair] for row in rows for pair in row], scales)
        if check != (687 if four else 994) or (number - den * check) % 1009:
            raise ArithmeticError(("ungrouped 1009 comparison", check))
    print(
        "raw partitions",
        sum(map(len, rows)),
        "active groups",
        sum(bool(sum(d for _, d in row)) for row in rows),
    )
    print("normalized sum", answer, "product", modulus, "numerator bound", limit)
    # Each zero scale takes its first Taylor coefficient.
    print("coefficient of pi", answer * prod(x for x in scales[:3] if x) / (2 * 158400))


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--four", action="store_true")
    main(parser.parse_args().four)
