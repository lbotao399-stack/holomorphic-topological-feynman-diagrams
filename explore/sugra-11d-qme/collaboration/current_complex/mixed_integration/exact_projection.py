"""Recover one projected rational partition by CRT with an absolute numerator bound."""

from fractions import Fraction
from math import factorial, isqrt, lcm, prod
from pathlib import Path
import argparse
import subprocess
import sys
from tempfile import TemporaryDirectory
from single_projection import cases, residual_det
from detector_factors import SINGLE_WEIGHTS
from detector_node import single_fields

ROOT = Path(__file__).resolve().parent


def bound(order):
    # Separate resistance degree seven, joint degree eight and conductance degree seven.
    path = lcm(
        *(
            factorial(5 + i) * factorial(5 + j) // factorial(5) ** 2
            for i in range(8)
            for j in range(8)
            if i + j <= 8
        )
    )
    angular = prod(range(13, 27, 2))
    abs_weights = sum(
        abs(int(row.split()[0])) for row in SINGLE_WEIGHTS.read_text().splitlines()
    )
    counter = 0
    total = 0
    for i, br, sc, mask, weight in cases():
        (mom, u0), columns = single_fields(br, 0, 0, 0, sc, mask)
        det = residual_det(order, br, columns)
        if not det:
            continue
        counter += 1
        u1 = single_fields(br, 31, 0, 0, sc, mask)[0][1]
        # Each edge coordinate is bounded by the total variation of that momentum component.
        limits = [sum(abs(p[j]) for p in mom) for j in range(5)]
        maximum = max(limits)
        ordinary = 1
        for p, us in zip(mom[:9], zip(u0, u1)):
            norm = 0
            for u in us:
                trace = sum(abs(x) * y for x, y in zip(u, limits))
                rank = max(
                    max(map(abs, u)) * sum(map(abs, p)),
                    max(map(abs, p)) * sum(map(abs, u)),
                )
                norm = max(norm, trace + rank)
            ordinary *= norm
        pg = mom[9]
        source = sum(map(abs, pg)) + max(map(abs, pg))
        sink = 30 * maximum * max(map(abs, pg))
        endpoints = 2 * (2 * maximum**3) * source + 2 * maximum**2 * sink
        total += abs(weight * det) * ordinary * endpoints
    denominator = 1536 * path * angular
    # Five independent fiber indices remain after multiplying the two paths.
    numerator_bound = path * angular * factorial(5) ** 2 * 5**5 * abs_weights * total
    return denominator, numerator_bound, counter


def primes():
    value = 6_000_001
    while True:
        if all(value % j for j in range(3, isqrt(value) + 1, 2)):
            yield value
        value += 2


def compile_solver(out, p):
    subprocess.run(
        [
            "c++",
            "-O3",
            "-std=c++17",
            f"-DTWIST_RESIDUE_PRIME={p}",
            "-I/opt/homebrew/opt/gmp/include",
            "-L/opt/homebrew/opt/gmp/lib",
            str(ROOT / "native_partition.cpp"),
            "-lgmpxx",
            "-lgmp",
            "-framework",
            "Accelerate",
            "-o",
            str(out),
        ],
        check=True,
    )


def run(a):
    order = tuple(map(int, a.order))
    den, maximum, count = bound(order)
    needed = 0
    mod = 1
    for p in primes():
        needed += 1
        mod *= p
        if mod > 2 * maximum:
            break
    print(
        "common denominator",
        den,
        "integer numerator bound",
        maximum,
        "primes",
        needed,
        "cases",
        count,
        flush=True,
    )
    if a.bound_only:
        return
    checkpoint = ROOT / f"single_exact_{a.order}.tsv"
    data = {}
    if checkpoint.exists():
        for line in checkpoint.read_text().splitlines():
            p, v = map(int, line.split())
            data[p] = v
    modulus = 1
    number = 0
    with TemporaryDirectory(prefix="twist166-crt-") as temp:
        solver = Path(temp) / "partition"
        for p in primes():
            if p not in data:
                compile_solver(solver, p)
                pending = ROOT / f"single_exact_{a.order}_{p}.partial.tsv"
                cmd = [
                    sys.executable,
                    str(ROOT / "single_projection.py"),
                    str(solver),
                    "--order",
                    a.order,
                    "--checkpoint",
                    str(pending),
                ]
                result = subprocess.run(cmd, check=True, text=True, capture_output=True)
                words = result.stdout.split()
                if words[-2:] != ["mod", str(p)] or int(words[4]) != count:
                    raise ValueError(result.stdout)
                data[p] = int(words[-3])
                with checkpoint.open("a") as f:
                    print(p, data[p], file=f)
                pending.unlink()
                print("residue", p, data[p], flush=True)
            number += modulus * ((data[p] * den - number) * pow(modulus, -1, p) % p)
            modulus *= p
            if modulus > 2 * maximum:
                break
    if 2 * number > modulus:
        number -= modulus
    if abs(number) > maximum:
        raise ArithmeticError("reconstruction exceeds proved bound")
    answer = Fraction(number, den)
    if a.verify1009 is not None and (number - a.verify1009 * den) % 1009:
        raise ArithmeticError("independent small-prime check failed")
    print(
        "exact normalized projected partition",
        answer,
        "CRT product",
        modulus,
        flush=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--order", default="0234516789")
    p.add_argument("--bound-only", action="store_true")
    p.add_argument("--verify1009", type=int)
    run(p.parse_args())
