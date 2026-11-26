"""Reconstruct the single planar cycle on a complete theta shape."""

import argparse
import fcntl
import os
import signal
import subprocess
from fractions import Fraction
from functools import cache
from math import ceil, factorial, gcd, lcm, prod
from pathlib import Path
from tempfile import TemporaryDirectory

from exact_projection import compile_solver, primes
from four_small_bound import RATIO as small_cycle_ratio, physical_bound
from single_planar_nodes import ROOT, PRIME, fields, records, run
from single_projection import detector_rows
from six_vector_166_exact import path_parts
from six_vector_crt import compile_general
from theta_denominators import SHAPES, component_denominator
from theta_spectral_bound import wick_bound


@cache
def polars():
    return [(w, fields(bits, eta)[1]) for w, bits, eta in detector_rows()]


@cache
def maximum_wick():
    return wick_bound(fields(0, 0)[0], polars())


def bound(shape, raw, symmetry):
    # Every complete external-label orbit cancels the graph symmetry divisor.
    _, den = component_denominator(shape, raw, 1)
    den //= gcd(den, 665280)
    orders = prod(factorial(n - 1) for n in shape)
    determinants = sum(abs(d) for _, d in records(shape, None))
    upper = Fraction(665280 * den * determinants * orders, 2 * raw) * maximum_wick()
    return den, ceil(upper)


@cache
def special_parts(order):
    return path_parts(order, fields(0, 0), polars())


def shape166(solver, prime):
    inputs = []
    shape = (1, 6, 6)
    for i, (order, det) in enumerate(records(shape, prime)):
        paths, weights, ends = special_parts(order)
        rows = [
            (i, 665280 * det % prime, 120),
            *paths,
            *[[x % prime for x in row] for row in weights],
            *ends,
        ]
        inputs.append("".join(" ".join(map(str, row)) + "\n" for row in rows))
    proc = subprocess.run(
        [solver],
        input="".join(inputs),
        text=True,
        capture_output=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
        check=True,
    )
    values = [tuple(map(int, row.split())) for row in proc.stdout.splitlines()]
    assert [i for i, _ in values] == list(range(len(inputs)))
    return sum(v for _, v in values) * pow(2 * 158400, -1, prime) % prime


def small_check(shape, raw):
    name = "".join(map(str, shape))
    path = ROOT / "single_planar_1009" / (name + ".tsv")
    lines = path.read_text().splitlines()
    assert (
        lines[0] == f"# 4+5 single planar four-cycle x 665280 reflection v1 {name} 1009"
    )
    values = {key: int(v) for key, v in (row.split() for row in lines[1:])}
    keys = {"".join(map(str, o)) for o, _ in records(shape, PRIME)}
    assert keys == set(values) and len(lines) == len(values) + 1
    return sum(values.values()) * pow(2 * raw, -1, PRIME) % PRIME


def solve(shape, raw, symmetry, args):
    name = "".join(map(str, shape))
    den, upper = bound(shape, raw, symmetry)
    print(name, "denominator", den, "integer numerator bound", upper, flush=True)
    if args.bound_only:
        return
    base = ROOT / "single_planar_crt"
    base.mkdir(exist_ok=True)
    with (
        (base / (name + ".tsv")).open("a+") as out,
        TemporaryDirectory(prefix="twist-single-ce-") as temp,
    ):
        fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.seek(0)
        header = f"# single planar CE cycle {name} x 665280 coefficient of pi"
        lines = out.read().splitlines()
        if lines and lines[0] != header:
            raise ValueError("incompatible residues")
        if not lines:
            print(header, file=out, flush=True)
        saved = {}
        for row in lines[1:]:
            p, v = map(int, row.split())
            if p in saved or not 0 <= v < p:
                raise ValueError("invalid residue")
            saved[p] = v
        modulus, number = 1, 0
        for p in primes():
            if p not in saved:
                solver = Path(temp) / "theta"
                if shape == (1, 6, 6):
                    compile_solver(solver, p)
                    saved[p] = shape166(solver, p)
                else:
                    compile_general(solver, p)
                    saved[p] = run(
                        argparse.Namespace(
                            solver=solver, prime=p, shape=name, directory=base / str(p)
                        )
                    )
                print(p, saved[p], file=out, flush=True)
                os.fsync(out.fileno())
            number += modulus * ((den * saved[p] - number) * pow(modulus, -1, p) % p)
            modulus *= p
            if modulus > 2 * upper:
                if number > modulus // 2:
                    number -= modulus
                assert (
                    abs(number) <= upper
                    and (number - den * small_check(shape, raw)) % PRIME == 0
                )
                answer = Fraction(number, den)
                print(
                    name,
                    "coefficient of pi",
                    answer,
                    "modulus bits",
                    modulus.bit_length(),
                    flush=True,
                )
                return answer


def total_bound():
    # mixed_probe_jets.py verifies this rational ratio on an integral pair
    # of transverse cycles. Applying both to the same single planar wave
    # gives the same guaranteed physical denominator for both evaluations.
    ratio = Fraction(
        -1029229447699740945723408289,
        114415887442730778653770948594609901093939178844554683860158369901482612,
    )
    rows = [bound(*data) for data in SHAPES]
    den = lcm(*(d for d, _ in rows))
    # Both calibrated cycles have integral denominators dividing den. If the
    # original integer numerator is r*n, its small-cycle value is n*r/RATIO.
    # The denominator of r/RATIO therefore also divides n.
    divisor = ratio.denominator * (ratio.denominator / small_cycle_ratio).denominator
    upper = min(
        sum(limit * (den // d) for d, limit in rows) // divisor,
        int(den * physical_bound() / divisor),
    )
    print(
        "full coefficient denominator divides",
        den,
        "scaled numerator divisible by",
        divisor,
        "remaining integer bound",
        upper,
        flush=True,
    )
    return den, divisor, upper


def total(args):
    den, divisor, upper = total_bound()
    if args.bound_only:
        return
    base = ROOT / "single_planar_crt"
    base.mkdir(exist_ok=True)
    modulus, number = 1, 0
    with (
        (base / "total.tsv").open("a+") as out,
        TemporaryDirectory(prefix="twist-ce-total-") as temp,
    ):
        fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.seek(0)
        header = "# five 4+5 shape coefficients of pi on the single planar cycle after 665280"
        lines = out.read().splitlines()
        if lines and lines[0] != header:
            raise ValueError("incompatible total")
        if not lines:
            print(header, file=out, flush=True)
        saved = {}
        for row in lines[1:]:
            prime, *values = map(int, row.split())
            if (
                len(values) != 5
                or prime in saved
                or any(not 0 <= x < prime for x in values)
            ):
                raise ValueError("invalid five-shape residue")
            saved[prime] = values
        for prime in primes():
            if prime not in saved:
                specialized = base / "166.tsv"
                available = {}
                if specialized.exists():
                    available = dict(
                        tuple(map(int, row.split()))
                        for row in specialized.read_text().splitlines()
                        if row and not row.startswith("#")
                    )
                solver = Path(temp) / "theta"
                if prime not in available:
                    compile_solver(solver, prime)
                    available[prime] = shape166(solver, prime)
                values = [available[prime]]
                compile_general(solver, prime)
                for shape, _, _ in SHAPES[1:]:
                    name = "".join(map(str, shape))
                    values.append(
                        run(
                            argparse.Namespace(
                                solver=solver,
                                prime=prime,
                                shape=name,
                                directory=base / str(prime),
                            )
                        )
                    )
                print(prime, *values, file=out, flush=True)
                os.fsync(out.fileno())
                saved[prime] = values
            residue = (den * sum(saved[prime]) * pow(divisor, -1, prime)) % prime
            number += modulus * ((residue - number) * pow(modulus, -1, prime) % prime)
            modulus *= prime
            if modulus * PRIME > 2 * upper:
                control = sum(small_check(shape, raw) for shape, raw, _ in SHAPES) % PRIME
                if modulus <= 2 * upper:
                    residue = den * control * pow(divisor, -1, PRIME) % PRIME
                    number += modulus * ((residue - number) * pow(modulus, -1, PRIME) % PRIME)
                    modulus *= PRIME
                if number > modulus // 2:
                    number -= modulus
                assert abs(number) <= upper
                assert (number * divisor - den * control) % PRIME == 0
                result = Fraction(number * divisor, den)
                print(
                    "full 4+5 coefficient of pi",
                    result,
                    "modulus bits",
                    modulus.bit_length(),
                    flush=True,
                )
                return result


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "shapes", nargs="+", choices=["166", "256", "346", "355", "445", "total"]
    )
    parser.add_argument("--bound-only", action="store_true")
    args = parser.parse_args()
    if "total" in args.shapes:
        total(args)
    for sh, raw, sym in SHAPES:
        if "".join(map(str, sh)) in args.shapes:
            solve(sh, raw, sym, args)
