"""Recover a complete six-vector shape amplitude with certified integer bounds."""

from argparse import ArgumentParser, Namespace
from fractions import Fraction
from math import ceil, factorial, gcd, prod
from pathlib import Path
from tempfile import TemporaryDirectory
import fcntl
import os
import signal
import subprocess

from exact_projection import primes
from six_vector_nodes import BASE, SHAPES, fields, fixture, records, run
from theta_denominators import component_denominator
from theta_spectral_bound import wick_bound


def bound(shape, raw, symmetry, conserving=False):
    # The automorphism group acts freely on the labelled placements. The
    # reversal factor in the native result cancels the final division by two.
    _, den = component_denominator(shape, raw, 1)
    den //= gcd(den, 60480)
    if conserving:
        from six_conserving_nodes import (
            fixture as balanced_fixture,
            records as balanced_records,
        )

        data, _, samples = balanced_fixture()
        geometries = [(data[0], samples)]
        jobs = balanced_records(shape)
    else:
        geometries = [(fields(0, i)[0][0], rows) for i, rows in fixture().items()]
        jobs = records(shape)
    det_sum = sum(abs(d) for _, d in jobs)
    order_count = prod(factorial(n - 1) for n in shape)
    total = sum(wick_bound(mom, rows) for mom, rows in geometries)
    upper = Fraction(60480 * den * det_sum * order_count, 2 * raw) * total
    return den, ceil(upper)


def compile_general(solver, p):
    subprocess.run(
        [
            "c++",
            "-O3",
            "-std=c++17",
            f"-DTWIST_RESIDUE_PRIME={p}",
            "-I/opt/homebrew/opt/gmp/include",
            "-L/opt/homebrew/opt/gmp/lib",
            str(BASE / "general_partition.cpp"),
            "-lgmpxx",
            "-lgmp",
            "-framework",
            "Accelerate",
            "-o",
            str(solver),
        ],
        check=True,
    )
    if int(subprocess.check_output([solver, "--prime"], text=True)) != p:
        raise RuntimeError("compiled modulus")


def small_residue(base, name, keys, raw, conserving):
    path = base / f"{name}.tsv"
    if not path.exists():
        return None
    with path.open() as control:
        small_rows = {}
        width = 1 if conserving else 2
        for line in control:
            if line.startswith("#"):
                continue
            k, *values = line.split()
            assert k in keys and k not in small_rows and len(values) == width
            small_rows[k] = tuple(map(int, values))
            assert all(0 <= value < 1009 for value in small_rows[k])
    if len(small_rows) < len(keys):
        return None
    return sum(map(sum, small_rows.values())) * pow(2 * raw, -1, 1009) % 1009


def solve(shape, raw, sym, a):
    name = "".join(map(str, shape))
    den, limit = bound(shape, raw, sym, a.conserving)
    print(
        name, "denominator", den, "numerator bit bound", limit.bit_length(), flush=True
    )
    if a.bound_only:
        return
    base = BASE / ("six_conserving_shapes" if a.conserving else "six_vector_shapes")
    base.mkdir(exist_ok=True)
    if a.conserving:
        from six_conserving_nodes import (
            records as balanced_records,
            run as balanced_run,
        )

        expected_jobs = balanced_records(shape)
    else:
        expected_jobs = records(shape)
    keys = {"".join(map(str, o)) for o, _ in expected_jobs}
    with (
        (base / (name + "_crt.tsv")).open("a+") as ledger,
        TemporaryDirectory(prefix="twist-six-general-crt-") as temp,
    ):
        fcntl.flock(ledger, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger.seek(0)
        saved = {}
        for line in ledger:
            if line.startswith("#"):
                continue
            p, x = map(int, line.split())
            if p in saved:
                raise ValueError("repeated prime")
            saved[p] = x
        if not saved and ledger.tell() == 0:
            print(
                "# six-vector",
                name,
                "physical coefficient divided by pi: prime, residue",
                file=ledger,
                flush=True,
            )
        number = 0
        modulus = 1
        for p in primes():
            if p not in saved:
                solver = Path(temp) / "theta"
                compile_general(solver, p)
                args = Namespace(solver=solver, limit=None, directory=base / str(p))
                if a.conserving:
                    saved[p] = balanced_run(shape, raw, solver, args.directory, p)
                else:
                    saved[p] = run(args, shape, raw, p)
                assert isinstance(saved[p], int)
                print(p, saved[p], sep="\t", file=ledger, flush=True)
                os.fsync(ledger.fileno())
            number += modulus * ((den * saved[p] - number) * pow(modulus, -1, p) % p)
            modulus *= p
            if modulus <= 2 * limit < 1009 * modulus:
                small = small_residue(base, name, keys, raw, a.conserving)
                if small is not None:
                    number += modulus * (
                        (den * small - number) * pow(modulus, -1, 1009) % 1009
                    )
                    modulus *= 1009
            if modulus > 2 * limit:
                if number > modulus // 2:
                    number -= modulus
                small = small_residue(base, name, keys, raw, a.conserving)
                assert abs(number) <= limit
                assert small is None or (number - den * small) % 1009 == 0
                assert all(
                    (number - den * val) % field == 0 for field, val in saved.items()
                )
                print(
                    name,
                    "exact coefficient divided by pi",
                    Fraction(number, den),
                    flush=True,
                )
                return Fraction(number, den)


def calibrated_table(options):
    from six_conserving_nodes import verify_exact
    from mixed_six_conserving import CLASS_RATIO

    total = CLASS_RATIO * verify_exact(BASE / "six_conserving_shapes")
    coefficients = [solve(*row, options) for row in SHAPES[:-1]]
    shape, raw, symmetry = SHAPES[-1]
    name = "".join(map(str, shape))
    coefficient = total - sum(coefficients)
    denominator, limit = bound(shape, raw, symmetry)
    integer = coefficient * denominator
    assert integer.denominator == 1 and abs(integer) <= limit
    keys = {"".join(map(str, order)) for order, _ in records(shape)}
    control = small_residue(BASE / "six_vector_shapes", name, keys, raw, False)
    assert coefficient.numerator * pow(coefficient.denominator, -1, 1009) % 1009 == control
    print(name, "by exact cohomology calibration and subtraction", coefficient)
    print("two-geometry five-shape sum divided by pi", total)


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    ap = ArgumentParser(description=__doc__)
    ap.add_argument("shapes", nargs="*")
    ap.add_argument("--bound-only", action="store_true")
    ap.add_argument(
        "--calibrated-table",
        action="store_true",
        help="recover the original four shapes and infer its last shape from the exact calibrated sum",
    )
    ap.add_argument(
        "--conserving",
        action="store_true",
        help="small integer cycle at one balanced geometry",
    )
    a = ap.parse_args()
    if a.calibrated_table:
        if a.shapes or a.bound_only or a.conserving:
            ap.error("use --calibrated-table by itself")
        calibrated_table(a)
        raise SystemExit
    if not a.shapes or set(a.shapes) - {"166", "256", "346", "355", "445"}:
        ap.error("supply one or more of 166, 256, 346, 355, 445")
    if "166" in a.shapes and not a.conserving:
        ap.error("use six_vector_166_exact.py for the two-geometry 166 ledger")
    for shape, raw, sym in SHAPES:
        if "".join(map(str, shape)) in a.shapes:
            solve(shape, raw, sym, a)
