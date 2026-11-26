"""Integrate the three-planar, six-transverse class at one balanced wave geometry."""

import argparse
from functools import cache
from fractions import Fraction
from pathlib import Path
from math import gcd, lcm
import signal
import sys

from complete_nodes import integrate_ledger, partitions
from ordered_detector import geometry
from ordered_shapes import determinant
from six_vector_nodes import BASE, X, Y
from theta_denominators import SHAPES, component_denominator

sys.path.insert(0, str(BASE.parent))
from mixed_six_conserving import WEIGHTS

DIRECTIONS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (1, 1, 0), (1, 0, 1))
SIGNS = (1, 1, -1, 1, -1, -1)


@cache
def fixture():
    transverse = [tuple(s * x for x in p) for s, p in zip(SIGNS, DIRECTIONS)]
    marked = tuple(-sum(c) for c in zip(*transverse))
    momentum = [X, Y, (0, 0, *marked), *((0, 0, *p) for p in transverse)]
    momentum.append(tuple(-sum(c) for c in zip(*momentum)))
    assert momentum[9] == (-1, -1, 0, 0, 0)
    columns = [X, Y, X, *((0, 0, *p) for p in transverse), (0, -1, 0, 0, 0)]
    samples = []
    for mask, weight in enumerate(WEIGHTS):
        if weight:
            polar = [Y, (-1, 0, 0, 0, 0), Y]
            for j, p in enumerate(DIRECTIONS):
                pivot = next(k for k, x in enumerate(p) if x)
                axis = (pivot + 1 + ((mask >> j) & 1)) % 3
                u = [0] * 3
                u[axis], u[pivot] = p[pivot], -p[axis]
                polar.append((0, 0, *u))
            polar.append(X)
            samples.append((weight, polar))
    return (momentum, samples[0][1]), columns, samples


def records(shape):
    _, columns, _ = fixture()
    for order in partitions(shape):
        det = determinant(columns, shape, order)
        if det:
            yield order, det


def payload(n, shape, order, det, prime):
    data, _, samples = fixture()
    return geometry(
        n,
        shape,
        order,
        data,
        -60480 * det * (-1) ** (9 - order.index(9)),
        samples,
        prime,
    )


def run(shape, raw, solver, directory, prime=1009):
    name = "".join(map(str, shape))
    jobs = list(records(shape))
    if solver is None:
        print(name, "supported partitions", len(jobs))
        return
    directory.mkdir(parents=True, exist_ok=True)
    total = integrate_ledger(
        solver,
        directory / (name + ".tsv"),
        f"# conserving six-vector cycle {prime} {name}: order, residue",
        jobs,
        lambda n, order, det: payload(n, shape, order, det, prime),
        prime,
    )
    result = total * pow(2 * raw, -1, prime) % prime
    print(name, "ordered sum", total, "coefficient divided by pi", result, flush=True)
    return result


@cache
def partition_keys(shape):
    return {"".join(map(str, order)) for order, _ in records(shape)}


def ledger_coefficient(directory, shape, raw, prime):
    name = "".join(map(str, shape))
    header, *entries = (directory / (name + ".tsv")).read_text().splitlines()
    assert header == f"# conserving six-vector cycle {prime} {name}: order, residue"
    rows = {}
    for entry in entries:
        key, value = entry.split()
        assert key not in rows and 0 <= int(value) < prime
        rows[key] = int(value)
    assert rows.keys() == partition_keys(shape)
    return sum(rows.values()) * pow(2 * raw, -1, prime) % prime


def verify_sum(directory, prime=1009):
    coefficients = []
    denominators = []
    for shape, raw, _ in SHAPES:
        coefficients.append(ledger_coefficient(directory, shape, raw, prime))
        physical = component_denominator(shape, raw, 1)[1]
        denominators.append(physical // gcd(60480, physical))
    denominator = lcm(*denominators)
    assert gcd(denominator, prime) == 1
    if prime == 1009:
        assert coefficients == [822, 669, 704, 50, 781]
    print("five rational coefficients modulo", prime, coefficients)
    print("total residue", sum(coefficients) % prime, "common denominator", denominator)
    return sum(coefficients) % prime


def verify_exact(directory):
    from exact_projection import primes
    from six_vector_crt import bound, small_residue
    from six_vector_nodes import records as original_records
    from mixed_six_conserving import CLASS_RATIO

    coefficients = []
    for shape, raw, symmetry in SHAPES:
        name = "".join(map(str, shape))
        denominator, limit = bound(shape, raw, symmetry, conserving=True)
        fields = primes()
        residues = {1009: ledger_coefficient(directory, shape, raw, 1009)}
        header, *lines = (directory / (name + "_crt.tsv")).read_text().splitlines()
        assert header.startswith(f"# six-vector {name} ")
        for line in lines:
            prime, value = map(int, line.split())
            assert prime == next(fields)
            residues[prime] = ledger_coefficient(
                directory / str(prime), shape, raw, prime
            )
            assert residues[prime] == value
        numerator, modulus = 0, 1
        for prime, value in residues.items():
            assert gcd(prime, denominator) == 1
            difference = (denominator * value - numerator) % prime
            numerator += modulus * (difference * pow(modulus, -1, prime) % prime)
            modulus *= prime
        if modulus <= 2 * limit:
            raise ValueError(f"{name}: prime product is insufficient for the signed numerator")
        if 2 * numerator > modulus:
            numerator -= modulus
        assert abs(numerator) <= limit
        coefficients.append(Fraction(numerator, denominator))
        print(name, coefficients[-1])
    total = sum(coefficients)
    original = []
    for shape, raw, _ in SHAPES:
        name = "".join(map(str, shape))
        keys = {"".join(map(str, order)) for order, _ in original_records(shape)}
        original.append(small_residue(BASE / "six_vector_shapes", name, keys, raw, False))
    assert original == [472, 106, 290, 904, 865]
    scaled = CLASS_RATIO * total
    assert scaled.numerator * pow(scaled.denominator, -1, 1009) % 1009 == sum(original) % 1009
    print("five-shape sum divided by pi", total)
    print("unrooted coefficient in units of g/pi^10", -total / 4096)
    return total


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--solver", type=Path)
    p.add_argument(
        "--shapes", nargs="+", default=["".join(map(str, s)) for s, _, _ in SHAPES]
    )
    p.add_argument("--directory", type=Path, default=BASE / "six_conserving_shapes")
    p.add_argument("--verify", action="store_true")
    p.add_argument("--verify-exact", action="store_true")
    a = p.parse_args()
    if a.verify_exact:
        verify_exact(a.directory)
        sys.exit(0)
    if a.verify:
        verify_sum(a.directory)
        sys.exit(0)
    for shape, raw, _ in SHAPES:
        if "".join(map(str, shape)) in a.shapes:
            run(shape, raw, a.solver, a.directory)
