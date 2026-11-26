"""Single planar CE cycle and reflected partition orbits for the mixed 4+5 current."""

import argparse
from functools import cache
from pathlib import Path
import signal
import subprocess
import sys
from fractions import Fraction
from math import factorial

from complete_nodes import ROOT, integrate_ledger, partitions
from detector_node import single_fields
from ordered_detector import geometry
from ordered_shapes import determinant, paths
from single_projection import detector_rows
from theta_denominators import SHAPES

PRIME = 1009
PLANAR = ((1, 0), (0, 1), (1, 1), (-1, -1))


def verify_planar():
    sys.path.insert(0, str(ROOT.parent))
    from block_homotopy import planar_dimensions
    from mixed_current_lift import lift
    from mixed_planar_soft import ordinary_wave_coefficient, scalar_boundaries

    def evaluate(cochain):
        total = Fraction(0)
        for word, coeff in cochain.items():
            if coeff:
                _, h, i = word[-1]
                j = h + 2 - i
                total += (
                    coeff
                    * ordinary_wave_coefficient(word, PLANAR[:3])
                    * Fraction((-1) ** (i + j), factorial(i) * factorial(j))
                )
        return total

    _, scalar, pairing, *_ = lift()
    boundaries = scalar_boundaries()
    assert len(boundaries) == 6 and not any(evaluate(b) for b in boundaries)
    assert evaluate(scalar) == Fraction(32, 9)
    assert pairing == 665280 * evaluate(scalar)
    assert planar_dimensions(6, 0) == {2: 1, 3: 1}
    print(
        "planar boundary values zero, class value 32/9, integer four-cycle factor 665280"
    )


def fields(bits, eta):
    (p, u), _ = single_fields(1, bits, 0, eta, (1, 1, 1, 1), 63)
    for j, (x, y) in enumerate(PLANAR):
        p[j] = (x, y, *p[j][2:])
        u[j] = (-y, x, 0, 0, 0)
    p[9] = tuple(-sum(c) for c in zip(*p[:9]))
    return p, u


@cache
def fixture(prime):
    return [
        (w % prime, fields(bits, eta)[1])
        for w, bits, eta in detector_rows()
        if w % prime
    ]


def reflected(order, shape):
    groups = [
        tuple(sorted(1 - j if j in (0, 1) else j for j in group))
        for group in paths(shape, order)
    ]
    if shape[0] == shape[1] and groups[0] > groups[1]:
        groups[0], groups[1] = groups[1], groups[0]
    if shape[1] == shape[2] and shape[0] != 1 and groups[1] > groups[2]:
        groups[1], groups[2] = groups[2], groups[1]
    return sum(groups, ())


def records(shape, prime, fold=True):
    columns = fields(0, 0)[0]
    for order in partitions(shape):
        other = reflected(order, shape)
        if fold and order > other:
            continue
        det = determinant(columns, shape, order)
        if prime is not None:
            det %= prime
        if det:
            other_det = determinant(columns, shape, other)
            if reflected(other, shape) != order or not (
                other_det if prime is None else other_det % prime
            ):
                raise ArithmeticError("reflection orbit")
            yield order, det * (1 + int(fold and other != order))


def payload(node, shape, order, det, prime):
    signed = 665280 * det * (-1) ** (9 - order.index(9))
    return geometry(node, shape, order, fields(0, 0), signed, fixture(prime), prime)


def run(args):
    prime = args.prime
    shape = tuple(map(int, args.shape))
    jobs = list(records(shape, prime))
    full = sum(1 + int(reflected(order, shape) != order) for order, _ in jobs)
    assert full == sum(1 for _ in records(shape, prime, False))
    print(
        args.shape,
        len(jobs),
        "reflection orbits in",
        full,
        "supported partitions",
        flush=True,
    )
    if not args.solver:
        return
    if (
        prime != PRIME
        and int(subprocess.check_output([args.solver, "--prime"], text=True)) != prime
    ):
        raise ValueError("incompatible native modulus")
    directory = args.directory or ROOT / ("single_planar_" + str(prime))
    directory.mkdir(exist_ok=True, parents=True)
    name = args.shape
    value = integrate_ledger(
        args.solver,
        directory / (name + ".tsv"),
        f"# 4+5 single planar four-cycle x 665280 reflection v1 {name} {prime}",
        jobs,
        lambda n, o, d: payload(n, shape, o, d, prime),
        prime,
    )
    scalar = next(den for sh, den, _ in SHAPES if sh == shape)
    result = value * pow(2 * scalar, -1, prime) % prime
    print(
        "shape",
        name,
        "weighted partition sum",
        value,
        "coefficient of pi",
        result,
        flush=True,
    )
    return result


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shape", choices=["166", "256", "346", "355", "445"])
    parser.add_argument("--solver")
    parser.add_argument("--prime", type=int, default=PRIME)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--verify-planar", action="store_true")
    args = parser.parse_args()
    if args.verify_planar:
        verify_planar()
    run(args)
