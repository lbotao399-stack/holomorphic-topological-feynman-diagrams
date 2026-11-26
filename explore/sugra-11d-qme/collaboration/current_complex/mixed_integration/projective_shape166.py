"""Sum S4 orbits of complete detector partitions with exact planar quadrature."""

import argparse
import os
import signal
import subprocess
from collections import defaultdict
from functools import cache
from itertools import combinations, product
from pathlib import Path

from detector_node import single_fields
from single_projection import effective_rows, input_geometry, residual_det

PRIME = 1009
ROOT = Path(__file__).resolve().parent
TABLE = "166 single-momentum 1009 projective-v1"


def balanced(x):
    return (x + PRIME // 2) % PRIME - PRIME // 2


@cache
def planar_nodes():
    omega = next(x for x in range(2, PRIME) if pow(x, 3, PRIME) == 1)
    imag = next(x for x in range(2, PRIME) if x * x % PRIME == PRIME - 1)
    nodes = []
    for j in range(4):
        z = pow(imag, j, PRIME)
        # On branch zero the two e1 waves alternate, of residual degree 3.
        n = pow(4 * (omega**2 - omega) * z**3, -1, PRIME)
        nodes.append((0, (1, 1, balanced(omega), balanced(z)), n * -1313280 % PRIME))
        # Project powers 2,3 in x2,x3. The remaining symmetric quadratic at
        # (1,i) loses both square terms and leaves i times the mixed coefficient.
        for k in range(3):
            t = pow(omega, k, PRIME)
            n = pow(12 * imag * t**2 * z**3, -1, PRIME)
            nodes.append(
                (1, (1, balanced(imag), balanced(t), balanced(z)), n * 17280 % PRIME)
            )
    return tuple(nodes)


def monomial_check():
    # Antisymmetrize only the reflection that the appropriate branch uses.
    for exp in product(range(8), repeat=4):
        if sum(exp) != 7:
            continue
        a, b, c, d = exp
        for br in (0, 1):
            if br == 0 and (a + c != 3 or b + d != 4):
                continue
            other = (c, b, a, d) if br == 0 else (b, a, c, d)
            # Swapping 0,1 in branch one reverses the marked Hamiltonian
            # direction, giving F(x1,x0,x2,-x3) = -F(x0,x1,x2,x3).
            sgn = -1 if not br else (-1) ** (d + 1)
            desired = int(exp == (1, 1, 2, 3)) + sgn * int(other == (1, 1, 2, 3))
            actual = 0
            for branch, scales, weight in planar_nodes():
                if branch == br:

                    def val(powers):
                        v = 1
                        for s, e in zip(scales, powers):
                            v = v * pow(s, e, PRIME) % PRIME
                        return v

                    actual += weight * (val(exp) + sgn * val(other))
            if actual % PRIME != desired * (-1313280 if br == 0 else 17280) % PRIME:
                raise AssertionError((br, exp, other, desired, actual))


def partition_orbits():
    groups = defaultdict(list)
    for left in combinations(range(9), 5):
        order = left + tuple(j for j in range(10) if j not in left)
        pure = "".join(str(j) for j in left if j >= 4)
        groups[pure].append(order)
    return groups


@cache
def geometry(node, mask):
    br, scales, _ = planar_nodes()[node]
    _, columns = single_fields(br, 0, 0, 0, scales, mask)
    rows = effective_rows(br, scales, mask, PRIME)
    return columns, bool(rows)


def tasks(orders):
    for order in orders:
        digits = "".join(map(str, order))
        for node, (br, sc, w) in enumerate(planar_nodes()):
            for m in range(63, 31, -1):
                cols, active = geometry(node, m)
                if active and residual_det(order, br, cols) % PRIME:
                    yield (digits, node, m), w * (-1) ** m.bit_count() % PRIME


def solve(args):
    monomial_check()
    orbits = partition_orbits()
    if args.check:
        counts = {k: sum(1 for _ in tasks(v)) for k, v in orbits.items()}
        print(
            "Monomial projection verified. Active nodes",
            sum(counts.values()),
            "in",
            sum(bool(v) for v in counts.values()),
            "orbits",
            sorted(counts.items(), key=lambda v: (v[1], v[0])),
        )
        return
    if args.all:
        totals = []
        for key, orders in sorted(
            orbits.items(), key=lambda item: (len(item[1]), item[0])
        ):
            args.orbit = key
            totals.append(solve_orbit(args, orders))
        normalized = sum(totals) * pow(2, -1, PRIME) % PRIME
        print(
            "full shape with graph symmetry",
            normalized,
            "coefficient of pi",
            normalized * pow(158400, -1, PRIME) % PRIME,
            "mod",
            PRIME,
            flush=True,
        )
        return
    if args.orbit not in orbits:
        raise ValueError("choose a valid transverse subset")
    return solve_orbit(args, orbits[args.orbit])


def solve_orbit(args, orders):
    if int(subprocess.check_output([args.solver, "--prime"], text=True)) != PRIME:
        raise ValueError("incompatible prime")
    pending = dict(tasks(orders))
    args.directory.mkdir(exist_ok=True, parents=True)
    checkpoint = args.directory / (args.orbit + ".tsv")
    import fcntl

    with checkpoint.open("a+") as out:
        fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.seek(0)
        hdr = "# " + TABLE + " " + args.orbit
        lines = out.read().splitlines()
        if lines and lines[0] != hdr:
            raise ValueError("incompatible quadrature checkpoint")
        if not lines:
            print(hdr, file=out, flush=True)
        total = 0
        for row in lines[1:]:
            digits, n, m, v = row.split()
            total = (total + pending.pop((digits, int(n), int(m))) * int(v)) % PRIME
        remaining = len(pending)
        if remaining:
            with subprocess.Popen(
                [args.solver],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                text=True,
                env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
            ) as proc:
                try:
                    for (digits, n, m), w in pending.items():
                        br, sc, _ = planar_nodes()[n]
                        order = tuple(map(int, digits))
                        job = n * 64 + m
                        payload = input_geometry(job, order, br, sc, m, 120, PRIME)
                        proc.stdin.write(payload)
                        proc.stdin.flush()
                        reply = proc.stdout.readline().split()
                        if len(reply) != 2 or int(reply[0]) != job:
                            raise RuntimeError((digits, n, m, reply))
                        v = int(reply[1])
                        total = (total + w * v) % PRIME
                        print(digits, n, m, v, file=out, flush=True)
                finally:
                    proc.stdin.close()
                    if proc.wait():
                        raise RuntimeError("native detector stopped")
        print(
            "orbit",
            args.orbit,
            "partitions",
            len(orders),
            "computed",
            remaining,
            "normalized coefficient",
            total,
            "mod",
            PRIME,
            flush=True,
        )
        return total


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("solver", nargs="?", default="")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--orbit")
    parser.add_argument("--all", action="store_true")
    parser.add_argument(
        "--directory", type=Path, default=ROOT / "single_shape_projective"
    )
    args = parser.parse_args()

    def interrupt(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    solve(args)
