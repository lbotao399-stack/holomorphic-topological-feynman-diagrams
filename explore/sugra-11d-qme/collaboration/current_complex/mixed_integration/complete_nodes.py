"""Sum complete path partitions with the two all-active geometries of the class test."""

import argparse
import fcntl
from itertools import combinations
from pathlib import Path
import os
import signal
import subprocess
from functools import cache

from detector_node import single_fields
from ordered_detector import geometry
from ordered_shapes import determinant
from single_projection import effective_rows, input_geometry
from theta_denominators import SHAPES

P = 1009
ROOT = Path(__file__).resolve().parent


def partitions(shape):
    labels = tuple(range(10))
    for first in combinations(labels, shape[0] - 1):
        other = tuple(i for i in labels if i not in first)
        for second in combinations(other, shape[1] - 1):
            third = tuple(i for i in other if i not in second)
            if shape == (1, 6, 6):
                if 9 not in third:
                    continue
            elif (shape[0] == shape[1] and first > second) or (
                shape[1] == shape[2] and second > third
            ):
                continue
            yield first + second + third


@cache
def fixture(s, prime=P):
    scales = (-38, -38, 29, s)
    data, cols = single_fields(1, 0, 0, 0, scales, 63)
    rows = [
        (w, single_fields(1, bits, 0, eta, scales, 63)[0][1])
        for w, bits, eta in effective_rows(1, scales, 63, prime)
    ]
    return scales, data, cols, rows


def records(shape, s, prime=P):
    scales, data, cols, rows = fixture(s, prime)
    for order in partitions(shape):
        det = determinant(cols, shape, order) % prime
        if det:
            yield order, det


def payload(case, shape, order, det, s, special, prime=P):
    scales, data, _, rows = fixture(s, prime)
    if special:
        return input_geometry(case, order, 1, scales, 63, 120, prime)
    signed = det * (-1) ** (9 - order.index(9))
    return geometry(case, shape, order, data, signed, rows, prime)


def run(a):
    prime = getattr(a, "prime", P)
    if (
        prime != P
        and int(subprocess.check_output([a.solver, "--prime"], text=True)) != prime
    ):
        raise ValueError("incompatible native modulus")
    shape = tuple(map(int, a.shape))
    scalar = next(d for sh, d, _ in SHAPES if sh == shape)
    special = shape == (1, 6, 6) and not a.general
    a.directory.mkdir(parents=True, exist_ok=True)
    path = a.directory / f"{a.shape}_{a.scale}.tsv"
    jobs = list(records(shape, a.scale, prime))
    if a.count:
        print(a.shape, a.scale, len(jobs))
        return
    value = integrate_ledger(
        a.solver,
        path,
        f"# full mixed-current class v1 {a.shape} {a.scale} {prime}",
        jobs,
        lambda n, order, det: payload(n, shape, order, det, a.scale, special, prime),
        prime,
    )
    print(
        a.shape,
        a.scale,
        "unweighted partition sum",
        value,
        "coefficient of pi",
        value * pow(2 * scalar, -1, prime) % prime,
    )


def integrate_ledger(solver, path, header, jobs, encode, prime):
    with path.open("a+") as ledger:
        fcntl.flock(ledger, fcntl.LOCK_EX | fcntl.LOCK_NB)
        ledger.seek(0)
        lines = ledger.read().splitlines()
        if lines and lines[0] != header:
            raise ValueError(path)
        if not lines:
            print(header, file=ledger, flush=True)
        saved = {}
        for row in lines[1:]:
            key, v = row.split()
            if key in saved:
                raise ValueError("duplicate")
            saved[key] = int(v)
        pending = [
            (n, order, det)
            for n, (order, det) in enumerate(jobs)
            if "".join(map(str, order)) not in saved
        ]
        # Reap and renew the native process between batches so its allocation
        # high-water mark is bounded during multi-hour integrations.
        for start in range(0, len(pending), 12):
            with subprocess.Popen(
                [solver],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                text=True,
                env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
            ) as proc:
                try:
                    for n, order, det in pending[start : start + 12]:
                        key = "".join(map(str, order))
                        proc.stdin.write(encode(n, order, det))
                        proc.stdin.flush()
                        row = proc.stdout.readline().split()
                        if len(row) != 2 or int(row[0]) != n:
                            raise RuntimeError((key, row))
                        saved[key] = int(row[1])
                        print(key, saved[key], file=ledger, flush=True)
                finally:
                    proc.stdin.close()
                    if proc.wait():
                        raise RuntimeError("native kernel stopped")
        if set(saved) != {"".join(map(str, order)) for order, _ in jobs}:
            raise ValueError("unexpected partition")
        return sum(saved.values()) % prime


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    p = argparse.ArgumentParser()
    p.add_argument("solver")
    p.add_argument("shape", choices=["166", "256", "346", "355", "445"])
    p.add_argument("scale", type=int, choices=(1, 2))
    p.add_argument("--prime", type=int, default=P)
    p.add_argument("--general", action="store_true")
    p.add_argument("--count", action="store_true")
    p.add_argument("--directory", type=Path, default=ROOT / "complete_shapes")
    run(p.parse_args())
