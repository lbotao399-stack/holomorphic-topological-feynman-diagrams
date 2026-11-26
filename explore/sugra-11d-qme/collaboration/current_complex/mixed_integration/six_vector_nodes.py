"""The planar-current test with six transverse vector fields on ordered theta partitions."""

from collections import defaultdict
from functools import cache
from pathlib import Path
import argparse
import fcntl
import signal
import os
import subprocess

from complete_nodes import partitions
from ordered_detector import geometry
from ordered_shapes import determinant
from theta_denominators import SHAPES

PRIME = 1009
BASE = Path(__file__).resolve().parent
TRANSVERSE = ((1, 0, 0), (1, 1, 1), (0, 1, 0), (1, 1, -1), (0, 0, 1), (1, 0, 1))
FUNCTION = ((1, 0, -1), (0, 1, 1))
X = (1, 0, 0, 0, 0)
Y = (0, 1, 0, 0, 0)


def fields(mask, choice):
    # The marked planar momentum is carried by its determinant column alone.
    momentum = [X, Y, (0, 0, *FUNCTION[choice]), *((0, 0, *p) for p in TRANSVERSE)]
    momentum.append(tuple(-sum(c) for c in zip(*momentum)))
    polar = [Y, (-1, 0, 0, 0, 0), Y]
    for j, p in enumerate(TRANSVERSE):
        axis = j // 2
        other = (axis + 1 + ((mask >> j) & 1)) % 3
        u = [0] * 3
        u[other], u[axis] = 1, -p[other]
        assert sum(v * w for v, w in zip(p, u)) == 0
        polar.append((0, 0, *u))
    polar.append(X)
    columns = [X, Y, X, *((0, 0, *p) for p in TRANSVERSE), (0, -1, 0, 0, 0)]
    return (momentum, polar), columns


@cache
def fixture():
    samples = defaultdict(list)
    for line in (BASE.parent / "mixed_six_weights.tsv").read_text().splitlines():
        weight, choice, mask = map(int, line.split())
        data, cols = fields(mask, choice)
        samples[choice].append((weight, data[1]))
    return samples


def records(shape):
    _, columns = fields(0, 0)
    for order in partitions(shape):
        d = determinant(columns, shape, order)
        if d:
            yield order, d


def payload(i, shape, order, det, choice, prime=PRIME):
    data, _ = fields(0, choice)
    # The 24-term planar cycle is -7 copies of the pure Hamiltonian
    # fixture. The Taylor-to-Hamiltonian multiplier is 3! 6! 2!.
    signed = det * (-60480) * (-1) ** (9 - order.index(9))
    return geometry(i, shape, order, data, signed, fixture()[choice], prime)


def run(args, shape, den, prime=PRIME):
    jobs = list(records(shape))
    name = "".join(map(str, shape))
    print(name, "determinant partitions", len(jobs), flush=True)
    if not args.solver:
        return
    args.directory.mkdir(exist_ok=True, parents=True)
    path = args.directory / f"{name}.tsv"
    with path.open("a+") as out:
        fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.seek(0)
        lines = out.read().splitlines()
        header = f"# six-vector class {prime} {name}: order, function 0, function 1"
        if lines and lines[0] != header:
            raise ValueError("incompatible checkpoint")
        if not lines:
            print(header, file=out, flush=True)
        saved = {}
        for line in lines[1:]:
            k, a, b = line.split()
            if k in saved:
                raise ValueError("duplicate partition")
            saved[k] = (int(a), int(b))
        pending = [
            (i, (order, det))
            for i, (order, det) in enumerate(jobs[: args.limit])
            if "".join(map(str, order)) not in saved
        ]
        for start in range(0, len(pending), 12):
            with subprocess.Popen(
                [args.solver],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                text=True,
                env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
            ) as proc:
                try:
                    for i, (order, det) in pending[start : start + 12]:
                        key = "".join(map(str, order))
                        if key in saved:
                            continue
                        values = []
                        for choice in (0, 1):
                            node = 2 * i + choice
                            proc.stdin.write(
                                payload(node, shape, order, det, choice, prime)
                            )
                            proc.stdin.flush()
                            row = proc.stdout.readline().split()
                            if len(row) != 2 or int(row[0]) != node:
                                raise ArithmeticError((key, choice, row))
                            values.append(int(row[1]))
                        saved[key] = values
                        print(key, *values, sep="\t", file=out, flush=True)
                finally:
                    proc.stdin.close()
                    if proc.wait():
                        raise RuntimeError("native integral failed")
        if len(saved) == len(jobs):
            assert set(saved) == {"".join(map(str, order)) for order, _ in jobs}
            totals = [sum(row[k] for row in saved.values()) % prime for k in (0, 1)]
            result = sum(totals) * pow(2 * den, -1, prime) % prime
            print(
                name,
                "order sums",
                totals,
                "coefficient divided by pi",
                result,
                flush=True,
            )
            return result


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shape", choices=["166", "256", "346", "355", "445"])
    parser.add_argument("--directory", type=Path, default=BASE / "six_vector_shapes")
    args = parser.parse_args()
    for shape, den, _ in SHAPES:
        if args.shape is None or args.shape == "".join(map(str, shape)):
            run(args, shape, den)
