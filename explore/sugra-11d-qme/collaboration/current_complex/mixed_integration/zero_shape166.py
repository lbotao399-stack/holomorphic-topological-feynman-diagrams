"""Integrate the mixed 166 detector, summing zero-wave assignments by determinant."""

import argparse
import fcntl
import os
import signal
from itertools import combinations
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

from projective_shape166 import (
    PRIME,
    ROOT,
    TABLE,
    geometry,
    monomial_check,
    partition_orbits,
    planar_nodes,
)
from single_projection import input_geometry, residual_det

HEADER = "# mixed 166 single-momentum determinant sums v1"


def groups(orders, node, mask):
    branch, _, _ = planar_nodes()[node]
    columns, active = geometry(node, mask)
    if not active:
        return {}
    zero = {4 + i for i in range(5) if not mask >> i & 1}
    grouped = defaultdict(list)
    for order in orders:
        det = residual_det(order, branch, columns) % PRIME
        if det:
            left = tuple(j for j in order[:5] if j not in zero)
            grouped[left].append((order, det))
    return grouped


def exterior_determinant(columns):
    """Finite-field determinant by row elimination on the column transpose."""
    matrix = [[v % PRIME for v in column] for column in columns]
    result = 1
    for i in range(10):
        pivot = next((j for j in range(i, 10) if matrix[j][i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            matrix[pivot], matrix[i] = matrix[i], matrix[pivot]
            result = -result
        factor = matrix[i][i]
        result = result * factor % PRIME
        inv = pow(factor, -1, PRIME)
        for j in range(i + 1, 10):
            scale = matrix[j][i] * inv % PRIME
            if scale:
                matrix[j][i:] = [
                    (a - scale * b) % PRIME
                    for a, b in zip(matrix[j][i:], matrix[i][i:])
                ]
    return result


def verify(solver):
    orders = sum(partition_orbits().values(), [])
    tested = shared = 0
    # Include a full-rank and a deficient-rank geometry, both planar branches,
    # and strata with two, three and four inactive vectors.
    with subprocess.Popen(
        [solver],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
    ) as proc:
        try:
            for node, mask in ((0, 58), (1, 60), (8, 56), (5, 38), (1, 33)):
                columns, active = geometry(node, mask)
                if not active:
                    continue
                zero = {4 + j for j in range(5) if not mask >> j & 1}
                br, _, _ = planar_nodes()[node]
                all_groups = groups(orders, node, mask)
                for k in range(6):
                    for left in combinations(sorted(set(range(9)) - zero), k):
                        if not 0 <= 5 - k <= len(zero):
                            continue
                        combined = []
                        for j, col in enumerate(columns):
                            first = col if j in left or j in zero else (0,) * 5
                            second = col if j not in left else (0,) * 5
                            combined.append((*first, *second))
                        expected = exterior_determinant(combined)
                        members = all_groups.get(left, [])
                        actual = sum(d for _, d in members) % PRIME
                        if expected != actual:
                            raise ArithmeticError((br, mask, left, actual, expected))
                        tested += 1
                for members in all_groups.values():
                    if len({d for _, d in members}) <= 1:
                        continue
                    wick = {
                        value(proc, order, node, mask) * pow(d, -1, PRIME) % PRIME
                        for order, d in members
                    }
                    if len(wick) != 1:
                        raise ArithmeticError((node, mask, members, wick))
                    shared += len(members)
                    break
        finally:
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError("native detector stopped")
    print("expanded determinants", tested, "equal-weight path values", shared)


def previous(directory):
    saved = {}
    for path in directory.glob("*.tsv"):
        with path.open() as handle:
            if handle.readline().strip() != f"# {TABLE} {path.stem}":
                raise ValueError(path)
            for line in handle:
                if not line.endswith("\n"):  # Another worker can append here.
                    break
                digits, n, m, v = line.split()
                saved[tuple(map(int, digits)), int(n), int(m)] = int(v)
    return saved


def value(proc, order, node, mask):
    br, scales, _ = planar_nodes()[node]
    idx = 64 * node + mask
    proc.stdin.write(input_geometry(idx, order, br, scales, mask, 120, PRIME))
    proc.stdin.flush()
    res = proc.stdout.readline().split()
    if len(res) != 2 or int(res[0]) != idx:
        raise RuntimeError((order, node, mask, res))
    return int(res[1])


def compute(args):
    monomial_check()
    if args.verify:
        return verify(args.solver)
    orders = sum(partition_orbits().values(), [])
    if args.orbit:
        orders = partition_orbits()[args.orbit]
    counts = Counter()
    batches = []
    nodes = range(len(planar_nodes())) if args.node is None else (args.node,)
    for n in nodes:
        for m in range(63, 31, -1):
            for members in groups(orders, n, m).values():
                zeros = 5 - (m & 31).bit_count()
                counts["original", zeros] += len(members)
                summed = sum(det for _, det in members) % PRIME
                if summed:
                    counts["remaining", zeros] += 1
                    batches.append((n, m, members, summed))
                if not args.orbit and zeros >= 4 and summed:
                    raise ArithmeticError("diagonal transverse determinant")
    if args.count:
        print([(i, counts["original", i], counts["remaining", i]) for i in range(6)])
        return
    if int(subprocess.check_output([args.solver, "--prime"], text=True)) != PRIME:
        raise ValueError("solver modulus")
    args.directory.mkdir(exist_ok=True, parents=True)
    saved = previous(args.projective)
    total = 0
    node_files = {}
    try:
        with subprocess.Popen(
            [args.solver],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
        ) as proc:
            try:
                for n, m, members, summed in batches:
                    if n not in node_files:
                        path = args.directory / f"{args.orbit or 'all'}-{n:02d}.tsv"
                        handle = path.open("a+")
                        try:
                            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                            handle.seek(0)
                            hdr = f"{HEADER} {args.orbit or 'all'} {n} {PRIME}"
                            lines = handle.read().splitlines()
                            if lines and lines[0] != hdr:
                                raise ValueError(path)
                            if not lines:
                                print(hdr, file=handle, flush=True)
                            for line in lines[1:]:
                                digits, mask, result = line.split()
                                saved[tuple(map(int, digits)), n, int(mask)] = int(
                                    result
                                )
                        except BaseException:
                            handle.close()
                            raise
                        node_files[n] = handle
                    available = {
                        saved[order, n, m] * pow(det, -1, PRIME) % PRIME
                        for order, det in members
                        if (order, n, m) in saved
                    }
                    if len(available) > 1:
                        raise ArithmeticError(
                            ("unequal normalized Wick tensors", n, m, members)
                        )
                    if available:
                        vk = available.pop()
                    else:
                        order, det = members[0]
                        raw = value(proc, order, n, m)
                        print(
                            "".join(map(str, order)),
                            m,
                            raw,
                            file=node_files[n],
                            flush=True,
                        )
                        vk = raw * pow(det, -1, PRIME) % PRIME
                    total += planar_nodes()[n][2] * (-1) ** m.bit_count() * summed * vk
            finally:
                proc.stdin.close()
                if proc.wait():
                    raise RuntimeError("native detector stopped")
    finally:
        for handle in node_files.values():
            handle.close()
    print("normalized sum", total % PRIME, "mod", PRIME)
    if args.node is None and not args.orbit:
        print(
            "coefficient of pi with graph symmetry",
            total * pow(2 * 158400, -1, PRIME) % PRIME,
        )


if __name__ == "__main__":

    def interrupt(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    p = argparse.ArgumentParser()
    p.add_argument("solver", nargs="?", default="")
    p.add_argument("--orbit")
    p.add_argument("--node", type=int, choices=range(16))
    p.add_argument("--count", action="store_true")
    p.add_argument("--verify", action="store_true")
    p.add_argument("--directory", type=Path, default=ROOT / "single_shape_zero")
    p.add_argument("--projective", type=Path, default=ROOT / "single_shape_projective")
    compute(p.parse_args())
