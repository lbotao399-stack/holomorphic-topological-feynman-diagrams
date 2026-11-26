"""Project one whole 166 partition onto the single-geometry mixed detector."""

import argparse
from collections import Counter, defaultdict
from itertools import product
from functools import cache
from math import comb
from pathlib import Path
import os
import subprocess
from detector_factors import SINGLE_WEIGHTS, key
from detector_node import single_fields
from exact_scalar import MASK, collect, graph_factors, multiply
from partition166 import determinant

PRIME = 1009


@cache
def determinant_table(order, branch):
    _, columns = single_fields(branch, 0, 0, 0, (1, 1, 1, 1), 63)
    table = []
    for i in range(5):
        for j in range(5):
            probe = list(columns)
            probe[3] = tuple(int(k == i) for k in range(5))
            probe[9] = tuple(int(k == j) for k in range(5))
            v = determinant(probe, order)
            if v:
                table.append((i, j, v))
    return table


def residual_det(order, branch, columns):
    return sum(
        v * columns[3][i] * columns[9][j]
        for i, j, v in determinant_table(order, branch)
    )


@cache
def detector_rows():
    return tuple(
        tuple(map(int, line.split()))
        for line in SINGLE_WEIGHTS.read_text().splitlines()
    )


def effective_rows(branch, scales, switches, prime):
    (mom, u0), _ = single_fields(branch, 0, 0, 0, scales, switches)
    # Row-reduce the external covectors. At rank four eta must occupy the
    # distinguished volume slot and every ordinary polarization pairs only
    # with one of the momentum covectors.
    basis = {}
    for p in mom:
        v = list(p)
        for k, row in sorted(basis.items()):
            q = v[k] % prime
            if q:
                v = [(a - q * b) % prime for a, b in zip(v, row)]
        pivot = next((i for i, x in enumerate(v) if x % prime), None)
        if pivot is None:
            continue
        inv = pow(v[pivot], -1, prime)
        basis[pivot] = [x * inv % prime for x in v]
    if len(basis) == 5:
        return detector_rows()
    if len(basis) < 4:
        return ()
    free = next(i for i in range(5) if i not in basis)
    normal = [0] * 5
    normal[free] = 1
    for k, row in sorted(basis.items(), reverse=True):
        normal[k] = -sum(x * y for x, y in zip(row, normal)) % prime
    eta = next((i for i in range(3) if normal[2 + i]), None)
    if eta is None:
        return ()
    factor_eta = pow(normal[2 + eta], -1, prime)
    alternatives = single_fields(branch, 31, 0, 0, scales, switches)[0][1]
    factors = []
    for j in range(9):
        groups = {}
        mapping = []
        for bit, u in enumerate((u0[j], alternatives[j])):
            v = tuple(sum(x * y for x, y in zip(q, u)) % prime for q in mom)
            first = next((a for a in v if a), 0)
            if not first:
                mapping.append((bit, 0))
                continue
            fingerprint = tuple(x * pow(first, -1, prime) % prime for x in v)
            ref, amplitude = groups.setdefault(fingerprint, (bit, first))
            mapping.append((ref, first * pow(amplitude, -1, prime) % prime))
        if j < 4 and not mapping[0][1]:
            return ()
        if j >= 4:
            factors.append(mapping)
    out = Counter()
    for w, bits, g in detector_rows():
        weight = w * normal[2 + g] * factor_eta % prime
        new_bits = 0
        for j, options in enumerate(factors):
            bit, scale = options[(bits >> j) & 1]
            new_bits |= bit << j
            weight = weight * scale % prime
        out[new_bits] += weight
    return tuple((w % prime, b, eta) for b, w in out.items() if w % prime)


def input_geometry(case, order, branch, scales, switches, limit, prime=PRIME, rows=None):
    data, columns = single_fields(branch, 0, 0, 0, scales, switches)
    det = residual_det(order, branch, columns)
    if not det:
        return ""
    polar = [dict(), dict()]
    weights = Counter()
    effective = effective_rows(branch, scales, switches, prime) if rows is None else rows
    for w, bits, eta in effective:
        indices = []
        for h, group in enumerate((order[:5], order[5:])):
            k = key(group, bits, eta)
            if k not in polar[h]:
                values = single_fields(branch, bits, 0, eta, scales, switches)[0][1]
                polar[h][k] = (len(polar[h]), [values[l] for l in group])
            indices.append(polar[h][k][0])
        weights[tuple(indices)] += w
    stream = [(case, det if effective else 0, limit)]
    for h, group in enumerate((order[:5], order[5:])):
        stream.extend((group, [v for l in group for v in data[0][l]], [len(polar[h])]))
        for _, us in polar[h].values():
            stream.append([v for u in us for v in u])
    for i in range(len(polar[0])):
        stream.append([weights[i, j] % prime for j in range(len(polar[1]))])
    for flow, _, (l, r) in graph_factors(data, order, ()):
        ends = defaultdict(dict)
        for (i, j, k), a in l.items():
            for (h, s, t), b in r.items():
                if i == h:
                    collect(ends[5 * j + s, 5 * k + t], multiply(a, b), flow)
        rows = [
            (js, kt, q & MASK, q >> 20, w)
            for (js, kt), p in ends.items()
            for q, w in p.items()
            if w
        ]
        stream.append((len(rows),))
        stream.extend(rows)
    return "".join(" ".join(str(v) for v in line) + "\n" for line in stream)


def cases():
    for node, (i, j, k) in enumerate(product(range(2), range(3), range(4))):
        scales = 1, 2 * i - 1, 2 * j - 2, 2 * k - 3
        centered = (-1) ** (6 - i - j - k) * 2 * comb(2, j) * comb(3, k)
        for branch in (0, 1):
            for masks in range(63, -1, -1):
                index = (node * 2 + branch) * 64 + masks
                factor = (
                    centered
                    * (17280 if branch else -1313280)
                    * (-1) ** masks.bit_count()
                )
                yield index, branch, scales, masks, factor


def run(args):
    order = tuple(map(int, args.order))
    if (
        sorted(order) != list(range(10))
        or 9 not in order[5:]
        or not 1 <= args.limit <= 120
    ):
        raise ValueError("invalid partition or path order count")
    prime = int(subprocess.check_output([args.solver, "--prime"], text=True))
    header = f"# {args.order} {args.limit} {prime} {SINGLE_WEIGHTS.name}"
    saved = {}
    if args.checkpoint:
        if args.checkpoint.exists():
            lines = args.checkpoint.read_text().splitlines()
            if lines[:1] != [header]:
                raise ValueError("incompatible checkpoint")
            for line in lines[1:]:
                i, v = map(int, line.split())
                saved[i] = v
        else:
            args.checkpoint.write_text(header + "\n")
    proc = subprocess.Popen(
        [args.solver],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
    )
    answer = active = 0
    try:
        for i, br, sc, mask, w in cases():
            if args.case is not None and i != args.case:
                continue
            data, columns = single_fields(br, 0, 0, 0, sc, mask)
            if not residual_det(order, br, columns):
                continue
            active += 1
            if i not in saved:
                proc.stdin.write(
                    input_geometry(i, order, br, sc, mask, args.limit, prime)
                )
                proc.stdin.flush()
                row = proc.stdout.readline().split()
                if len(row) != 2 or int(row[0]) != i:
                    raise RuntimeError(f"case {i}, solver status {proc.poll()}")
                saved[i] = int(row[1])
                if args.checkpoint:
                    with args.checkpoint.open("a") as f:
                        print(i, saved[i], file=f)
            answer = (answer + w * saved[i]) % prime
            if args.case is not None:
                print("case", i, "detector", saved[i], "mod", prime)
    finally:
        proc.stdin.close()
        if proc.wait():
            raise RuntimeError(f"solver exit {proc.returncode}")
    if args.case is None:
        print(
            "partition",
            args.order,
            "surviving cases",
            active,
            "orders",
            args.limit**2,
            "projected normalized coefficient",
            answer * pow(1536, -1, prime) % prime,
            "mod",
            prime,
        )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("solver")
    p.add_argument("--order", default="0234516789")
    p.add_argument("--limit", type=int, default=120)
    p.add_argument("--case", type=int)
    p.add_argument("--checkpoint", type=Path)
    run(p.parse_args())
