"""Project the 94-row mixed detector on a single ordered graph of any shape."""

import argparse
import os
import subprocess
from collections import Counter
from functools import cache
from pathlib import Path

from detector_node import single_fields
from exact_scalar import MASK, fixture, graph_factors
from ordered_shapes import determinant, paths
from single_projection import cases, effective_rows
from theta_denominators import SHAPES

PRIME = 1009


@cache
def determinant_table(shape, order, br):
    _, cols = single_fields(br, 0, 0, 0, (1, 1, 1, 1), 63)
    table = []
    for i in range(5):
        for j in range(5):
            probe = list(cols)
            probe[3] = tuple(int(k == i) for k in range(5))
            probe[9] = tuple(int(k == j) for k in range(5))
            d = determinant(probe, shape, order)
            if d:
                table.append((i, j, d))
    return table


def geometry(i, shape, order, data, det, weighted_polars, prime=PRIME):
    if not det % prime or not weighted_polars:
        return ""
    p, pol = data
    lines = [(i, det % prime, *shape)]
    groups = paths(shape, order)
    polar = [{} for _ in groups]
    weights = Counter()
    for w, u in weighted_polars:
        choices = []
        for h, group in enumerate(groups):
            v = tuple(u[j] for j in group)
            choices.append(polar[h].setdefault(v, len(polar[h])))
        weights[tuple(choices)] += w
    for h, group in enumerate(groups):
        lines.append(group)
        lines.append(tuple(a for j in group for a in p[j]))
        lines.append((len(polar[h]),))
        for row in polar[h]:
            lines.append(tuple(x for v in row for x in v))
    rows = [(*r, w % prime) for r, w in weights.items() if w % prime]
    lines.append((len(rows),))
    lines.extend(rows)
    for s, _, ends in graph_factors(data, order, (), True, True, shape):
        lines.append((s,))
        for e in ends:
            rows = [
                (25 * x + 5 * y + z, k & MASK, k >> 20, w % prime)
                for (x, y, z), q in e.items()
                for k, w in q.items()
                if w % prime
            ]
            lines.append((len(rows),))
            lines.extend(rows)
    return "".join(" ".join(map(str, row)) + "\n" for row in lines)


def dataset(i, shape, order, br, sc, m):
    data, cols = single_fields(br, 0, 0, 0, sc, m)
    d = sum(
        c * cols[3][j] * cols[9][k] for j, k, c in determinant_table(shape, order, br)
    )
    if not d % PRIME:
        return ""
    d *= (-1) ** (9 - order.index(9))
    weights = [
        (w, single_fields(br, bit, 0, eta, sc, m)[0][1])
        for w, bit, eta in effective_rows(br, sc, m, PRIME)
    ]
    return geometry(i, shape, order, data, d, weights)


def run(args):
    labels = ("0146723589", "0146723589", "0146723589", "6901472358", "0691472358")
    chosen = tuple(map(int, args.shape))
    order = tuple(
        map(
            int,
            args.order
            or labels[next(i for i, (s, _, _) in enumerate(SHAPES) if s == chosen)],
        )
    )
    if sorted(order) != list(range(10)):
        raise ValueError("expected ten distinct labels")
    with subprocess.Popen(
        [args.solver],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
    ) as proc:

        def eval_node(payload, i):
            if not payload:
                return 0
            proc.stdin.write(payload)
            proc.stdin.flush()
            reply = proc.stdout.readline().split()
            if len(reply) != 2 or int(reply[0]) != i:
                raise RuntimeError((i, reply))
            return int(reply[1])

        try:
            if args.fixture:
                data = fixture()
                d = determinant(data[0], chosen, order) * (-1) ** (9 - order.index(9))
                value = eval_node(
                    geometry(0, chosen, order, data, d, [(1, data[1])]), 0
                )
                print("fixture", args.shape, "det", d, "value", value, "mod", PRIME)
                return
            directory = args.directory
            directory.mkdir(exist_ok=True, parents=True)
            path = directory / (args.shape + "_" + "".join(map(str, order)) + ".tsv")
            import fcntl

            with path.open("a+") as out:
                fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
                hdr = (
                    "# mixed_single_weights.tsv ordered "
                    + args.shape
                    + " "
                    + "".join(map(str, order))
                    + " 1009"
                )
                out.seek(0)
                existing = out.read().splitlines()
                if existing and existing[0] != hdr:
                    raise ValueError("incompatible checkpoint")
                if not existing:
                    print(hdr, file=out, flush=True)
                saved = dict(tuple(map(int, row.split())) for row in existing[1:])
                total = 0
                for i, br, sc, m, w in cases():
                    if args.case is not None and i != args.case:
                        continue
                    if i not in saved:
                        saved[i] = eval_node(dataset(i, chosen, order, br, sc, m), i)
                        print(i, saved[i], file=out, flush=True)
                    total = (total + w * saved[i]) % PRIME
                if args.case is None:
                    print(
                        "complete ordered projected mean",
                        args.shape,
                        "".join(map(str, order)),
                        total * pow(1536, -1, PRIME) % PRIME,
                        "mod",
                        PRIME,
                    )
                else:
                    print("case", args.case, saved[args.case], "mod", PRIME)
        finally:
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError("native integrator failed")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("solver")
    p.add_argument("shape", choices=["".join(map(str, s)) for s, _, _ in SHAPES])
    p.add_argument("--fixture", action="store_true")
    p.add_argument("--order")
    p.add_argument("--case", type=int)
    p.add_argument(
        "--directory",
        type=Path,
        default=Path(__file__).resolve().parent / "ordered_detector",
    )
    run(p.parse_args())
