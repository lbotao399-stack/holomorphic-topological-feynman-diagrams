"""Two exact external geometries for the complete five-shape class detector."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import prod
import os
from pathlib import Path
import subprocess

from detector_node import single_fields
from single_projection import input_geometry, residual_det
from theta_denominators import compositions
from ordered_shapes import determinant

ROOT = Path(__file__).resolve().parent
P = 1009
# Each value includes the full 94-row transverse detector and the signed
# determinant with its eight pure wave scales divided out. Sources and graph
# symmetries enter once per shape when combining the two weighted geometries.
GEOMETRIES = (
    (1, (-38, -38, 29, 1), 63, Fraction(160, 10469)),
    (1, (-38, -38, 29, 2), 63, Fraction(-80, 10469)),
)


def verify():
    import sys

    parent = str(ROOT.parent)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    from mixed_graph_kernel import PLANAR, invariants, pure, verify_planar

    verify_planar()
    result = []
    for i in range(3):
        order = [j for j in range(3) if j != i] + [i, 3]
        forms = [pure(d, PLANAR[1][j]) for d, j in zip((3, 3, 4, 4), order)]
        result.append(
            tuple((-1) ** (2 - i) * Fraction(x, 20736) for x in invariants(forms))
        )
    assert result == [(-1, -78), (-1, -78), (16, -48)]
    s = (-38, -38, 29)
    vec = tuple(sum(t * r[j] for t, r in zip(s, result)) for j in range(2))
    assert vec == (540, 4536) and tuple(Fraction(-5, 54) * x for x in vec) == (
        -50,
        -420,
    )
    assert -Fraction(5 * 20736, 54 * prod(s) * 6) == GEOMETRIES[1][-1]

    pairs = list(combinations(range(4), 2))
    checked = 0
    for multiplicity in compositions(7, 6):
        degree = [
            sum(n for p, n in zip(pairs, multiplicity) if j in p) for j in range(4)
        ]
        # Each pure Hamiltonian wave contributes one polarization and one
        # extracted determinant column. The mark contributes its polarization.
        power = tuple(d - e for d, e in zip(degree, (2, 2, 2, 1)))
        if min(power[:3]) < 1 or power[3] < 0:
            continue

        def sample(branch, scales, extract=False):
            q = PLANAR[branch]
            v = 0
            for perm in permutations(range(3)):
                o = (*perm, 3)
                sign = (-1) ** sum(
                    x > y for k, x in enumerate(perm) for y in perm[k + 1 :]
                )
                e = [0] * 4
                for j in range(4):
                    e[o[j]] = power[j]
                if extract and tuple(e) != (1, 1, 2, 3):
                    continue
                invari = prod(
                    (q[o[i]][0] * q[o[j]][1] - q[o[i]][1] * q[o[j]][0]) ** n
                    for (i, j), n in zip(pairs, multiplicity)
                )
                v += sign * invari * prod(t**a for t, a in zip(scales, e))
            return v

        expected = -1313280 * sample(0, (1, 1, 1, 1), True) + 17280 * sample(
            1, (1, 1, 1, 1), True
        )
        value = sum(w * sample(b, sc) for b, sc, m, w in GEOMETRIES)
        assert expected == value, (multiplicity, expected, value)
        checked += 1
    print(
        "quartic-placement invariant vectors",
        result,
        "planar bracket monomials",
        checked,
        flush=True,
    )


def field(q):
    return q.numerator * pow(q.denominator, -1, P) % P


def shape166(args):
    import fcntl

    if int(subprocess.check_output([args.solver, "--prime"], text=True)) != P:
        raise ValueError("expected 1009")
    file = args.checkpoint
    hdr = "# complete-class two-geometry theta166 1009 mixed_single_weights.tsv"
    with (
        file.open("a+") as out,
        subprocess.Popen(
            [args.solver],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
        ) as proc,
    ):
        fcntl.flock(out, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.seek(0)
        lines = out.read().splitlines()
        if lines and lines[0] != hdr:
            raise ValueError("incompatible checkpoint")
        if not lines:
            print(hdr, file=out, flush=True)
        saved = {}
        for r in lines[1:]:
            s, a, v = r.split()
            saved[s, int(a)] = int(v)
        sums = Counter()
        try:
            for group in combinations(range(9), 5):
                order = group + tuple(j for j in range(10) if j not in group)
                digits = "".join(map(str, order))
                for br, sc, m, _ in GEOMETRIES:
                    a = sc[-1]
                    if (digits, a) not in saved:
                        data, cols = single_fields(br, 0, 0, 0, sc, m)
                        det = residual_det(order, br, cols)
                        assert det == determinant(cols, (1, 6, 6), order)
                        if det % P:
                            proc.stdin.write(
                                input_geometry(a, order, br, sc, m, 120, P)
                            )
                            proc.stdin.flush()
                            row = proc.stdout.readline().split()
                            if len(row) != 2 or int(row[0]) != a:
                                raise RuntimeError((digits, a, row))
                            v = int(row[1])
                        else:
                            v = 0
                        saved[digits, a] = v
                        print(digits, a, v, file=out, flush=True)
                    sums[a] = (sums[a] + saved[digits, a]) % P
            for a, v in sums.items():
                print(
                    "166 scale",
                    a,
                    "normalized sum",
                    v,
                    "coefficient of pi",
                    v * pow(2 * 158400, -1, P) % P,
                    flush=True,
                )
            result = (
                sum(field(w) * sums[sc[-1]] for b, sc, m, w in GEOMETRIES)
                * pow(2 * 158400, -1, P)
                % P
            )
            print(
                "166 contribution to two-geometry complete class sum",
                result,
                "mod",
                P,
                flush=True,
            )
        finally:
            proc.stdin.close()
            if proc.wait():
                raise RuntimeError("partition solver stopped")


def complete_shapes(directory):
    from complete_nodes import records
    from theta_denominators import SHAPES

    complete = True
    totals = Counter()
    for shape, denominator, _ in SHAPES:
        name = "".join(map(str, shape))
        for _, scales, _, _ in GEOMETRIES:
            scale = scales[-1]
            path = directory / f"{name}_{scale}.tsv"
            expected = {"".join(map(str, order)) for order, _ in records(shape, scale)}
            lines = path.read_text().splitlines() if path.exists() else []
            if (
                lines
                and lines[0] != f"# full mixed-current class v1 {name} {scale} {P}"
            ):
                raise ValueError(f"incompatible ledger {path}")
            values = {}
            for line in lines[1:]:
                order, value = line.split()
                value = int(value)
                if order in values or order not in expected or not 0 <= value < P:
                    raise ValueError(f"invalid partition {path} {line}")
                values[order] = value
            if len(values) != len(expected):
                print(
                    name,
                    scale,
                    "completed partitions",
                    len(values),
                    "of",
                    len(expected),
                )
                complete = False
                continue
            value = sum(values.values()) * pow(2 * denominator, -1, P) % P
            print(
                name, scale, "all partitions", len(values), "coefficient of pi", value
            )
            totals[scale] += value
    if complete:
        print(
            "all five shapes, projected coefficient of pi",
            sum(field(w) * totals[scales[-1]] for _, scales, _, w in GEOMETRIES) % P,
        )
    return complete


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("solver", nargs="?")
    p.add_argument("--checkpoint", type=Path, default=ROOT / "total_projection_166.tsv")
    p.add_argument("--complete-shapes", type=Path)
    args = p.parse_args()
    verify()
    if args.solver:
        shape166(args)
    if args.complete_shapes is not None:
        complete_shapes(args.complete_shapes)
