"""Recover the 166 contribution to the six-vector class using exact CRT bounds."""

import argparse
from collections import Counter, defaultdict
from functools import cache
from fractions import Fraction
from math import ceil
from pathlib import Path
from tempfile import TemporaryDirectory
import fcntl
import os
import signal
import subprocess
from exact_projection import compile_solver, primes
from linear_shape166 import bound
from six_vector_nodes import fields, records, fixture, BASE
from exact_scalar import graph_factors, collect, multiply, MASK


@cache
def parts(order, f):
    data, _ = fields(0, f)
    return path_parts(order, data, fixture()[f])


def path_parts(order, data, weighted_polars):
    pol = [{}, {}]
    w = Counter()
    for weight, us in weighted_polars:
        ix = []
        for h, g in enumerate((order[:5], order[5:])):
            q = tuple(us[j] for j in g)
            ix.append(pol[h].setdefault(q, len(pol[h])))
        w[tuple(ix)] += weight
    lines = []
    for h, g in enumerate((order[:5], order[5:])):
        lines.extend((g, tuple(a for j in g for a in data[0][j]), (len(pol[h]),)))
        for qs in pol[h]:
            lines.append(tuple(t for q in qs for t in q))
    ends = []
    for flow, _, (left, right) in graph_factors(data, order, ()):
        z = defaultdict(dict)
        for (i, j, k), a in left.items():
            for (h, t, l), b in right.items():
                if h == i:
                    collect(z[5 * j + t, 5 * k + l], multiply(a, b), flow)
        rr = [
            (x, y, k & MASK, k >> 20, c)
            for (x, y), poly in z.items()
            for k, c in poly.items()
            if c
        ]
        ends.extend(((len(rr),), *rr))
    return (
        lines,
        [[w[i, j] for j in range(len(pol[1]))] for i in range(len(pol[0]))],
        ends,
    )


def run(solver, prime, limit=28):
    shape = (1, 6, 6)
    out = []
    for order, det in list(records(shape))[:limit]:
        for f in range(2):
            a, b, c = parts(order, f)
            lines = [
                (len(out), det * (-60480) % prime, 120),
                *a,
                *[[v % prime for v in row] for row in b],
                *c,
            ]
            out.append("".join(" ".join(map(str, line)) + "\n" for line in lines))
    r = subprocess.run(
        [solver],
        input="".join(out),
        capture_output=True,
        text=True,
        check=True,
        env={**os.environ, "VECLIB_MAXIMUM_THREADS": "1"},
    )
    rows = [tuple(map(int, x.split())) for x in r.stdout.splitlines()]
    assert [i for i, v in rows] == list(range(len(out)))
    return [v for i, v in rows]


def numerator_bound():
    # The physical denominator bound is the first row of theta_denominators.
    # bound() estimates every order and normalized integral before the shape factor.
    den = 293772917487329280000
    groups = [list(records((1, 6, 6)))]
    limit = Fraction(0)
    for f, samples in fixture().items():
        (mom, _), _ = fields(0, f)
        for w, polar in samples:
            normalized, integer = bound(groups, mom, polar)
            limit += abs(w) * Fraction(integer, normalized)
    return den, ceil(den * limit * Fraction(60480, 2 * 158400))


def main(args):
    den, maximum = numerator_bound()
    print(
        "common denominator",
        den,
        "numerator bit bound",
        maximum.bit_length(),
        flush=True,
    )
    if args.bound_only:
        return
    path = BASE / "six_vector_shapes" / "166_crt.tsv"
    with (
        path.open("a+") as file,
        TemporaryDirectory(prefix="twist-six-vector-crt-") as temp,
    ):
        fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        file.seek(0)
        saved = dict(tuple(map(int, l.split())) for l in file if not l.startswith("#"))
        if not saved and file.tell() == 0:
            print(
                "# six-vector 166 physical coefficient divided by pi: prime, residue",
                file=file,
                flush=True,
            )
        modulus = 1
        number = 0
        for i, p in enumerate(primes()):
            if i >= args.max_primes:
                return
            if p not in saved:
                solver = Path(temp) / "partition"
                compile_solver(solver, p)
                vals = run(solver, p)
                saved[p] = sum(vals) * pow(2 * 158400, -1, p) % p
                print(p, saved[p], sep="\t", file=file, flush=True)
                print("prime", p, "residue", saved[p], flush=True)
            number += modulus * ((saved[p] * den - number) * pow(modulus, -1, p) % p)
            modulus *= p
            if modulus > 2 * maximum:
                if number > modulus // 2:
                    number -= modulus
                assert abs(number) <= maximum and (number - 472 * den) % 1009 == 0
                print(
                    "coefficient divided by pi",
                    Fraction(number, den),
                    "primes",
                    i + 1,
                    flush=True,
                )
                return


if __name__ == "__main__":

    def interrupt(sig, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    ap = argparse.ArgumentParser()
    ap.add_argument("--bound-only", action="store_true")
    ap.add_argument("--max-primes", type=int, default=100)
    main(ap.parse_args())
