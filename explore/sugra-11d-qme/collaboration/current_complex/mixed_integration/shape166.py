"""Sum every nonzero unordered partition of the scalar fixture with checkpointed exact rows."""

import argparse
from contextlib import redirect_stdout
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryFile
from partition166 import determinant, export, fixture


def solve_partition(solver, data, order, limit=120):
    with TemporaryFile(mode="w+") as source:
        with redirect_stdout(source):
            export(data, order, limit=limit)
        source.seek(0)
        result = subprocess.run(
            [solver], stdin=source, capture_output=True, text=True, check=True
        )
    match = re.fullmatch(
        r"orders (\d+), signed Wick ([\d/-]+), determinant-weighted ([\d/-]+)\n",
        result.stdout,
    )
    if not match or int(match[1]) != limit**2:
        raise ValueError(result.stdout)
    raw, weighted = Fraction(match[2]), Fraction(match[3])
    if weighted != determinant(data[0], order) * raw:
        raise ValueError("determinant mismatch")
    return weighted


def compute(solver, file):
    data = fixture()
    existing = {}
    if file.exists():
        for line in file.read_text().splitlines():
            order, value = line.split()
            if order in existing:
                raise ValueError("duplicate partition")
            existing[order] = Fraction(value)
    total = Fraction()
    seen = set()
    for b in combinations(range(9), 5):
        order = b + tuple(x for x in range(10) if x not in b)
        if not determinant(data[0], order):
            continue
        word = "".join(map(str, order))
        seen.add(word)
        if word not in existing:
            existing[word] = solve_partition(solver, data, order)
            with file.open("a") as out:
                print(word, existing[word], file=out)
        total += existing[word]
    if seen != set(existing):
        raise ValueError("unexpected checkpoint row")
    print(len(seen), "unordered partitions, determinant-weighted sum", total)
    print(
        "after path reversal and scalar kernel, coefficient of pi", total / (2 * 158400)
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("solver")
    parser.add_argument("checkpoint", type=Path)
    args = parser.parse_args()
    compute(args.solver, args.checkpoint)
