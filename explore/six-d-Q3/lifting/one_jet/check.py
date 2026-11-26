"""Exact sl2 relative CE matrices for one holomorphic derivative."""

from collections import Counter, defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quadratic import exact_rank
from sl2_zero_jets import BRACKET, basis, differential as base_differential

WEIGHTS = {"c": (0, 0), "x": (1, 0), "y": (0, 1), "b": (1, 1)}


def plus(vec, index):
    return tuple(val + (i == index) for i, val in enumerate(vec))


def derivative(mon):
    x, y, bs, kind, index = mon
    row = Counter()
    for target, coeff in base_differential((x, y, bs)).items():
        row[target + (kind, index)] += coeff
    sign = (-1) ** len(bs)
    for i, j, structure in BRACKET[index]:
        coeff = sign * structure
        if kind == "x":
            row[plus(x, j), y, bs, "c", i] += coeff
        elif kind == "y":
            row[x, plus(y, j), bs, "c", i] += coeff
        elif kind == "b":
            row[x, plus(y, j), bs, "x", i] += coeff
            row[plus(x, i), y, bs, "y", j] += coeff
            if j not in bs:
                inversions = sum(k > j for k in bs)
                new_bs = tuple(sorted(bs + (j,)))
                row[x, y, new_bs, "c", i] -= coeff * (-1) ** inversions
    return {m: n for m, n in row.items() if n}


def chains(p, q, charge):
    spaces = defaultdict(list)
    for kind, (dp, dq) in WEIGHTS.items():
        for index in range(3):
            for degree in range(4):
                for x, y, bs in basis(p-dp, q-dq, degree, charge-(1-index)):
                    length = 1 + sum(x) + sum(y) + degree
                    spaces[length].append((x, y, bs, kind, index))
    return spaces


def cohomology(p, q, charge):
    spaces = chains(p, q, charge)
    ranks = {}
    for length, monomials in spaces.items():
        rows = []
        for mon in monomials:
            row = derivative(mon)
            square = Counter()
            for target, coeff in row.items():
                for key, val in derivative(target).items():
                    square[key] += coeff * val
            assert not any(square.values())
            rows.append(row)
        ranks[length] = exact_rank(rows)
    return Counter({k: len(v) - ranks[k] - ranks.get(k-1, 0) for k, v in spaces.items()})


def predicted(p, q):
    result = Counter()
    for dp, dq in WEIGHTS.values():
        a, b = p-dp, q-dq
        if min(a, b) < 0:
            continue
        if (a+b) % 2 == 1:
            result[1+a+b] += 1
        if min(a, b) >= 1 and a+b >= 4 and (a+b) % 2 == 0:
            result[a+b] += 1
    return result


if __name__ == "__main__":
    total = 0
    for p in range(6):
        for q in range(6):
            singlets = cohomology(p, q, 0)
            singlets.subtract(cohomology(p, q, 1))
            assert singlets == predicted(p, q), (p, q, singlets, predicted(p, q))
            total += sum(singlets.values())
    print(f"36 bidegrees, {total} relative classes; exact ranks and all differential squares passed")
