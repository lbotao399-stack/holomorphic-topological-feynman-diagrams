"""Lift the one-mark transverse current through the planar divergence."""

from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from functools import reduce

from mixed_planar import basis, diff, graph_test, lead, linear


def conditions(extra_relations=None):
    # dx carries the first planar covector in the dual-divergence chain and
    # dy the second. The raising operator takes dy to dx with coefficient one.
    relations = {}
    columns = []

    def add(row):
        row = {k: Fraction(v) for k, v in row.items() if v}
        columns.append(row.copy())
        key, val = lead(row, relations)
        if key is not None:
            relations[key] = {k: v / val for k, v in row.items()}

    words = []
    for kind, weight in (("dx", -1), ("dy", 1)):
        states = basis(5, weight)
        words.extend((kind, w) for w in states[4])
        for w in states[5]:
            add({(kind, k): v for k, v in diff(w).items()})
        for w in basis(5, weight - 2)[4]:
            row = {(kind, k): v for k, v in linear(w, "raise").items()}
            if kind == "dy":
                row["dx", w] = 1
            add(row)

    assert (len(words), len(columns), len(relations)) == (224, 289, 215)
    if extra_relations is not None:
        for row in extra_relations(words):
            add(row)

    target = Counter()
    for w, c in graph_test().items():
        for kind in ("dx", "dy"):
            for k, v in linear(w, kind).items():
                target[kind, k] += c * v

    # Reduce all pivots, including those following a free coordinate.
    def remainder(row):
        row = Counter(row)
        for k in sorted(relations):
            coeff = row[k]
            if coeff:
                for j, v in relations[k].items():
                    row[j] -= coeff * v
        return {k: v for k, v in row.items() if v}

    return words, columns, relations, target, remainder


def lift(extra_relations=None):
    words, columns, relations, target, remainder = conditions(extra_relations)
    rem = remainder(target)
    if not rem:
        raise ArithmeticError("dual divergence is a boundary")
    free = min(rem)
    current = {w: remainder({w: 1}).get(free, 0) for w in words}
    unit = lcm(*(c.denominator for c in current.values()))
    current = {k: int(v * unit) for k, v in current.items() if v}
    common = reduce(gcd, current.values())
    current = {k: v // common for k, v in current.items()}
    if current[min(current)] < 0:
        current = {k: -v for k, v in current.items()}
    for col in columns:
        assert sum(current.get(k, 0) * v for k, v in col.items()) == 0

    scalar = {}
    for w in basis(6, 0)[4]:
        v = sum(
            current.get((kind, k), 0) * n
            for kind in ("dx", "dy")
            for k, n in linear(w, kind).items()
        )
        if v:
            scalar[w] = v
    for w in basis(6, 0)[5]:
        assert sum(v * scalar.get(k, 0) for k, v in diff(w).items()) == 0
    for w in basis(6, -2)[4]:
        assert sum(v * scalar.get(k, 0) for k, v in linear(w, "raise").items()) == 0
    pairing = sum(v * scalar.get(w, 0) for w, v in graph_test().items())
    assert pairing == sum(v * current.get(k, 0) for k, v in target.items())
    if extra_relations is None:
        assert (len(current), len(scalar), pairing) == (114, 148, 2365440)
    return current, scalar, pairing, len(words), len(relations), len(columns)


if __name__ == "__main__":
    current, scalar, pairing, count, rank, checks = lift()
    print("planar coordinates", count, "relations", checks, "rank", rank)
    print("current and divergence nonzero coordinates", len(current), len(scalar))
    print("divergence on the integer four-cycle", pairing)
