"""A planar Hamiltonian three-cochain for the six-transverse-vector sector."""

from collections import Counter
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
from mixed_planar import basis, diff, lead, linear, subtract


def invariant_cycles(profile):
    # The two components of the chain boundary are CE and the raising operator.
    rows = {}
    cycles = []
    for k, weight in (("dx", -1), ("dy", 1)):
        for w in basis(5, weight)[3]:
            if tuple(i[1] for i in w) != profile:
                continue
            row = {(0, k, u): Fraction(v) for u, v in diff(w).items() if v}
            row.update({(1, k, u): Fraction(v) for u, v in linear(w, "raise").items()})
            if k == "dy":
                row[1, "dx", w] = Fraction(1)
            chain = {(k, w): Fraction(1)}
            while row:
                j = min(row)
                v = row[j]
                if j not in rows:
                    rows[j] = (
                        {i: a / v for i, a in row.items()},
                        {i: a / v for i, a in chain.items()},
                    )
                    break
                r, c = rows[j]
                subtract(row, r, v)
                subtract(chain, c, v)
            else:
                unit = lcm(*(x.denominator for x in chain.values()))
                scale = reduce(gcd, (int(x * unit) for x in chain.values()))
                if chain[min(chain)] < 0:
                    scale = -scale
                cycles.append({i: int(x * unit) // scale for i, x in chain.items()})
    return cycles


def stream():
    # Scalar cochains at Euler weight four, modulo closure and linear invariance.
    relations = {}
    cols = [diff(w) for w in basis(4, 0)[4]]
    cols += [linear(w, "raise") for w in basis(4, -2)[3]]
    for col in cols:
        row = {k: Fraction(v) for k, v in col.items() if v}
        k, a = lead(row, relations)
        if k is not None:
            relations[k] = {i: v / a for i, v in row.items()}

    def remainder(row):
        for k in sorted(relations):
            if row.get(k):
                subtract(row, relations[k], row[k])
        return row

    cycles = []
    images = []
    for profile in ((1, 4, 0), (1, 3, 1)):
        for cycle in invariant_cycles(profile):
            image = Counter()
            for (kind, w), weight in cycle.items():
                for k, v in linear(w, "dy" if kind == "dx" else "dx").items():
                    image[k] += (-1 if kind == "dx" else 1) * weight * v
            rem = remainder({k: Fraction(v) for k, v in image.items() if v})
            if rem:
                cycles.append(cycle)
                images.append(rem)
    free = min(images[0])
    vals = {w: remainder({w: Fraction(1)}).get(free, 0) for w in basis(4, 0)[3]}
    unit = lcm(*(x.denominator for x in vals.values()))
    coeff = {w: int(x * unit) for w, x in vals.items() if x}
    common = reduce(gcd, coeff.values())
    if coeff[min(coeff)] < 0:
        common = -common
    coeff = {w: x // common for w, x in coeff.items()}
    assert (len(cols), len(relations), len(vals), len(coeff)) == (97, 54, 56, 36)
    assert all(sum(coeff.get(w, 0) * c for w, c in col.items()) == 0 for col in cols)
    current = {}
    for kind, weight in (("dx", -1), ("dy", 1)):
        for w in basis(5, weight)[3]:
            coord = sum(
                coeff.get(k, 0) * v
                for k, v in linear(w, "dy" if kind == "dx" else "dx").items()
            ) * (-1 if kind == "dx" else 1)
            if coord:
                current[kind, w] = coord
    for w in basis(6, 0)[3]:
        assert (
            sum(
                current.get((k, u), 0) * v
                for k in ("dx", "dy")
                for u, v in linear(w, k).items()
            )
            == 0
        )
    for kind, weight in (("dx", -1), ("dy", 1)):
        for w in basis(5, weight)[4]:
            assert sum(current.get((kind, u), 0) * v for u, v in diff(w).items()) == 0
        for w in basis(5, weight - 2)[3]:
            assert (
                sum(
                    current.get((kind, u), 0) * v for u, v in linear(w, "raise").items()
                )
                + (current.get(("dx", w), 0) if kind == "dy" else 0)
                == 0
            )
    pairings = [sum(current.get(k, 0) * v for k, v in c.items()) for c in cycles]
    assert len(current) == 120 and [len(c) for c in cycles] == [24, 30]
    assert pairings == [-2520, 2520]
    fixture = ("dx", ((0, 1, 3), (0, 4, 0), (1, 0, 2)))
    assert current[fixture] == 360 and pairings[0] == -7 * current[fixture]
    print("Hamiltonian scalar and current supports", len(coeff), len(current))
    print("cycle pairings", pairings, "one planar fixture", current[fixture])
    return coeff, current, cycles


if __name__ == "__main__":
    stream()
