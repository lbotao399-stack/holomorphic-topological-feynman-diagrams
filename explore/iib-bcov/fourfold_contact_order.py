"""Exact contact order of the primary fourfold vertex on conserved polynomial arcs."""

from functools import cache
from itertools import combinations
from math import inf
from random import Random

from flint import fmpq, fmpq_poly

from worldsheet_primary import fourfold_primary


def valuation(polynomial):
    return next((i for i, coefficient in enumerate(polynomial) if coefficient), inf)


def minor_orders(points):
    """Rows are momenta. Entries are scalars or ascending polynomial coefficients."""
    if len(points) < 5 or any(len(p) != 4 for p in points):
        raise ValueError("Expected at least five four-component momenta")

    def polynomial(entry):
        if isinstance(entry, fmpq_poly):
            return entry
        if isinstance(entry, (tuple, list)):
            return fmpq_poly([fmpq(str(x)) for x in entry])
        return fmpq_poly([fmpq(str(entry))])

    rows = [tuple(map(polynomial, row)) for row in points]
    if any(map(sum, zip(*rows))):
        raise ValueError("Momentum must be conserved identically")

    @cache
    def determinant(slots, columns):
        if not slots:
            return fmpq_poly([1])
        return sum(((-1)**i * rows[slots[0]][j]
                    * determinant(slots[1:], columns[:i] + columns[i+1:])
                    for i, j in enumerate(columns) if rows[slots[0]][j]), fmpq_poly())

    return {s: valuation(determinant(s, (0, 1, 2, 3)))
            for s in combinations(range(len(points)), 4)}


def fourfold_order(points):
    orders = minor_orders(points)
    minimum = min(orders.values())
    if minimum == inf:
        return inf
    per_slot = [min(v for slots, v in orders.items() if i in slots)
                for i in range(len(points))]
    genus = len(points) - 3
    return 2*genus*minimum + 2*sum(v-minimum for v in per_slot)


def covering_order(orders, count):
    """Independent minimum over genus many four-row bases jointly covering all slots."""
    bases = [(sum(1 << i for i in slots), order) for slots, order in orders.items()
             if order < inf]
    best = {0: 0}
    for _ in range(count-3):
        following = {}
        for mask, order in bases:
            for covered, cost in best.items():
                union = mask | covered
                following[union] = min(following.get(union, inf), cost+order)
        best = following
    return 2*best.get((1 << count)-1, inf)


def interpolated_order(genus, points, bound):
    """Interpolate both exact graph period channels at integers up to a proven degree bound."""
    differences = []
    for x in range(bound+1):
        momenta = [[str(a(x)) for a in row] for row in points]
        differences.append(fourfold_primary(genus, momenta))
    result = [fmpq_poly(), fmpq_poly()]
    binomial = fmpq_poly([1])
    for j in range(bound+1):
        for i, coefficient in enumerate(differences[0]):
            result[i] += fmpq(str(coefficient)) * binomial
        differences = [(b[0]-a[0], b[1]-a[1]) for a, b in zip(differences, differences[1:])]
        binomial *= fmpq_poly([-j, 1]) / (j+1)
    return min(map(valuation, result))


def check():
    rng = Random(24971)
    t = fmpq_poly([0, 1])
    for genus in range(2, 8):
        # Normalize a coordinate basis, add genus-2 soft waves, then drop the
        # leading target rank by two independent powers of t.
        raw = [tuple(fmpq_poly([int(i == j)]) for j in range(4)) for i in range(4)]
        for j in range(genus-2):
            raw.append(tuple(t**(1+j%3)*fmpq_poly([rng.randrange(1, 4), rng.randrange(-1, 2)])
                             for _ in range(4)))
        raw.append(tuple(-sum(col) for col in zip(*raw)))
        scaled = [tuple(a*factor for a, factor in zip(row, (1, 1, t, t*t))) for row in raw]
        expected = 6*genus + 2*sum(1+j%3 for j in range(genus-2))
        assert fourfold_order(scaled) == expected == covering_order(minor_orders(scaled), len(raw))
        if genus <= 3:
            # Coordinate factors add exactly 6g to the degree. Only one
            # extra input at genus 3, of degree <=2, can enter each minor.
            bound = 6*genus + (4*genus if genus == 3 else 0)
            assert interpolated_order(genus, scaled, bound) == expected
        # Neither deficient rank nor an identically zero slot admits a covering term.
        deficient = [(*row[:3], 0) for row in raw]
        assert fourfold_order(deficient) == inf
        zero_slot = [*raw[:-2], tuple(a+b for a,b in zip(*raw[-2:])), (0,0,0,0)]
        assert fourfold_order(zero_slot) == inf
        print("genus", genus, "fourfold order", expected)


if __name__ == "__main__":
    check()
