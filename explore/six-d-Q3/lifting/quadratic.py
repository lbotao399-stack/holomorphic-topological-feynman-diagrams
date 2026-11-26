"""Exact local CE degree-two dimensions and odd-root contractions on C^{3|2}."""

from fractions import Fraction
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "s2code"))
from ops import L, Theory
from salg import EPS, TH1, TH2, add, gen, lderiv, mul, scale, sderiv


def monomials(weight):
    return product(*(range(min(1, d) + 1) if i >= 3 else range(d + 1)
                     for i, d in enumerate(weight)))


def subtract(a, b):
    return tuple(i - j for i, j in zip(a, b))


def parity(a):
    return (a[3] + a[4]) % 2


def multiply(a, b):
    weight = tuple(i + j for i, j in zip(a, b))
    if max(weight[3:]) > 1:
        return 0, weight
    return (-1) ** (a[4] * b[3]), weight


def exact_rank(rows):
    pivots = {}
    for row in rows:
        row = {key: Fraction(value) for key, value in row.items() if value}
        while row:
            lead = min(row)
            coeff = row[lead]
            if lead not in pivots:
                pivots[lead] = {key: value / coeff for key, value in row.items()}
                break
            for key, value in pivots[lead].items():
                updated = row.get(key, 0) - coeff * value
                if updated:
                    row[key] = updated
                else:
                    row.pop(key, None)
    return len(pivots)


def ce_dimension(weight):
    pairs = [(a, subtract(weight, a)) for a in monomials(weight)
             if max(subtract(weight, a)[3:]) <= 1]
    index = {pair: i for i, pair in enumerate(pairs)}
    rows = []
    for a, b in pairs:
        row = {index[a, b]: 1}
        j = index[b, a]
        row[j] = row.get(j, 0) + (-1) ** (parity(a) * parity(b))
        rows.append(row)
    for a in monomials(weight):
        rest = subtract(weight, a)
        for b in monomials(rest):
            c = subtract(rest, b)
            if max(c[3:]) > 1:
                continue
            row = {}
            for x, y, z in ((a, b, c), (b, c, a), (c, a, b)):
                sign, xy = multiply(x, y)
                if sign:
                    j = index[xy, z]
                    row[j] = row.get(j, 0) + sign * (-1) ** (parity(x) * parity(z))
            rows.append(row)
    return len(pairs) - exact_rank(rows)


def differential_dimension(weight):
    def exists(degree):
        return min(degree) >= 0 and max(degree[3:]) <= 1

    forms = sum(exists(tuple(x - (j == i) for j, x in enumerate(weight)))
                for i in range(5))
    return forms - int(any(weight) and exists(weight))


def check_roots():
    # The root uses native component functional derivatives at theta; the
    # inserted C is at two independent exterior generators eta.
    eta = [(0, 3, 0, 0, 0), (0, 4, 0, 0, 0)]
    eta1, eta2 = map(gen, eta)
    C = add(add(L("c", 0), mul(eta1, L("g1", 0))),
            add(mul(eta2, L("g2", 0)), mul(mul(eta1, eta2), L("b", 0))))
    delta = {}
    for t in (1, 2, 3, 4):
        delta = add(delta, sderiv(C, (t, 0, 0, 0, 0)))
    left = add(gen(TH1), eta1, -1)
    right = add(gen(TH2), eta2, -1)
    assert delta == mul(left, right)
    assert not mul(delta, delta)
    for odd in eta:
        derivative = lderiv(delta, odd)
        assert derivative and not mul(delta, derivative)
    assert mul(delta, lderiv(lderiv(delta, eta[0]), eta[1]))

    theory = Theory(3, EPS)
    for field in ("g1", "g2"):
        quadratic = {}
        for j in range(3):
            quadratic = add(quadratic, mul(L(field, j), L("b", j)))
        assert quadratic and not theory.Q0(quadratic) and not theory.Q1(quadratic)


if __name__ == "__main__":
    check_roots()
    dimensions = []
    for alpha in product(range(4), repeat=3):
        if sum(alpha) > 3:
            continue
        for beta in product(range(3), repeat=2):
            weight = alpha + beta
            dim = ce_dimension(weight)
            assert dim == differential_dimension(weight), (weight, dim)
            dimensions.append(dim)
    print(f"{len(dimensions)} multidegrees, {sum(dimensions)} independent quadratic classes")
    print("native odd delta, both first derivatives, two mixed scalar-antifield cocycles passed")
