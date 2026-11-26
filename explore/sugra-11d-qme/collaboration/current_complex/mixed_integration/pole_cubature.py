"""Seventh-adic factors on the bounded joint degree-eight part of theta_166."""

from fractions import Fraction
from functools import cache
from collections import defaultdict
from math import factorial, prod
from cubature import compositions
from pole_weights import value


def residue(w):
    return w.numerator * pow(w.denominator, -1, 7) % 7


@cache
def simplex(s):
    nodes = defaultdict(int)
    for i in range(s + 1):
        r = s - i
        denominator = 6 + 2 * r
        weight = residue(
            Fraction(
                7 * (-1) ** i * factorial(5) * denominator ** (2 * s + 1),
                4**s * factorial(i) * factorial(6 + 2 * s - i),
            )
        )
        for a in compositions(r, 6):
            bary = tuple((1 + 2 * k) * pow(denominator, -1, 7) % 7 for k in a)
            nodes[bary] = (nodes[bary] + weight) % 7
    return {node: w for node, w in nodes.items() if w}


def joint():
    nodes = defaultdict(int)
    for layer, sgn in [(4, 1), (3, -1)]:
        for r in range(1, layer):
            for a, wa in simplex(r).items():
                for b, wb in simplex(layer - r).items():
                    nodes[a, b] = (nodes[a, b] + sgn * wa * wb) % 7
    return {node: w for node, w in nodes.items() if w}


def conductance():
    w = {
        (b, c): pow(b * c, 3, 7) * (1 - b - c) % 7
        for b in range(1, 7)
        for c in range(1, 7)
    }
    w = {nodes: v for nodes, v in w.items() if v}
    for d in range(11):
        for i in range(d + 1):
            assert sum(
                v * pow(b, i, 7) * pow(c, d - i, 7) for (b, c), v in w.items()
            ) % 7 == residue(7 * value(i, d - i))
    return w


def paths_check():
    nodes = joint()
    # All barycentric exponents up to separate degrees 7 and joint degree 8.
    # Orbit representatives under the independent S6 actions suffice.
    cases = [
        [(d, a) for a in compositions(d, 6) if tuple(sorted(a)) == a] for d in range(8)
    ]

    def one(x, e):
        return prod(pow(v, j, 7) for v, j in zip(x, e)) % 7

    @cache
    def exact(e):
        return residue(
            Fraction(7 * factorial(5) * prod(map(factorial, e)), factorial(5 + sum(e)))
        )

    for p in range(8):
        for q in range(min(7, 8 - p) + 1):
            for _, a in cases[p]:
                for _, b in cases[q]:
                    result = (
                        sum(v * one(x, a) * one(y, b) for (x, y), v in nodes.items())
                        % 7
                    )
                    assert result == exact(a) * exact(b) % 7
    return len(nodes)


if __name__ == "__main__":
    import sys

    part = sys.argv[1] if len(sys.argv) > 1 else "check"
    if part == "check":
        print(
            "bounded resistance nodes",
            paths_check(),
            "conductance nodes",
            len(conductance()),
        )
    elif part == "resistance":
        for (a, b), w in joint().items():
            print(*a, *b, w)
    elif part == "conductance":
        for (x, y), w in conductance().items():
            print(x, y, w)
    else:
        raise ValueError(part)
