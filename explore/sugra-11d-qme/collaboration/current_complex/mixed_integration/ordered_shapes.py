"""Exact three-path Schwinger integral for a labeled graph with one one-form."""

from collections import Counter
from fractions import Fraction
from functools import cache
from math import prod
from pathlib import Path
import sys


from exact_scalar import (
    MASK,
    collect,
    components,
    fixture,
    graph_factors,
    multiply,
    resistance,
    splits,
)
from theta_denominators import SHAPES, angular


def paths(shape, order):
    start = 0
    out = []
    for size in shape:
        out.append(tuple(order[start : start + size - 1]))
        start += size - 1
    return out


def numerator(data, shape, order, lie):
    result = {}
    for sign, chain, (left, right) in graph_factors(
        data, order, range(3), True, lie, shape
    ):
        # Contract one path matrix at a time into the left endpoint tensor.
        tensor = left
        for slot, matrix in chain.items():
            new = {}
            for key, coeff in tensor.items():
                for j, poly in enumerate(matrix[key[slot]]):
                    if poly:
                        idx = key[:slot] + (j,) + key[slot + 1 :]
                        collect(new.setdefault(idx, {}), multiply(coeff, poly))
            tensor = {k: {e: v for e, v in p.items() if v} for k, p in new.items()}
        for key, p in tensor.items():
            if key in right:
                collect(result, multiply(p, right[key]), sign)
    return {k: v for k, v in result.items() if v}


@cache
def ternary_splits(key):
    degree = sum(components(key))
    return [
        (a, b, key - a - b, (da, db, degree - da - db), ca * cb)
        for a, da, ca in splits(key)
        for b, db, cb in splits(key - a)
    ]


def expectation(poly, shape, order, momenta):
    waves = paths(shape, order)
    means = [resistance(tuple(momenta[i] for i in path)) for path in waves]
    # Each angular average depends only on the numbers of copies assigned
    # to each path, independent of their component labels.
    grouped = Counter()
    for key, v in poly.items():
        for ba, bb, bc, bcounts, u in ternary_splits(key & MASK):
            for ca, cb, cc, ccounts, t in ternary_splits(key >> 20):
                exps = (ba + ca, bb + cb, bc + cc)
                if any(h and not path for h, path in zip(exps, waves)):
                    continue
                grouped[exps, bcounts, ccounts] += v * u * t
    total = Fraction()
    for (exps, b, c), v in grouped.items():
        if v:
            radial = prod(k(a) for k, a in zip(means, exps))
            total += v * radial * angular(shape, b, c)
    return total


def determinant(data, shape, order):
    # The signed planar path order is already included by keeping columns in
    # external label order. Rows are (A+B, B+C) in the paper.
    current = str(Path(__file__).resolve().parents[1])
    if current not in sys.path:
        sys.path.insert(0, current)
    from mixed_graph_kernel import det, lift

    wave = paths(shape, order)
    return det(
        tuple(
            lift(data[i], next(h for h, row in enumerate(wave) if i in row))
            for i in range(10)
        )
    )


def run():
    data = fixture()
    labels = ("0146723589", "0146723589", "0146723589", "6901472358", "0691472358")
    for (shape, _, _), seq in zip(SHAPES, labels):
        order = tuple(map(int, seq))
        signed = determinant(data[0], shape, order) * (-1) ** (9 - order.index(9))
        poly = numerator(data, shape, order, True)
        answer = expectation(poly, shape, order, data[0])
        print(
            "shape",
            shape,
            "terms",
            len(poly),
            "signed determinant",
            signed,
            "normalized Wick",
            answer,
            "weighted",
            signed * answer,
            flush=True,
        )
        if shape == (1, 6, 6):
            from exact_scalar import polynomials, expectation as independent

            assert poly == polynomials(data, order, lie=True)
            assert (
                answer
                == independent(poly, data[0], order)
                == Fraction(4170228116417, 25758708393600)
            )


if __name__ == "__main__":
    run()
