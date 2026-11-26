"""S2 kernel of m parallel three-edge paths and the marked edge. No fitted constants."""

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
from math import factorial

ZERO = (0, 0, 0)


def step(power, axis, n=1):
    return tuple(x + (n if i == axis else 0) for i, x in enumerate(power))


def plus(a, b):
    return tuple(x + y for x, y in zip(a, b))


def crossed(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def combine(result, source, factor):
    if factor:
        for power, value in source.items():
            result[power] = result.get(power, 0) + value * factor
    return result


def clean(values):
    return {k: v for k, v in values.items() if v}


def determinant_product(pairs):
    values = {ZERO: 1}
    for p, q in pairs:
        result = defaultdict(int)
        for power, value in values.items():
            for i, c in enumerate(crossed(p, q)):
                result[step(power, i)] += value * c
        values = clean(result)
    return values


def divided_derivative(values, vector, order):
    result = defaultdict(int)
    for power, value in values.items():
        for i, v in enumerate(vector):
            if power[i] and v:
                result[step(power, i, -1)] += F(power[i] * v, order) * value
    return clean(result)


def transfer(values, p, q, local):
    # local[mu][a,b] multiplies w_i^|mu| L^mu. Its triangle moment is exact.
    result = {mu: {} for mu in local}
    first = values
    j = 0
    while first:
        second = first
        k = 0
        while second:
            for mu, entries in local.items():
                weight = sum(
                    F(c, (k + b + 1) * (j + k + a + b + 2))
                    for (a, b), c in entries.items()
                )
                combine(result[mu], second, weight * factorial(sum(mu)))
            k += 1
            second = divided_derivative(second, q, k)
        j += 1
        first = divided_derivative(first, p, j)
    return {mu: clean(value) for mu, value in result.items()}


def current_factors(p, q, insertions):
    terms = {(ZERO, 0, 0): 1}
    for edge, axis in insertions:
        constant = (0, p[axis], p[axis] + q[axis])[edge]
        result = defaultdict(int)
        for (mu, a, b), v in terms.items():
            result[step(mu, axis), a, b] += v
            result[mu, a, b] += constant * v
            result[mu, a + 1, b] -= p[axis] * v
            result[mu, a, b + 1] -= q[axis] * v
        terms = clean(result)
    pieces = defaultdict(dict)
    for (mu, a, b), c in terms.items():
        pieces[mu][a, b] = c
    return pieces


def kernel(source, pairs, powers=None):
    """Path orientation A->B. Ordered edge 0 is the mark, then three edges per path."""
    m = len(pairs)
    if not m:
        raise ValueError("at least one path required")
    insertions = [[] for _ in pairs]
    marked = ZERO
    bound = ZERO
    for (edge, axis), n in (powers or {}).items():
        if edge not in range(3 * m + 1) or axis not in range(3) or n < 0:
            raise ValueError("invalid current power")
        if edge == 0:
            marked = step(marked, axis, n)
        else:
            i, j = divmod(edge - 1, 3)
            insertions[i].extend([(j, axis)] * n)
            bound = step(bound, axis, n)
    branches = [current_factors(p, q, ins) for (p, q), ins in zip(pairs, insertions)]
    base = determinant_product(pairs)
    total = F()
    for target in product(*(range(c + 1) for c in bound)):
        initial = {plus(power, plus(target, marked)): v for power, v in base.items()}
        states = {ZERO: initial}
        for (p, q), entries in zip(pairs, branches):
            new = {}
            for used, poly in states.items():
                possible = {
                    mu: value
                    for mu, value in entries.items()
                    if all(x + y <= z for x, y, z in zip(used, mu, target))
                }
                for mu, pol in transfer(poly, p, q, possible).items():
                    combine(new.setdefault(plus(used, mu), {}), pol, 1)
            states = new
        last = states.get(target, {})
        value = sum(
            c * source[0] ** i * source[1] ** j * source[2] ** k
            for (i, j, k), c in last.items()
        )
        total += (
            F(factorial(2 + sum(marked)), factorial(m + 2 + sum(marked) + sum(target)))
            * value
        )
    return (-1) ** m * total
