"""Scalar f coefficients in units of (16 pi^2)^(-m). Equations (74),(80)."""

from fractions import Fraction as F
from functools import cache
from math import factorial as fac


def add(*xs):
    return tuple(map(sum, zip(*xs)))


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def positive(*xs):
    return all(min(x) >= 0 for x in xs)


def times(n, p):
    return (n * p[0], n * p[1])


ONE = (1, 1)


def moment(r, k, l):
    if not positive(r, k, l):
        return 0
    return F(
        fac(r[0]) * fac(r[1]) * fac(1 + sum(k)) * fac(sum(l)),
        fac(1 + sum(r)) * fac(k[0]) * fac(k[1]) * fac(l[0]) * fac(l[1]),
    )


@cache
def one(r, t, k, l):
    if not positive(r, t, k, l) or add(k, l, ONE) != add(r, t):
        return F()
    result = F()
    for j in range(1, t[0] + 1):
        rest = (t[0] - j, t[1])
        tail = (l[0] + 1 - j, l[1])
        result += moment(add(r, rest), k, tail) - moment(rest, sub(k, r), tail)
    for j in range(1, t[1] + 1):
        result -= moment(add(r, (0, t[1] - j)), k, (l[0] - t[0], l[1] + 1 - j))
    return result


@cache
def coefficient(r, t, outputs, reverse_steps=False):
    m = len(outputs) - 1
    if m < 1:
        raise ValueError("at least two output modes are required")
    if (
        not positive(r, t, *outputs)
        or min(sum(r), sum(t)) < m
        or add(*outputs, times(m, ONE)) != add(r, t)
    ):
        return F()
    if m == 1:
        return one(r, t, *outputs)
    axis = int(not t[0]) if not reverse_steps else int(bool(t[1]))
    step = ((1, 0), (0, 1))[axis]
    reduced = sub(t, step)
    last = outputs[-1]
    l = sub(add(outputs[-2], last, ONE), step)
    h = sub(sub(add(reduced, step), ONE), last)
    return (
        coefficient(r, reduced, (*outputs[:-1], sub(last, step)), reverse_steps)
        - coefficient(r, reduced, (*outputs[:-2], l), reverse_steps)
        * one(l, step, outputs[-2], last)
        + coefficient(r, h, outputs[:-1], reverse_steps) * one(reduced, step, h, last)
    )
