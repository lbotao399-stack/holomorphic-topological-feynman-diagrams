"""Test the direct transverse fan reduction on n <= 2 from the paper's coordinates."""

from functools import cache
from math import comb, factorial as fac
from itertools import product
import sympy as sp
from sympy.polys.rings import ring
from fractions import Fraction as F


@cache
def coordinates(n):
    R, *r = ring("x,y,u,v", sp.ZZ)
    x, y, u, v = r
    if n == 1:
        return R, (1 - u, u), (v, 1 - v), R.one
    if n == 2:
        return (
            R,
            (1 - x + x * (1 - u) * (1 - y), x * (1 - u) * y, x * u),
            (v * (1 - y), v * y, 1 - v),
            x * y * v,
        )
    raise ValueError("the replay has coordinates for n=1,2")


@cache
def axis_coefficient(n, a, b, k):
    R, x, y, _ = coordinates(n)
    if min(a, b) < 0 or a + b != sum(k):
        return R.zero
    val = R.zero
    for s in product(*(range(m + 1) for m in k)):
        if sum(s) != a:
            continue
        z = R.one
        for xi, yi, ki, si in zip(x, y, k, s):
            z *= xi**si * yi ** (ki - si)
        val += z * (
            fac(a) * fac(b) // prod_factorials((*s, *(ki - si for ki, si in zip(k, s))))
        )
    return val


def prod_factorials(nums):
    result = 1
    for i in nums:
        result *= fac(i)
    return result


def graph(r, t, outputs):
    n = len(outputs) - 1
    R, x, y, measure = coordinates(n)
    value = R.zero
    for j in range(n + 1):
        contracted = (r[0] - n + j, r[1] - j, t[0] - j, t[1] - n + j)
        if min(contracted) < 0:
            continue
        factor = (
            (-1) ** j
            * comb(n, j)
            * prod_factorials((*r, *t))
            // prod_factorials(contracted)
        )
        pol = measure * axis_coefficient(
            n, contracted[0], contracted[2], tuple(k[0] for k in outputs)
        )
        pol *= axis_coefficient(
            n, contracted[1], contracted[3], tuple(k[1] for k in outputs)
        )
        value += pol * factor
    return -sum(
        F(int(c) * fac(u) * fac(v), (1 + x) * (1 + y) * fac(2 + u + v))
        for (x, y, u, v), c in value.items()
    )


from coefficients import ONE, add, sub, positive, times


def splits(upper):
    return product(range(upper[0] + 1), range(upper[1] + 1))


def choose2(upper, lower):
    return comb(upper[0], lower[0]) * comb(upper[1], lower[1])


@cache
def beta(a, b):
    return F(fac(a) * fac(b), fac(a + b + 1))


@cache
def triangle(a, b, c, d):
    # u^a (1-u)^b v^c (1-v)^d. Integrate v first.
    return sum(
        (-1) ** j * comb(d, j) * beta(a, b + c + j + 1) / (c + j + 1)
        for j in range(d + 1)
    )


def transferred(r, t, outputs):
    n = len(outputs) - 1
    if n < 1:
        raise ValueError("at least two output modes are required")
    if not positive(r, t, *outputs) or add(*outputs, times(n, ONE)) != add(r, t):
        return F()
    suffix = [(0, 0)] * (n + 2)
    for j in range(n, -1, -1):
        suffix[j] = add(outputs[j], suffix[j + 1])

    @cache
    def transfer(j, a, b):
        if not positive(sub(add(a, b), suffix[j])):
            return F()
        last = outputs[j]
        value = F()
        if j == n:
            for d in splits(a):
                e = sub(last, d)
                if positive(e, sub(b, e)):
                    value += (
                        choose2(a, d)
                        * choose2(b, e)
                        * triangle(
                            sum(d), sum(sub(a, d)), n - 1 + sum(sub(b, e)), sum(e)
                        )
                    )
            return value
        for d in splits(a):
            carried = sub(add(d, b), suffix[j + 1])
            pushed = sub(carried, last)
            if not positive(pushed):
                continue
            value += (
                choose2(a, d)
                * choose2(carried, last)
                * beta(n - j + sum(d), sum(sub(a, d)))
                * beta(1 + sum(last), j - 1 + sum(pushed))
                * transfer(j + 1, d, b)
            )
        return value

    value = F()
    for j in range(n + 1):
        a = sub(r, (n - j, j))
        b = sub(t, (j, n - j))
        if not positive(a, b):
            continue
        numerator = (-1) ** j * comb(n, j) * prod_factorials((*r, *t))
        value -= F(numerator, prod_factorials((*a, *b))) * transfer(1, a, b)
    return value
