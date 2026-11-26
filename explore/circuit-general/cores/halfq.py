"""Helpers for half-integer q (odd dimension 2q, T = 1 kernels).

Values are dicts  monomial -> Fraction,  monomial = sorted tuple of (symbol, power) with
symbols 'sqrtpi', 'L2', 'z2', 'z3', ...  ('sqrtpi', 2k) = pi^k.
"""
from fractions import Fraction as F
from math import factorial

import prop


def gamma(x):
    """Gamma(x) for x a positive integer or half-integer: (Fraction, sqrtpi power)"""
    x = F(x)
    if x.denominator == 1:
        if x <= 0:
            raise ZeroDivisionError('Gamma pole')
        return F(factorial(int(x) - 1)), 0
    m = int(x - F(1, 2))
    if m >= 0:
        return F(factorial(2 * m), 4 ** m * factorial(m)), 1
    # Gamma(1/2 - k) = (-4)^k k! sqrt(pi) / (2k)!
    k = -m
    return F((-4) ** k * factorial(k), factorial(2 * k)), 1


def G2q(a, b, q):
    """Gamma(q-a)Gamma(q-b)Gamma(a+b-q)/(Gamma(a)Gamma(b)Gamma(2q-a-b)) as (Fraction, sqrtpi power)"""
    num = [gamma(q - a), gamma(q - b), gamma(a + b - q)]
    den = [gamma(a), gamma(b), gamma(2 * q - a - b)]
    r = F(1)
    k = 0
    for v, p in num:
        r *= v
        k += p
    for v, p in den:
        r /= v
        k -= p
    return r, k


def mono_mul(m, sym_pow):
    d = dict(m)
    for s, p in sym_pow:
        d[s] = d.get(s, 0) + p
    return tuple(sorted((s, p) for s, p in d.items() if p))


def scale_value(val, r, k):
    out = {}
    for m, v in val.items():
        mm = mono_mul(m, (('sqrtpi', k),)) if k else m
        out[mm] = out.get(mm, 0) + v * r
    return {m: v for m, v in out.items() if v}


def normalize(val):
    """zeta(2)^k -> pi^(2k)/6^k ; returns dict monomial -> Fraction"""
    out = {}
    for m, v in val.items():
        d = dict(m)
        k = d.pop('z2', 0)
        if k:
            d['sqrtpi'] = d.get('sqrtpi', 0) + 4 * k
            v = v / F(6) ** k
        mm = tuple(sorted((s, p) for s, p in d.items() if p))
        out[mm] = out.get(mm, 0) + v
    return {m: v for m, v in out.items() if v}


def wheel_value(spoke, rim, q):
    """projective W_n period at half-integer q by the momentum-space reducer (spoke cut)"""
    n = len(spoke)
    edges = [(0, i + 1) for i in range(n)] + [(i + 1, (i + 1) % n + 1) for i in range(n)]
    powers = list(spoke) + list(rim)
    nonint = [i for i, p in enumerate(powers) if F(p).denominator != 1]
    if len(nonint) > 1:
        raise ValueError('more than one half-integer power')
    cut = nonint[0] if nonint else 0
    return normalize(prop.cut_period(edges, powers, q, cut=cut))
