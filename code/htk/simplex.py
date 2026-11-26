"""Exact Dirichlet integrals over products of simplices.

    int_{Delta^m} prod_{j=0..m} b_j^{a_j} db_1 ... db_m = prod_j Gamma(a_j + 1) / Gamma(m + 1 + sum_j a_j),
    Delta^m = {b_j > 0, sum_j b_j = 1},  a_j > -1 integer or half-integer.

With h half-integer exponents the value is a rational number times pi^{floor(h/2)}.
A number sum_k q_k pi^k is represented by PiNumber({k: q_k}).
"""
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial

__all__ = ["gamma_half", "dirichlet", "PiNumber"]


@lru_cache(maxsize=None)
def gamma_half(x):
    """Gamma(x) for x > 0 integer or half-integer, as (q, h): Gamma(x) = q * pi^(h/2), h in {0, 1}."""
    x = Fr(x)
    if x <= 0:
        raise ValueError("Gamma at %s: divergent Dirichlet integral" % x)
    if x.denominator == 1:
        return Fr(factorial(int(x) - 1)), 0
    if x.denominator != 2:
        raise ValueError("only integer and half-integer exponents are supported")
    m = (x.numerator - 1) // 2
    return Fr(factorial(2 * m), 4 ** m * factorial(m)), 1


@lru_cache(maxsize=None)
def dirichlet(exps):
    """exps: tuple (a_0, ..., a_m) of Fractions.  Returns (q, k) with integral = q * pi^k."""
    q, h = Fr(1), 0
    tot = Fr(0)
    for a in exps:
        g, hh = gamma_half(Fr(a) + 1)
        q *= g
        h += hh
        tot += Fr(a) + 1
    g, hh = gamma_half(tot)
    q /= g
    h -= hh
    assert h % 2 == 0
    return q, h // 2


class PiNumber:
    """sum_k q_k pi^k with rational q_k."""
    __slots__ = ("c",)

    def __init__(self, c=None):
        if c is None:
            c = {}
        elif not isinstance(c, dict):
            c = {0: Fr(c)}
        self.c = {int(k): Fr(v) for k, v in c.items() if v != 0}

    def is_rational(self):
        return all(k == 0 for k in self.c)

    def rational(self):
        if not self.is_rational():
            raise ValueError("%r is not rational" % self)
        return self.c.get(0, Fr(0))

    def coeff(self, k):
        return self.c.get(k, Fr(0))

    def __add__(self, o):
        o = o if isinstance(o, PiNumber) else PiNumber(o)
        out = dict(self.c)
        for k, v in o.c.items():
            out[k] = out.get(k, 0) + v
        return PiNumber(out)

    __radd__ = __add__

    def __neg__(self):
        return PiNumber({k: -v for k, v in self.c.items()})

    def __sub__(self, o):
        o = o if isinstance(o, PiNumber) else PiNumber(o)
        return self + (-o)

    def __rsub__(self, o):
        return PiNumber(o) - self

    def __mul__(self, o):
        if isinstance(o, PiNumber):
            out = {}
            for k, v in self.c.items():
                for l, w in o.c.items():
                    out[k + l] = out.get(k + l, 0) + v * w
            return PiNumber(out)
        return PiNumber({k: v * Fr(o) for k, v in self.c.items()})

    __rmul__ = __mul__

    def __truediv__(self, o):
        return PiNumber({k: v / Fr(o) for k, v in self.c.items()})

    def __eq__(self, o):
        try:
            o = o if isinstance(o, PiNumber) else PiNumber(o)
        except (TypeError, ValueError):
            return NotImplemented
        return self.c == o.c

    def __hash__(self):
        return hash(tuple(sorted(self.c.items())))

    def __float__(self):
        from math import pi
        return float(sum(float(v) * pi ** k for k, v in self.c.items()))

    def sympy(self):
        import sympy as sp
        return sum((sp.Rational(v.numerator, v.denominator) * sp.pi ** k for k, v in self.c.items()),
                   sp.Integer(0))

    def __repr__(self):
        if not self.c:
            return "0"
        parts = []
        for k in sorted(self.c):
            v = self.c[k]
            parts.append(str(v) if k == 0 else "%s*pi" % v if k == 1 else "%s*pi^%d" % (v, k))
        return " + ".join(parts)
