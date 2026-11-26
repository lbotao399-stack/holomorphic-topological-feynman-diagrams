"""Sparse ranks over F_p and over Q (exact).  Rows are dicts {column_key: int or Fraction}.

rank_p(rows, p)            rank over F_p
rank_Q(rows)               rank over Q (fraction-free integer elimination, exact)
SpanP / SpanQ              incremental row spaces: add(row) -> bool, dim()
"""
from fractions import Fraction

P1 = 2147483647          # 2^31 - 1
P2 = 2147483629


class SpanP:
    """Row space over F_p.  piv[c] = reduced row whose minimal column is c (leading coeff 1)."""

    def __init__(self, p=P1):
        self.p = p
        self.piv = {}

    def reduce(self, row):
        p = self.p
        r = {}
        for k, v in row.items():
            v %= p
            if v:
                r[k] = v
        piv = self.piv
        while r:
            c = min(r)
            pr = piv.get(c)
            if pr is None:
                return r
            a = r[c]
            for k, v in pr.items():
                w = (r.get(k, 0) - a * v) % p
                if w:
                    r[k] = w
                else:
                    r.pop(k, None)
        return r

    def add(self, row):
        r = self.reduce(row)
        if not r:
            return False
        c = min(r)
        inv = pow(r[c], self.p - 2, self.p)
        self.piv[c] = {k: (v * inv) % self.p for k, v in r.items()}
        return True

    def dim(self):
        return len(self.piv)


class SpanQ:
    """Row space over Q with integer rows (primitive, leading coefficient positive)."""

    def __init__(self):
        self.piv = {}

    @staticmethod
    def _prim(r):
        from math import gcd
        g = 0
        for v in r.values():
            g = gcd(g, v)
            if g == 1:
                break
        if g > 1:
            r = {k: v // g for k, v in r.items()}
        c = min(r)
        if r[c] < 0:
            r = {k: -v for k, v in r.items()}
        return r

    def reduce(self, row):
        r = {}
        for k, v in row.items():
            if isinstance(v, Fraction):
                if v.denominator != 1:
                    raise ValueError("SpanQ expects integer rows")
                v = v.numerator
            if v:
                r[k] = v
        piv = self.piv
        while r:
            c = min(r)
            pr = piv.get(c)
            if pr is None:
                return self._prim(r)
            a = r[c]
            b = pr[c]
            # r <- b*r - a*pr  (clears column c)
            new = {}
            for k, v in r.items():
                new[k] = b * v
            for k, v in pr.items():
                w = new.get(k, 0) - a * v
                if w:
                    new[k] = w
                else:
                    new.pop(k, None)
            r = self._prim(new) if new else new
        return r

    def add(self, row):
        r = self.reduce(row)
        if not r:
            return False
        self.piv[min(r)] = r
        return True

    def dim(self):
        return len(self.piv)


def rank_p(rows, p=P1):
    S = SpanP(p)
    for r in rows:
        S.add(r)
    return S.dim()


def rank_Q(rows):
    S = SpanQ()
    for r in rows:
        S.add(r)
    return S.dim()


def make_span(field):
    """field: 'Q' or a prime p."""
    if field == 'Q':
        return SpanQ()
    return SpanP(field)
