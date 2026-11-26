"""wick.py -- exact forward-mode derivatives (Dual) and form-valued Gaussian polynomials (Poly).

Dual(value, grad): grad is a dict {coordinate index: d value / d coordinate}; d(Dual) = one-form.
Poly: dict {monomial (sorted tuple of Gaussian variable indices): form}; E(P) = sum_m <x^m> form_m.
"""
from fractions import Fraction as Fr
from functools import lru_cache
from tb import wedge, fadd, fscale


class Dual:
    __slots__ = ("v", "g")

    def __init__(self, v, g=None):
        self.v = Fr(v)
        self.g = dict(g) if g else {}

    @staticmethod
    def var(v, i):
        return Dual(v, {i: Fr(1)})

    def _c(self, o):
        return o if isinstance(o, Dual) else Dual(o)

    def __add__(self, o):
        o = self._c(o)
        g = dict(self.g)
        for k, x in o.g.items():
            g[k] = g.get(k, 0) + x
        return Dual(self.v + o.v, g)
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, {k: -x for k, x in self.g.items()})

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        g = {k: x * o.v for k, x in self.g.items()}
        for k, x in o.g.items():
            g[k] = g.get(k, 0) + self.v * x
        return Dual(self.v * o.v, g)
    __rmul__ = __mul__

    def inv(self):
        return Dual(1 / self.v, {k: -x / self.v ** 2 for k, x in self.g.items()})

    def __truediv__(self, o):
        return self * self._c(o).inv()

    def __rtruediv__(self, o):
        return self._c(o) * self.inv()

    def form(self):
        """the one-form d(self) as a dict {bitmask: coefficient}"""
        return {1 << k: x for k, x in self.g.items() if x != 0}


def dlog(x):
    return {m: c / x.v for m, c in x.form().items()}


# ------------------------------------------------------------- form-valued polynomials

def padd(P, Q, c=1):
    out = {k: dict(v) for k, v in P.items()}
    for m, f in Q.items():
        out[m] = fadd(out.get(m, {}), f, c)
        if not out[m]:
            del out[m]
    return out


def pmul(P, Q):
    out = {}
    for m1, f1 in P.items():
        for m2, f2 in Q.items():
            w = wedge(f1, f2)
            if not w:
                continue
            m = tuple(sorted(m1 + m2))
            out[m] = fadd(out.get(m, {}), w)
            if not out[m]:
                del out[m]
    return out


def pconst(form):
    return {(): dict(form)} if form else {}


def pvar(i, form):
    """x_i * form"""
    return {(i,): dict(form)} if form else {}


def pscale(P, c):
    return {m: fscale(f, c) for m, f in P.items() if c != 0}


class Gauss:
    def __init__(self, mean, cov):
        self.mean = [Fr(x) for x in mean]
        self.cov = [[Fr(x) for x in r] for r in cov]

        @lru_cache(maxsize=None)
        def mom(tup):
            if not tup:
                return Fr(1)
            first, rest = tup[0], tup[1:]
            tot = self.mean[first] * mom(rest)
            for k in range(len(rest)):
                c = self.cov[first][rest[k]]
                if c != 0:
                    tot += c * mom(rest[:k] + rest[k + 1:])
            return tot
        self.mom = mom

    def E(self, P):
        out = {}
        for m, f in P.items():
            x = self.mom(tuple(sorted(m)))
            if x != 0:
                out = fadd(out, f, x)
        return out


def block_gauss(blocks):
    """direct sum of independent Gaussian blocks [(mean, cov), ...]; returns Gauss and offsets"""
    mean, offs = [], []
    n = sum(len(m) for m, c in blocks)
    cov = [[Fr(0)] * n for _ in range(n)]
    o = 0
    for m, c in blocks:
        offs.append(o)
        mean += list(m)
        for i in range(len(m)):
            for j in range(len(m)):
                cov[o + i][o + j] = Fr(c[i][j])
        o += len(m)
    return Gauss(mean, cov), offs
