"""Truncated Laurent series in eps over Q with precision tracking, Gamma channels,
and the eps-expansion of Gamma monomials in Q[gamma, zeta(2), zeta(3), ...].

A value is a dict  channel -> Series,  channel = tuple(sorted((k, e))) meaning
prod_k Gamma(1 + k eps)^e.  Every Gamma(n + k eps) is rewritten as a rational
series times Gamma(1 + k eps).
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial

ORDER = 14          # highest eps power kept (absolute)


class Series:
    """sum_{i} c[i] eps^(val+i) + O(eps^prec)."""
    __slots__ = ("val", "c", "prec")

    def __init__(self, val, c, prec):
        self.val = val
        self.c = c
        self.prec = prec

    @staticmethod
    def const(x, prec=None):
        p = ORDER + 1 if prec is None else prec
        x = F(x)
        if x == 0:
            return Series(p, [], p)
        return Series(0, [x] + [F(0)] * max(0, p - 1), p)

    def trim(self):
        c, v = self.c, self.val
        i = 0
        while i < len(c) and c[i] == 0:
            i += 1
        c = c[i:]
        v += i
        p = min(self.prec, ORDER + 1)
        c = c[:max(0, p - v)]
        if not c:
            return Series(p, [], p)
        return Series(v, c, p)

    def coeff(self, n):
        if n >= self.prec:
            raise ArithmeticError("coefficient beyond precision: %d >= %d" % (n, self.prec))
        i = n - self.val
        return self.c[i] if 0 <= i < len(self.c) else F(0)

    def is_zero(self):
        return not any(self.c)

    def __add__(self, o):
        p = min(self.prec, o.prec)
        v = min(self.val, o.val)
        if v >= p:
            return Series(p, [], p)
        c = [F(0)] * (p - v)
        for i, x in enumerate(self.c):
            j = self.val + i - v
            if j < p - v:
                c[j] += x
        for i, x in enumerate(o.c):
            j = o.val + i - v
            if j < p - v:
                c[j] += x
        return Series(v, c, p).trim()

    def __neg__(self):
        return Series(self.val, [-x for x in self.c], self.prec)

    def __sub__(self, o):
        return self + (-o)

    def scale(self, x):
        x = F(x)
        if x == 0:
            return Series(self.prec, [], self.prec)
        return Series(self.val, [x * y for y in self.c], self.prec)

    def __mul__(self, o):
        if not self.c or not o.c:
            p = min(self.prec + (o.val if o.c else o.prec), o.prec + (self.val if self.c else self.prec))
            p = min(p, ORDER + 1)
            return Series(p, [], p)
        v = self.val + o.val
        p = min(self.prec + o.val, o.prec + self.val, ORDER + 1)
        n = p - v
        if n <= 0:
            return Series(p, [], p)
        c = [F(0)] * n
        a, b = self.c, o.c
        for i, x in enumerate(a):
            if i >= n:
                break
            if x == 0:
                continue
            lim = min(len(b), n - i)
            for j in range(lim):
                y = b[j]
                if y:
                    c[i + j] += x * y
        return Series(v, c, p).trim()

    def shift(self, m):
        """multiply by eps^m"""
        return Series(self.val + m, list(self.c), self.prec + m).trim()

    def inverse(self):
        s = self.trim()
        if not s.c:
            raise ZeroDivisionError("inverse of zero series")
        # s = eps^v (c0 + c1 eps + ...); known to prec -> relative precision r = prec - v
        v = s.val
        r = s.prec - v
        c = s.c + [F(0)] * max(0, r - len(s.c))
        inv = [F(1) / c[0]]
        for n in range(1, r):
            acc = F(0)
            for j in range(1, n + 1):
                if j < len(c) and c[j]:
                    acc += c[j] * inv[n - j]
            inv.append(-acc / c[0])
        return Series(-v, inv, -v + r).trim()

    def __repr__(self):
        terms = ["%s*e^%d" % (x, self.val + i) for i, x in enumerate(self.c) if x]
        return "(" + " + ".join(terms) + " + O(e^%d))" % self.prec


def linear(c0, c1):
    """the exact series c0 + c1 eps"""
    p = ORDER + 1
    c0, c1 = F(c0), F(c1)
    if c0 != 0:
        return Series(0, [c0, c1] + [F(0)] * (p - 2), p)
    if c1 != 0:
        return Series(1, [c1] + [F(0)] * (p - 2), p)
    raise ZeroDivisionError("zero linear factor")


HALF = 1000      # channel key HALF + k  stands for Gamma(1/2 + k eps)


@lru_cache(maxsize=None)
def gamma_ratio(n, k):
    """Gamma(n + k eps) / Gamma(base + k eps) as a Series, base = 1 (n integer) or 1/2
    (n half-integer).  k = 0: Gamma(n) / Gamma(base)."""
    n = F(n)
    if n.denominator == 1:
        n = int(n)
        if k == 0:
            if n <= 0:
                raise ZeroDivisionError("Gamma pole at integer %d" % n)
            return Series.const(factorial(n - 1))
        s = Series.const(1)
        if n >= 1:
            for j in range(1, n):
                s = s * linear(j, k)
        else:
            den = Series.const(1)
            for j in range(n, 1):
                den = den * linear(j, k)
            s = den.inverse()
        return s
    if n.denominator != 2:
        raise ValueError("Gamma argument must be integer or half-integer")
    # n = m + 1/2
    m = int(n - F(1, 2))
    s = Series.const(1)
    if m >= 0:
        for j in range(m):
            s = s * linear(F(1, 2) + j, k)
    else:
        den = Series.const(1)
        for j in range(m, 0):
            den = den * linear(F(1, 2) + j, k)
        s = den.inverse()
    return s


def value_const(x):
    return {(): Series.const(x)}


def vmul(A, B):
    out = {}
    for ca, sa in A.items():
        for cb, sb in B.items():
            d = dict(ca)
            for k, e in cb:
                d[k] = d.get(k, 0) + e
            key = tuple(sorted((k, e) for k, e in d.items() if e))
            s = sa * sb
            if key in out:
                out[key] = out[key] + s
            else:
                out[key] = s
    return out


def vadd(A, B, scale=1):
    out = dict(A)
    for k, s in B.items():
        t = s.scale(scale) if scale != 1 else s
        out[k] = out[k] + t if k in out else t
    return out


def vscale_series(A, s):
    return {k: v * s for k, v in A.items()}


def gamma_value(n, k, power=1):
    """Gamma(n + k eps)^power as a value"""
    n = F(n)
    r = gamma_ratio(n, k)
    if power == -1:
        r = r.inverse()
    elif power != 1:
        raise ValueError
    if n.denominator == 1:
        key = ((k, power),) if k != 0 else ()
    else:
        key = ((HALF + k, power),)
    return {key: r}


# ---------- eps-expansion of Gamma monomials ----------
# polynomials in symbols: 'g' (Euler gamma), 'z2','z3','z5','z7','z9','z11','z13'
EVEN = {4: F(2, 5), 6: F(8, 35), 8: F(24, 175), 10: F(32, 385)}   # zeta(2m)/zeta(2)^m


def pmul(a, b):
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            d = dict(ma)
            for s, e in mb:
                d[s] = d.get(s, 0) + e
            key = tuple(sorted(d.items()))
            out[key] = out.get(key, 0) + ca * cb
    return {k: v for k, v in out.items() if v}


def padd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + s * v
    return {k: v for k, v in out.items() if v}


def zeta_poly(n):
    if n % 2:
        return {(("z%d" % n, 1),): F(1)}
    if n == 2:
        return {(("z2", 1),): F(1)}
    return {(("z2", n // 2),): EVEN[n]}


@lru_cache(maxsize=None)
def gamma_expansion(channel, depth):
    """prod Gamma(1+k eps)^e * prod Gamma(1/2+k eps)^e as list of polynomials p_0..p_depth.
    Gamma(1/2+x) = sqrt(pi) exp(-(g + 2 L2) x + sum_{n>=2} (-1)^n (2^n-1) zeta(n) x^n / n)."""
    full = [(k, e) for k, e in channel if k < HALF // 2]
    half = [(k - HALF, e) for k, e in channel if k >= HALF // 2]
    logc = [dict() for _ in range(depth + 1)]
    S1 = sum(e * k for k, e in full)
    T1 = sum(e * k for k, e in half)
    if depth >= 1:
        d1 = {}
        if S1 + T1:
            d1[(("g", 1),)] = F(-(S1 + T1))
        if T1:
            d1[(("L2", 1),)] = F(-2 * T1)
        logc[1] = d1
    for n in range(2, depth + 1):
        Sn = sum(e * k ** n for k, e in full) + (2 ** n - 1) * sum(e * k ** n for k, e in half)
        if Sn:
            logc[n] = {m: c * F((-1) ** n * Sn, n) for m, c in zeta_poly(n).items()}
    # exp of series with zero constant term
    out = [dict() for _ in range(depth + 1)]
    nh = sum(e for k, e in half)
    out[0] = {((("sqrtpi", nh),) if nh else ()): F(1)}
    # E' = L' E  ->  n E_n = sum_{j=1}^n j L_j E_{n-j}
    for n in range(1, depth + 1):
        acc = {}
        for j in range(1, n + 1):
            if logc[j] and out[n - j]:
                acc = padd(acc, pmul({k: v * j for k, v in logc[j].items()}, out[n - j]))
        out[n] = {k: v / n for k, v in acc.items()}
    return tuple(tuple(sorted(o.items())) for o in out)


def finite_part(value, check_poles=True):
    """eps^0 coefficient of sum_channel series * Gamma-monomial, as a polynomial dict.
    Also returns the pole polynomials for checking."""
    total = {}
    poles = {}
    minval = min((s.val for s in value.values() if s.c), default=0)
    for ch, s in value.items():
        if not s.c:
            continue
        if any(k == 999999 for k, e in ch):
            if s.val < 1:
                raise ArithmeticError("finite part depends on an unknown t1star term (val %d)" % s.val)
            continue
        if s.prec < 1:
            raise ArithmeticError("insufficient precision %d in channel %s" % (s.prec, ch))
        depth = max(0, -s.val)
        exp = gamma_expansion(ch, depth)
        for i, x in enumerate(s.c):
            p = s.val + i
            if p > 0:
                break
            # contributes to eps^(p + j) with j from exp; need p + j = 0 and poles p+j<0
            for j in range(0, depth + 1):
                tgt = p + j
                if tgt > 0:
                    break
                poly = dict(exp[j])
                if not poly:
                    continue
                term = {k: v * x for k, v in poly.items()}
                if tgt == 0:
                    total = padd(total, term)
                else:
                    poles.setdefault(tgt, {})
                    poles[tgt] = padd(poles[tgt], term)
    poles = {k: v for k, v in poles.items() if v}
    if check_poles and poles:
        raise ArithmeticError("uncancelled poles: %s" % poles)
    return total


def as_zeta_dict(poly):
    """map polynomial to {weight-label: Fraction} for monomials 1, z3, z5, z7, z3^2, ...;
    raises if gamma or zeta(2) survive."""
    out = {}
    for m, c in poly.items():
        if any(s in ("g", "z2") for s, e in m):
            raise ArithmeticError("gamma or zeta(2) survived: %s" % (poly,))
        key = "*".join("%s^%d" % (s, e) if e > 1 else s for s, e in m) or "1"
        out[key] = out.get(key, 0) + c
    return {k: v for k, v in out.items() if v}
