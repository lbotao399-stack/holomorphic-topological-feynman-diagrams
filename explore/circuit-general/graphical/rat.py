"""Rational functions in x (= z) and y (= zbar) with denominators x^a y^b (x-1)^c (y-1)^d (x-y)^e.

R = (num, den): num = {(i, j): Fraction} polynomial, den = (a, b, c, d, e) >= 0; reduced form
(num not divisible by x, y, x-1, y-1, x-y wherever the exponent is > 0) is unique.
"""
from fractions import Fraction as F

ZERO_DEN = (0, 0, 0, 0, 0)


# ---------------- polynomials {(i,j): c} ----------------
def padd(p, q, c=1):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, 0) + c * v
        if not out[k]:
            del out[k]
    return out


def pmul(p, q):
    out = {}
    for (i1, j1), v1 in p.items():
        for (i2, j2), v2 in q.items():
            k = (i1 + i2, j1 + j2)
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v}


def pscal(p, c):
    if not c:
        return {}
    return {k: v * c for k, v in p.items()}


_FACT = {0: {(1, 0): 1}, 1: {(0, 1): 1}, 2: {(1, 0): 1, (0, 0): -1}, 3: {(0, 1): 1, (0, 0): -1},
         4: {(1, 0): 1, (0, 1): -1}}
_POW = {}


def factor_pow(f, k):
    key = (f, k)
    if key not in _POW:
        p = {(0, 0): 1}
        for _ in range(k):
            p = pmul(p, _FACT[f])
        _POW[key] = p
    return _POW[key]


def div_x(p):
    if any(i == 0 for (i, j) in p):
        return None
    return {(i - 1, j): v for (i, j), v in p.items()}


def div_y(p):
    if any(j == 0 for (i, j) in p):
        return None
    return {(i, j - 1): v for (i, j), v in p.items()}


def _div_lin_x(p, a):
    """divide by (x - a) as polynomial in x with coefficients in Q[y]; None if remainder"""
    # group by j
    byj = {}
    for (i, j), v in p.items():
        byj.setdefault(j, {})[i] = v
    out = {}
    for j, coeffs in byj.items():
        deg = max(coeffs)
        # synthetic division
        q = {}
        carry = 0
        for i in range(deg, -1, -1):
            c = coeffs.get(i, 0) + carry
            if i == 0:
                if c:
                    return None
            else:
                if c:
                    q[i - 1] = c
                carry = c * a
        for i, v in q.items():
            out[(i, j)] = v
    return out


def _swap(p):
    return {(j, i): v for (i, j), v in p.items()}


def div_xm1(p):
    return _div_lin_x(p, 1)


def div_ym1(p):
    r = _div_lin_x(_swap(p), 1)
    return None if r is None else _swap(r)


def div_xmy(p):
    """divide by (x - y)"""
    # substitute x = y + t: p(y+t, y) ; divisible iff no t^0 part. Do division via ordering in x.
    # Use long division in x over Q[y]: leading x-degree term
    rem = dict(p)
    q = {}
    while rem:
        dx = max(i for (i, j) in rem)
        if dx == 0:
            return None
        # take all terms with x-degree dx: c(y) x^dx ; quotient term c(y) x^(dx-1)
        lead = {(i, j): v for (i, j), v in rem.items() if i == dx}
        for (i, j), v in lead.items():
            q[(i - 1, j)] = q.get((i - 1, j), 0) + v
            # subtract v x^(i-1) y^j (x - y) = v x^i y^j - v x^(i-1) y^(j+1)
            rem[(i, j)] = rem.get((i, j), 0) - v
            if not rem[(i, j)]:
                del rem[(i, j)]
            k = (i - 1, j + 1)
            rem[k] = rem.get(k, 0) + v
            if not rem[k]:
                del rem[k]
    return {k: v for k, v in q.items() if v}


_DIVS = (div_x, div_y, div_xm1, div_ym1, div_xmy)


def peval(p, x, y):
    return sum(v * x ** i * y ** j for (i, j), v in p.items())


class R:
    __slots__ = ('num', 'den')

    def __init__(self, num, den=ZERO_DEN, reduce=True):
        self.num = {k: (v if isinstance(v, F) else F(v)) for k, v in num.items() if v}
        self.den = tuple(den)
        if reduce:
            self._reduce()

    def _reduce(self):
        if not self.num:
            self.den = ZERO_DEN
            return
        den = list(self.den)
        num = self.num
        for f in range(5):
            while den[f] > 0:
                q = _DIVS[f](num)
                if q is None:
                    break
                num = q
                den[f] -= 1
        self.num = num
        self.den = tuple(den)

    @staticmethod
    def const(c):
        return R({(0, 0): F(c)})

    @staticmethod
    def mono(i=0, j=0, den=ZERO_DEN, c=1):
        return R({(i, j): F(c)}, den)

    def iszero(self):
        return not self.num

    def __eq__(self, other):
        return self.den == other.den and self.num == other.num

    def __hash__(self):
        return hash((self.den, tuple(sorted(self.num.items()))))

    def __add__(self, other):
        if not other.num:
            return self
        if not self.num:
            return other
        if self.den == other.den:
            return R(padd(self.num, other.num), self.den)
        den = tuple(max(a, b) for a, b in zip(self.den, other.den))
        return R(padd(self._lift(den), other._lift(den)), den)

    def __sub__(self, other):
        return self + other.scale(-1)

    def __neg__(self):
        return self.scale(-1)

    def _lift(self, den):
        p = self.num
        for f in range(5):
            k = den[f] - self.den[f]
            if k:
                p = pmul(p, factor_pow(f, k))
        return p

    def scale(self, c):
        if not c:
            return R({})
        return R({k: v * c for k, v in self.num.items()}, self.den, reduce=False)

    def __mul__(self, other):
        if isinstance(other, R):
            if not self.num or not other.num:
                return R({})
            return R(pmul(self.num, other.num), tuple(a + b for a, b in zip(self.den, other.den)))
        return self.scale(other)

    __rmul__ = __mul__

    def mulpoly(self, p):
        return R(pmul(self.num, p), self.den)

    def divfactor(self, f, k=1):
        d = list(self.den)
        d[f] += k
        return R(self.num, d)

    def swap(self):
        """f(x,y) -> f(y,x)"""
        a, b, c, d, e = self.den
        num = _swap(self.num)
        if e % 2:
            num = {k: -v for k, v in num.items()}
        return R(num, (b, a, d, c, e), reduce=False)

    def dx(self):
        # d/dx (N / D) = N_x / D - N * sum_f a_f f_x / f / D
        a, b, c, d, e = self.den
        N = self.num
        Nx = {(i - 1, j): v * i for (i, j), v in N.items() if i}
        newden = (a + 1, b, c + 1, d, e + 1)
        # N_x * x (x-1) (x-y) - N * [a (x-1)(x-y) + c x (x-y) + e x (x-1)]
        t = pmul(Nx, pmul(factor_pow(0, 1), pmul(factor_pow(2, 1), factor_pow(4, 1))))
        if a:
            t = padd(t, pmul(N, pmul(factor_pow(2, 1), factor_pow(4, 1))), -a)
        if c:
            t = padd(t, pmul(N, pmul(factor_pow(0, 1), factor_pow(4, 1))), -c)
        if e:
            t = padd(t, pmul(N, pmul(factor_pow(0, 1), factor_pow(2, 1))), -e)
        return R(t, newden)

    def dy(self):
        return self.swap().dx().swap()

    def evaluate(self, x, y):
        a, b, c, d, e = self.den
        return peval(self.num, x, y) / (x ** a * y ** b * (x - 1) ** c * (y - 1) ** d * (x - y) ** e)

    def __repr__(self):
        return 'R(%s / x^%d y^%d (x-1)^%d (y-1)^%d (x-y)^%d)' % ((self.num,) + self.den)

    # ----- partial fractions in x (no (x-y) poles allowed) -----
    def pf_x(self):
        """self = sum_k Px[k] x^k + sum_{j>=1} A[j] x^{-j} + sum_{j>=1} B[j] (x-1)^{-j}
        with coefficients rational in y (R objects with x-degree 0). Requires den[4] == 0.
        returns (poly, A, B) dicts k -> R."""
        a, b, c, d, e = self.den
        assert e == 0, 'pf_x with (x-y) pole'
        # num(x,y) / (x^a (x-1)^c) * 1/(y^b (y-1)^d)
        yden = (0, b, 0, d, 0)
        byj = {}
        for (i, j), v in self.num.items():
            byj.setdefault(i, {})
            byj[i][(0, j)] = v
        # numerator as polynomial in x with coefficients in Q[y]: coef[i] = poly in y
        coef = {i: p for i, p in byj.items()}
        poly, A, B = {}, {}, {}
        # Laurent at x=0 part: we do it by generic algorithm on univariate-in-x with poly-in-y coefficients
        # N(x)/(x^a (x-1)^c) : write N = Q x^a (x-1)^c + rem; use expansions
        # Step 1: principal part at x = 0 of order a: expand 1/(x-1)^c around 0
        # (x-1)^{-c} = (-1)^c (1-x)^{-c} = (-1)^c sum_n binom(n+c-1, c-1) x^n
        from math import comb
        degN = max(coef) if coef else 0

        def ycoef(i):
            return coef.get(i, {})
        if a:
            for j in range(1, a + 1):
                # coefficient of x^{-j}: sum over i + n = a - j of N_i * s_n, s_n series of (x-1)^{-c}
                acc = {}
                for i in range(0, a - j + 1):
                    n = a - j - i
                    s = (-1) ** c * (comb(n + c - 1, c - 1) if c else (1 if n == 0 else 0))
                    if s and ycoef(i):
                        acc = padd(acc, ycoef(i), s)
                if acc:
                    A[j] = acc
        if c:
            # expand around x = 1: x = 1 + t ; x^{-a} = sum_n binom(-a, n) t^n ; N(1+t)
            Nt = {}
            for i, p in coef.items():
                for k in range(i + 1):
                    cc = comb(i, k)
                    if cc:
                        Nt[k] = padd(Nt.get(k, {}), p, cc)
            for j in range(1, c + 1):
                acc = {}
                for k in range(0, c - j + 1):
                    n = c - j - k
                    # binom(-a, n) = (-1)^n binom(a+n-1, n)
                    s = (-1) ** n * (comb(a + n - 1, n) if a else (1 if n == 0 else 0))
                    if s and Nt.get(k):
                        acc = padd(acc, Nt[k], s)
                if acc:
                    B[j] = acc
        # polynomial part: N - sum principal parts * x^a (x-1)^c, divided exactly
        remnum = {(i, j): v for (i, j), v in self.num.items()}
        # subtract principal parts: A[j] x^{-j}: multiply by x^a (x-1)^c -> A[j] x^{a-j} (x-1)^c
        for j, r in A.items():
            p = pmul({(a - j, 0): 1}, factor_pow(2, c))
            remnum = padd(remnum, pmul(r, p), -1)
        for j, r in B.items():
            p = pmul({(0, 0): 1}, pmul(factor_pow(0, a), factor_pow(2, c - j)))
            remnum = padd(remnum, pmul(r, p), -1)
        # now remnum divisible by x^a (x-1)^c
        q = remnum
        for _ in range(a):
            q = div_x(q) if q else q
            assert q is not None
        for _ in range(c):
            q = div_xm1(q) if q else q
            assert q is not None
        if q:
            byi = {}
            for (i, j), v in q.items():
                byi.setdefault(i, {})[(0, j)] = v
            for i, p in byi.items():
                poly[i] = R(p, yden)
        A = {j: R(p, yden) for j, p in A.items()}
        B = {j: R(p, yden) for j, p in B.items()}
        return poly, A, B


def x_power(k):
    """x^k as R (k may be negative)"""
    return R({(k, 0): 1}) if k >= 0 else R({(0, 0): 1}, (-k, 0, 0, 0, 0))


def y_power(k):
    return R({(0, k): 1}) if k >= 0 else R({(0, 0): 1}, (0, -k, 0, 0, 0))


def xm1_power(k):
    return R(factor_pow(2, k)) if k >= 0 else R({(0, 0): 1}, (0, 0, -k, 0, 0))


def ym1_power(k):
    return R(factor_pow(3, k)) if k >= 0 else R({(0, 0): 1}, (0, 0, 0, -k, 0))


def w_power(k):
    return R(factor_pow(4, k)) if k >= 0 else R({(0, 0): 1}, (0, 0, 0, 0, -k))


ONE = R({(0, 0): 1})
X = R({(1, 0): 1})
Y = R({(0, 1): 1})
