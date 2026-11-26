"""Single-valued polylogarithms P_w(z) in letters {0,1} (Brown; conventions of Schnetz arXiv:1302.6445).

Holomorphic L_w(z):  dL_{wa}/dz = L_w/(z-a),  L_{0^k} = log^k z / k!,  generating series L = sum L_w w,
dL = L (x0/z + x1/(z-1)).  Z = sum zeta_w w = L(1) (shuffle-regularized), convergent words start with 1, end with 0.
SV:  P(z) = Ltilde_{x0,x1'}(zbar) L_{x0,x1}(z),   Z_{-x0,-x1'} x1' Z_{-x0,-x1'}^{-1} = Z x1 Z^{-1},
     dP/dz = P (x0/z + x1/(z-1)),   dP/dzbar = (x0/zbar + x1'/(zbar-1)) P.
Functions: SVF = sum_{word, mzv-monomial} R(x,y) * mono * P_word(z),  x = z, y = zbar.
"""
from fractions import Fraction as F
from functools import lru_cache

import mzv
from mzv import shuffle as _shuffle
from rat import R, ONE, ZERO_DEN, factor_pow, pmul

import os
WMAX = int(os.environ.get('SV_WMAX', '9'))


def words(n):
    if n == 0:
        return ['']
    return [w + a for w in words(n - 1) for a in '01']


def all_words(maxlen):
    out = []
    for n in range(maxlen + 1):
        out += words(n)
    return out


# ---------------- noncommutative series {word: mzv element} ----------------
def s_add(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = mzv.add(out.get(k, {}), v, c)
        if not out[k]:
            del out[k]
    return out


def s_mul(a, b, maxlen):
    out = {}
    for u, cu in a.items():
        for v, cv in b.items():
            if len(u) + len(v) > maxlen:
                continue
            w = u + v
            out[w] = mzv.add(out.get(w, {}), mzv.mul(cu, cv))
    return {k: v for k, v in out.items() if v}


ONE_S = {'': {(): F(1)}}


# ---------------- regularized zeta_w ----------------
@lru_cache(maxsize=None)
def zeta_word(w):
    """regularized L_w(1) as MZV element (dict)"""
    if w == '':
        return {(): F(1)}
    if set(w) == {'0'} or set(w) == {'1'}:
        return {}
    # trailing 1s
    if w[-1] == '1':
        k = len(w) - len(w.rstrip('1'))
        v = w[:len(w) - k]
        out = {}
        for x, c in _shuffle(v, '1' * k).items():
            if x == w:
                assert c == 1
                continue
            out = mzv.add(out, zeta_word(x), -c)
        return out
    if w[0] == '0':
        k = len(w) - len(w.lstrip('0'))
        v = w[k:]
        out = {}
        for x, c in _shuffle('0' * k, v).items():
            if x == w:
                assert c == 1
                continue
            out = mzv.add(out, zeta_word(x), -c)
        return out
    # convergent: starts with 1, ends with 0
    d = w.count('1')
    s = mzv.comp(w[::-1])
    return mzv.scal(mzv.zeta(*s), (-1) ** d)


def Zseries(maxlen):
    out = {}
    for w in all_words(maxlen):
        z = zeta_word(w)
        if z:
            out[w] = z
    return out


def subst(series_of_letter, w, maxlen):
    """product over letters of w of series_of_letter[letter], truncated"""
    out = ONE_S
    for a in w:
        out = s_mul(out, series_of_letter[a], maxlen)
    return out


_X1P = {}


def x1prime(maxlen=WMAX):
    if maxlen in _X1P:
        return _X1P[maxlen]
    Z = Zseries(maxlen)
    Zinv = {w: mzv.scal(c, (-1) ** len(w)) for w, c in ((w[::-1], zeta_word(w)) for w in all_words(maxlen)) if c}
    # check Zinv: antipode
    rhs = s_mul(s_mul(Z, {'1': {(): F(1)}}, maxlen), Zinv, maxlen)
    x1p = {'1': {(): F(1)}}
    for m in range(2, maxlen + 1):
        # LHS = Zm x1p Zm^{-1}, Zm = Z_{-x0,-x1p}; compute with current x1p (length < m), read length-m part
        minus = {'0': {'0': {(): F(-1)}}, '1': {k: mzv.scal(v, -1) for k, v in x1p.items()}}
        Zm = {}
        for w, c in Z.items():
            for u, cu in subst(minus, w, m).items():
                Zm[u] = mzv.add(Zm.get(u, {}), mzv.mul(c, cu))
        Zm = {k: v for k, v in Zm.items() if v}
        # inverse of Zm (group-like? not nec. after substitution of non-primitive; use series inverse)
        Zminv = s_inverse(Zm, m)
        lhs = s_mul(s_mul(Zm, x1p, m), Zminv, m)
        for w in words(m):
            d = mzv.add(rhs.get(w, {}), lhs.get(w, {}), -1)
            if d:
                x1p[w] = mzv.add(x1p.get(w, {}), d)
    _X1P[maxlen] = x1p
    return x1p


def s_inverse(a, maxlen):
    """inverse of a series with constant term 1"""
    assert a.get('') == {(): F(1)}
    b = {k: v for k, v in a.items() if k}
    out = dict(ONE_S)
    term = dict(ONE_S)
    for k in range(1, maxlen + 1):
        term = s_mul(term, b, maxlen)
        if not term:
            break
        out = s_add(out, term, (-1) ** k)
    return out


_SUB = {}


def subcoef(maxlen=WMAX):
    """C[u][s] = [ subst(s reversed, x1 -> x1') ]_u  for |s| <= |u| <= maxlen"""
    if maxlen in _SUB:
        return _SUB[maxlen]
    x1p = x1prime(maxlen)
    letter = {'0': {'0': {(): F(1)}}, '1': x1p}
    memo = {'': ONE_S}
    for t in all_words(maxlen):
        if t:
            memo[t] = s_mul(memo[t[:-1]], letter[t[-1]], maxlen)
    C = {}
    for t, ser in memo.items():
        s = t[::-1]
        for u, c in ser.items():
            C.setdefault(u, {})[s] = c
    _SUB[maxlen] = C
    return C


def tensor(w):
    """P_w = sum T[(s, v)] L_s(zbar) L_v(z)"""
    if len(w) > WMAX:
        raise OverflowError('word length %d > SV_WMAX = %d' % (len(w), WMAX))
    C = subcoef()
    out = {}
    for k in range(len(w) + 1):
        u, v = w[:k], w[k:]
        for s, c in C.get(u, {}).items():
            out[(s, v)] = mzv.add(out.get((s, v), {}), c)
    return {k: v for k, v in out.items() if v}


@lru_cache(maxsize=None)
def value_at_1(w):
    """regularized P_w(1) as MZV element"""
    out = {}
    for (s, v), c in tensor(w).items():
        out = mzv.add(out, mzv.mul(c, mzv.mul(zeta_word(s), zeta_word(v))))
    return out


@lru_cache(maxsize=None)
def swap_word(w):
    """sigma P_w (z <-> zbar) = sum_s C[w][s] P_s"""
    if len(w) > WMAX:
        raise OverflowError('word length %d > SV_WMAX = %d' % (len(w), WMAX))
    C = subcoef()
    return {s: c for s, c in C.get(w, {}).items() if c}


@lru_cache(maxsize=None)
def dzb_word(w):
    """d/dzbar P_w = sum over (word, mzv, kind) with kind 0: 1/zbar, 1: 1/(zbar-1)  -> list of (u, c, a)"""
    if len(w) > WMAX:
        raise OverflowError('word length %d > SV_WMAX = %d' % (len(w), WMAX))
    out = []
    if not w:
        return ()
    if w[0] == '0':
        out.append((w[1:], {(): F(1)}, '0'))
    x1p = x1prime()
    for k in range(1, len(w) + 1):
        v = w[:k]
        if v in x1p:
            out.append((w[k:], x1p[v], '1'))
    return tuple(out)


# ---------------- SVF algebra ----------------
def _addto(t, w, m, r):
    if r.iszero():
        return
    d = t.setdefault(w, {})
    if m in d:
        s = d[m] + r
        if s.iszero():
            del d[m]
            if not d:
                del t[w]
        else:
            d[m] = s
    else:
        d[m] = r


class SVF:
    __slots__ = ('t',)

    def __init__(self, t=None):
        self.t = t if t is not None else {}

    @staticmethod
    def rat(r, w='', m=()):
        f = SVF()
        _addto(f.t, w, m, r)
        return f

    def copy(self):
        return SVF({w: dict(d) for w, d in self.t.items()})

    def iszero(self):
        return not self.t

    def terms(self):
        for w, d in self.t.items():
            for m, r in d.items():
                yield w, m, r

    def __add__(self, other):
        out = self.copy()
        for w, m, r in other.terms():
            _addto(out.t, w, m, r)
        return out

    def __sub__(self, other):
        return self + other.scale(-1)

    def scale(self, c):
        out = SVF()
        for w, m, r in self.terms():
            _addto(out.t, w, m, r.scale(c))
        return out

    def mulmzv(self, z):
        out = SVF()
        for w, m, r in self.terms():
            for m2, c in z.items():
                _addto(out.t, w, tuple(sorted(m + m2)), r.scale(c))
        return out

    def mulrat(self, s):
        out = SVF()
        for w, m, r in self.terms():
            _addto(out.t, w, m, r * s)
        return out

    def __mul__(self, other):
        if isinstance(other, R):
            return self.mulrat(other)
        out = SVF()
        for w1, m1, r1 in self.terms():
            for w2, m2, r2 in other.terms():
                rr = r1 * r2
                m = tuple(sorted(m1 + m2))
                for w, c in _shuffle(w1, w2).items():
                    _addto(out.t, w, m, rr.scale(c))
        return out

    def swap(self):
        out = SVF()
        for w, m, r in self.terms():
            rs = r.swap()
            for s, c in swap_word(w).items():
                for m2, cc in c.items():
                    _addto(out.t, s, tuple(sorted(m + m2)), rs.scale(cc))
        return out

    def dz(self):
        out = SVF()
        for w, m, r in self.terms():
            _addto(out.t, w, m, r.dx())
            if w:
                a = w[-1]
                _addto(out.t, w[:-1], m, r.divfactor(0 if a == '0' else 2))
        return out

    def dzb(self):
        out = SVF()
        for w, m, r in self.terms():
            _addto(out.t, w, m, r.dy())
            for u, c, a in dzb_word(w):
                rr = r.divfactor(1 if a == '0' else 3)
                for m2, cc in c.items():
                    _addto(out.t, u, tuple(sorted(m + m2)), rr.scale(cc))
        return out

    def maxweight(self):
        return max((len(w) + mzv.weight(m) for w, m, r in self.terms()), default=-1)

    def has_w_poles(self):
        return any(r.den[4] > 0 for w, m, r in self.terms())

    def int_z(self):
        out = SVF()
        for w, m, r in self.terms():
            for (u, rr) in _int_z_term(w, r):
                _addto(out.t, u, m, rr)
        return out

    def int_zb(self):
        out = SVF()
        for w, m, r in self.terms():
            for (u, m2, rr) in _int_zb_term(w, r):
                _addto(out.t, u, tuple(sorted(m + m2)), rr)
        return out

    def evaluate(self, z, dps=30):
        import mpmath as mp
        tot = mp.mpc(0)
        zb = mp.conj(z)
        for w, m, r in self.terms():
            rv = r.num and _reval(r, z, zb)
            tot += rv * mzv.numeric({m: F(1)}, dps) * Pnum(w, z, dps)
        return tot

    def __repr__(self):
        return 'SVF(%d words, %d terms)' % (len(self.t), sum(len(d) for d in self.t.values()))


def _reval(r, x, y):
    a, b, c, d, e = r.den
    num = sum((v.numerator / __import__('mpmath').mpf(v.denominator)) * x ** i * y ** j for (i, j), v in r.num.items())
    return num / (x ** a * y ** b * (x - 1) ** c * (y - 1) ** d * (x - y) ** e)


class Obstruction(Exception):
    pass


def _int_z_term(w, r):
    """holomorphic SV primitive of r * P_w; returns list of (word, R)"""
    if r.den[4] > 0:
        raise Obstruction('int_z: (z - zbar) pole in integrand')
    poly, A, B = r.pf_x()
    out = []
    if 1 in A:
        out.append((w + '0', A[1]))
    if 1 in B:
        out.append((w + '1', B[1]))
    # rational primitive rho of the remaining part
    rho = R({})
    for k, c in poly.items():
        rho = rho + c * R({(k + 1, 0): F(1, k + 1)})
    for j, c in A.items():
        if j >= 2:
            rho = rho + c * R({(0, 0): F(-1, j - 1)}, (j - 1, 0, 0, 0, 0))
    for j, c in B.items():
        if j >= 2:
            rho = rho + c * R({(0, 0): F(-1, j - 1)}, (0, 0, j - 1, 0, 0))
    if rho.iszero():
        return out
    out.append((w, rho))
    if w:
        a = w[-1]
        rr = rho.divfactor(0 if a == '0' else 2).scale(-1)
        out += _int_z_term(w[:-1], rr)
    return out


def _int_zb_term(w, r):
    """antiholomorphic SV primitive of r * P_w; returns list of (word, mzv-monomial, R)"""
    if r.den[4] > 0:
        raise Obstruction('int_zb: (z - zbar) pole in integrand')
    poly, A, B = r.swap().pf_x()  # in swapped variables: y <-> x
    poly = {k: c.swap() for k, c in poly.items()}
    A = {k: c.swap() for k, c in A.items()}
    B = {k: c.swap() for k, c in B.items()}
    out = []
    if 1 in A:
        out += _prim_letter(w, '0', A[1])
    if 1 in B:
        out += _prim_letter(w, '1', B[1])
    rho = R({})
    for k, c in poly.items():
        rho = rho + c * R({(0, k + 1): F(1, k + 1)})
    for j, c in A.items():
        if j >= 2:
            rho = rho + c * R({(0, 0): F(-1, j - 1)}, (0, j - 1, 0, 0, 0))
    for j, c in B.items():
        if j >= 2:
            rho = rho + c * R({(0, 0): F(-1, j - 1)}, (0, 0, 0, j - 1, 0))
    if rho.iszero():
        return out
    out.append((w, (), rho))
    for u, c, a in dzb_word(w):
        rr = rho.divfactor(1 if a == '0' else 3).scale(-1)
        for m2, cc in c.items():
            for (u2, m3, r3) in _int_zb_term(u, rr.scale(cc)):
                out.append((u2, tuple(sorted(m2 + m3)), r3))
    return out


def _prim_letter(w, a, coef):
    """antiholomorphic primitive of coef(x) * P_w / (zbar - a), coef independent of y"""
    out = [(a + w, (), coef)]
    for u, c, b in dzb_word(a + w):
        if u == w and b == a and c == {(): F(1)}:
            continue  # main term
        # correction: c * P_u / (zbar - b) must be subtracted (integrated)
        for m2, cc in c.items():
            for (u2, m3, r3) in _prim_letter(u, b, coef.scale(-cc)):
                out.append((u2, tuple(sorted(m2 + m3)), r3))
    return out


# ---------------- numerics ----------------
_PCACHE = {}


def Lnum_series(w, N):
    """L_w(z) = sum_j log(z)^j * sum_n c[j][n] z^n   (exact Fractions)"""
    key = (w, N)
    if key in _PCACHE:
        return _PCACHE[key]
    if w == '':
        res = {0: {0: F(1)}}
    else:
        f = Lnum_series(w[:-1], N)
        a = w[-1]
        res = {}

        def addc(j, n, c):
            if c:
                res.setdefault(j, {})
                res[j][n] = res[j].get(n, 0) + c
        if a == '0':
            for j, ser in f.items():
                for n, c in ser.items():
                    if n == 0:
                        addc(j + 1, 0, c / (j + 1))
                    else:
                        _intlog(j, n - 1, c, addc)
        else:
            for j, ser in f.items():
                for n, c in ser.items():
                    for m in range(0, N - n):
                        _intlog(j, n + m, -c, addc)
    _PCACHE[key] = res
    return res


def _intlog(j, N, c, addc):
    """int_0^z log^j(t) t^N dt = sum_k (-1)^k j!/(j-k)! log^{j-k} z * z^{N+1}/(N+1)^{k+1}"""
    fac = 1
    for k in range(j + 1):
        if k:
            fac *= (j - k + 1)
        addc(j - k, N + 1, c * (-1) ** k * fac / F(N + 1) ** (k + 1))


def Lnum(w, z, N=60):
    import mpmath as mp
    ser = Lnum_series(w, N)
    lz = mp.log(z)
    tot = mp.mpc(0)
    for j, s in ser.items():
        tot += lz ** j * mp.fsum(mp.mpf(c.numerator) / c.denominator * z ** n for n, c in s.items())
    return tot


def Pnum(w, z, dps=30, N=None):
    import mpmath as mp
    if N is None:
        N = int(dps * 3.4 / max(1e-9, -mp.log10(abs(z)))) + 10 if abs(z) < 1 else 200
    zb = mp.conj(z)
    tot = mp.mpc(0)
    for (s, v), c in tensor(w).items():
        tot += mzv.numeric(c, dps) * Lnum(s, zb, N) * Lnum(v, z, N)
    return tot


# ---------------- S3 transformations ----------------
@lru_cache(maxsize=None)
def trans_1mz(w):
    """P_w(1-z) = sum_v m[v] P_v(z),  m[v] MZV"""
    if w == '':
        return {'': {(): F(1)}}
    prev = trans_1mz(w[:-1])
    ab = '1' if w[-1] == '0' else '0'
    out = {}
    for v, c in prev.items():
        out[v + ab] = mzv.add(out.get(v + ab, {}), c)
    const = {}
    for v, c in out.items():
        const = mzv.add(const, mzv.mul(c, value_at_1(v)), -1)
    if const:
        out[''] = mzv.add(out.get('', {}), const)
    return {k: v for k, v in out.items() if v}


@lru_cache(maxsize=None)
def trans_inv(w):
    """P_w(1/z) = sum_v m[v] P_v(z)"""
    if w == '':
        return {'': {(): F(1)}}
    prev = trans_inv(w[:-1])
    out = {}
    for v, c in prev.items():
        if w[-1] == '0':
            out[v + '0'] = mzv.add(out.get(v + '0', {}), c, -1)
        else:
            out[v + '1'] = mzv.add(out.get(v + '1', {}), c)
            out[v + '0'] = mzv.add(out.get(v + '0', {}), c, -1)
    out = {k: v for k, v in out.items() if v}
    const = dict(value_at_1(w))
    for v, c in out.items():
        const = mzv.add(const, mzv.mul(c, value_at_1(v)), -1)
    if const:
        out[''] = mzv.add(out.get('', {}), const)
    return {k: v for k, v in out.items() if v}
