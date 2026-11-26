"""Moebius substitutions, local expansions at 0, 1, infinity, and integration over C of SV functions.

expand0(F, D): F pole-free in (z - zbar); returns {(l, m, mb): mzv} with F = sum c log^l(z zbar) z^m zbar^mb,
total degree m + mb <= D.  Expansion at s=1 uses u = 1-z, at infinity v = 1/z (log|v|^2 = -log|z|^2).
"""
from fractions import Fraction as F
from math import comb

import mzv
import sv
from rat import R, padd, pmul, factor_pow
from sv import SVF, _addto


# ---------------- Moebius substitutions on R ----------------
def r_1mx(r):
    """r(1-x, 1-y)"""
    a, b, c, d, e = r.den
    num = {}
    for (i, j), v in r.num.items():
        # (1-x)^i (1-y)^j
        for k in range(i + 1):
            for l in range(j + 1):
                key = (k, l)
                num[key] = num.get(key, 0) + v * comb(i, k) * comb(j, l) * (-1) ** (k + l)
    sign = (-1) ** (a + b + c + d + e)
    num = {k: v * sign for k, v in num.items() if v}
    return R(num, (c, d, a, b, e))


def r_inv(r):
    """r(1/x, 1/y)"""
    a, b, c, d, e = r.den
    if not r.num:
        return R({})
    Dx = max(i for (i, j) in r.num)
    Dy = max(j for (i, j) in r.num)
    sx = a + c + e
    sy = b + d + e
    num = {(Dx - i + sx, Dy - j + sy): v * (-1) ** (c + d + e) for (i, j), v in r.num.items()}
    # remaining denominator x^Dx y^Dy (x-1)^c (y-1)^d (x-y)^e
    return R(num, (Dx, Dy, c, d, e))


def svf_1mz(f):
    out = SVF()
    for w, m, r in f.terms():
        rr = r_1mx(r)
        for v, c in sv.trans_1mz(w).items():
            for m2, cc in c.items():
                _addto(out.t, v, tuple(sorted(m + m2)), rr.scale(cc))
    return out


def svf_inv(f):
    out = SVF()
    for w, m, r in f.terms():
        rr = r_inv(r)
        for v, c in sv.trans_inv(w).items():
            for m2, cc in c.items():
                _addto(out.t, v, tuple(sorted(m + m2)), rr.scale(cc))
    return out


# ---------------- Laurent expansion at 0 ----------------
def _rat_series(r, D):
    """Laurent series of r (no (x-y) pole) at 0: {(m, mb): Fraction}, total degree <= D"""
    a, b, c, d, e = r.den
    assert e == 0
    # (x-1)^{-c} = (-1)^c sum binom(n+c-1, c-1) x^n
    out = {}
    mind = min(i + j for (i, j) in r.num) - a - b
    nmax = D - mind
    if nmax < 0:
        return {}
    sx = [(-1) ** c * (comb(n + c - 1, c - 1) if c else int(n == 0)) for n in range(nmax + 1)]
    sy = [(-1) ** d * (comb(n + d - 1, d - 1) if d else int(n == 0)) for n in range(nmax + 1)]
    for (i, j), v in r.num.items():
        base = i + j - a - b
        for n1 in range(nmax + 1):
            if base + n1 > D or not sx[n1]:
                continue
            for n2 in range(nmax + 1 - n1):
                if base + n1 + n2 > D:
                    break
                if not sy[n2]:
                    continue
                k = (i - a + n1, j - b + n2)
                out[k] = out.get(k, 0) + v * sx[n1] * sy[n2]
    return {k: v for k, v in out.items() if v}


_LS = {}


def _L_series(v, D):
    """L_v(x) = sum_j log^j x * sum_n c x^n, n <= D: returns {(j, n): Fraction}"""
    key = (v, D)
    if key not in _LS:
        ser = sv.Lnum_series(v, max(D, 1) + 1)
        out = {}
        for j, s in ser.items():
            for n, c in s.items():
                if n <= D and c:
                    out[(j, n)] = c
        _LS[key] = out
    return _LS[key]


_PS = {}


def _P_series(w, D):
    """P_w = sum c lx^j ly^k x^n y^nb ; returns {(j, k, n, nb): mzv} total degree <= D"""
    key = (w, D)
    if key in _PS:
        return _PS[key]
    out = {}
    for (s, v), c in sv.tensor(w).items():
        Ls = _L_series(s, D)
        Lv = _L_series(v, D)
        for (k, nb), cs in Ls.items():
            for (j, n), cv in Lv.items():
                if n + nb > D:
                    continue
                kk = (j, k, n, nb)
                out[kk] = mzv.add(out.get(kk, {}), mzv.scal(c, cs * cv))
    out = {k: v for k, v in out.items() if v}
    _PS[key] = out
    return out


def expand0(f, D):
    """expansion of pole-free SVF f at 0: {(l, m, mb): mzv}, total degree <= D"""
    acc = {}
    for w, m, r in f.terms():
        mind = min(i + j for (i, j) in r.num) - r.den[0] - r.den[1]
        rs = _rat_series(r, D)
        if not rs:
            continue
        Dp = D - mind
        ps = _P_series(w, Dp)
        for (j, k, n, nb), c in ps.items():
            cm = mzv.mul(c, {m: F(1)})
            for (p, pb), rv in rs.items():
                if p + n + pb + nb > D:
                    continue
                key = (j, k, p + n, pb + nb)
                acc[key] = mzv.add(acc.get(key, {}), mzv.scal(cm, rv))
    # convert lx^j ly^k into Lambda^l = log(x y)^l : coefficient of Lambda^l = coefficient of lx^l ly^0
    out = {}
    bad = {}
    for (j, k, mm, mb), c in acc.items():
        if not c:
            continue
        if k == 0:
            out[(j, mm, mb)] = c
    # consistency: coefficient of lx^j ly^k must be binom(j+k, j) * out[(j+k)]
    for (j, k, mm, mb), c in acc.items():
        if not c or k == 0:
            continue
        expect = mzv.scal(out.get((j + k, mm, mb), {}), comb(j + k, j))
        if mzv.add(c, expect, -1):
            bad[(j, k, mm, mb)] = c
    if bad:
        raise ArithmeticError('expansion not single-valued: %s' % list(bad.items())[:3])
    return {k: v for k, v in out.items() if v}


def homog_divide(series, lin, times):
    """divide a series {(l, m, mb): c} (complete in total degree) `times` times by (x - y) if lin == 1
    (or by (y - x) if lin == -1); returns (quotient, remainders) where remainders lists nonzero
    obstructions (exactness failures)."""
    cur = series
    rems = []
    for _ in range(times):
        # group by (l, total degree)
        groups = {}
        for (l, m, mb), c in cur.items():
            groups.setdefault((l, m + mb), {})[m] = c
        new = {}
        for (l, D), coef in groups.items():
            ms = sorted(coef)
            q = {}
            prev = {}
            for mm in range(ms[0], ms[-1] + 1):
                cm = coef.get(mm, {})
                # q_m = q_{m-1} - c_m   (for division by (x - y)), coefficient x^m y^(D-1-m)
                qm = mzv.add(prev, cm, -1)
                if mm == ms[-1]:
                    if qm:
                        rems.append(((l, mm, D - mm), qm))
                    break
                if qm:
                    q[mm] = qm
                prev = qm
            for mm, c in q.items():
                c2 = c if lin == 1 else mzv.scal(c, -1)
                new[(l, mm, D - 1 - mm)] = c2
        cur = new
    return cur, rems


# ---------------- integration over C ----------------
def integrate_C(phi):
    """int_C phi d^2z / pi  (d^2z = dRe dIm) for pole-free SVF phi with integrable singularities"""
    H = phi.int_z()
    if H.has_w_poles():
        raise sv.Obstruction('integrate_C: (z-zbar) poles')
    e0 = expand0(H, 0)
    e1 = expand0(svf_1mz(H), 0)
    ei = expand0(svf_inv(H), 2)
    c0 = e0.get((0, 0, -1), {})
    c1 = mzv.scal(e1.get((0, 0, -1), {}), -1)
    ci = ei.get((0, 0, 1), {})
    # divergence checks: terms with m+mb+1 <= 0 and mb = m-1 other than (0,0,-1) must vanish
    for tag, e, cond in (('0', e0, lambda m, mb: m + mb + 1 < 0), ('1', e1, lambda m, mb: m + mb + 1 < 0)):
        for (l, m, mb), c in e.items():
            if mb == m - 1 and c and (cond(m, mb) or (m == 0 and l > 0)):
                raise ArithmeticError('integrate_C: divergent term at %s: %s' % (tag, (l, m, mb)))
    for (l, m, mb), c in ei.items():
        # in v = 1/z: z^M zb^(M-1) <-> v^-M vb^(1-M); divergent if M + (M-1) + 1 > 0 i.e. v-exponent m=-M<0
        if c and mb == m + 1 and (m < 0 or (m == 0 and l > 0)):
            raise ArithmeticError('integrate_C: divergent term at infinity: %s' % ((l, m, mb),))
    return mzv.add(mzv.add(ci, c0, -1), c1, -1)
