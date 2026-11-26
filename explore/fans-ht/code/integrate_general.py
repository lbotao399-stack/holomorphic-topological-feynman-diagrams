"""Exact integration of the (1,1) fan density in S1 coordinates for general sources:

    K(lam; 0) = int_{(0,1)^{2n-2} x Delta^2} 2^{-n/2} Pi(r) Haf[L_i M_j] d^{2n} r ,

    Pi = v_n^{-1} prod_{k<n} v_k^{1/2} (1-v_k)^{k/2-1} prod_{k=1}^n u_k^{-(k+1)/2},
    L_i = sum_{w<=i} lam_w / P_w,   P_w = prod_{l<w} u_l,
    M_j = sum_{k=j+1}^n P_k^2 c_k (L_k - L_j) + P_n u_n v_n lam_{n+1},   c_k = u_k (1-u_k) R_k,
    R_k = prod_{m=k}^{n-1} u_{m+1}/(1-v_m),  R_n = 1.

Polynomials: dict {exponent tuple: Fraction}; variables (u_1..u_n, U_1..U_n, v_1..v_n, V_1..V_n,
lam_1..lam_{n+1}) with U = 1-u, V = 1-v; exponents of u, U, v, V are stored doubled (half-integers).
Integration: int_0^1 x^a (1-x)^b = B(a+1, b+1) (k < n);  int_{Delta^2} u^a (1-u)^b v^c du dv
= B(a+1, b+c+2)/(c+1).  Values: q * pi^(j/2) with exact q.
"""
import sys
import math
from fractions import Fraction as Fr
from functools import lru_cache
import itertools


def gamma_half(x2):
    """Gamma(x2/2) for a positive integer x2, as (q, k): value = q * sqrt(pi)^k."""
    if x2 <= 0:
        raise ValueError("non-positive Gamma argument %s/2" % x2)
    if x2 % 2 == 0:
        return Fr(math.factorial(x2 // 2 - 1)), 0
    m = (x2 - 1) // 2            # Gamma(m + 1/2) = (2m)! sqrt(pi) / (4^m m!)
    return Fr(math.factorial(2 * m), 4 ** m * math.factorial(m)), 1


def beta2(a2, b2):
    """B(a2/2, b2/2) as (q, k)."""
    qa, ka = gamma_half(a2)
    qb, kb = gamma_half(b2)
    qc, kc = gamma_half(a2 + b2)
    return qa * qb / qc, ka + kb - kc


class Setup:
    def __init__(self, n):
        self.n = n
        self.nv = 4 * n + n + 1
        self.iu = lambda k: k - 1
        self.iU = lambda k: n + k - 1
        self.iv = lambda k: 2 * n + k - 1
        self.iV = lambda k: 3 * n + k - 1
        self.il = lambda w: 4 * n + w - 1

    def mono(self, **kw):
        e = [0] * self.nv
        for name, items in kw.items():
            idx = getattr(self, "i" + name)
            for k, p in items.items():
                e[idx(k)] += p
        return tuple(e)


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, 0) + va * vb
    return {k: v for k, v in out.items() if v != 0}


def padd(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def build(n, sources=None):
    """L_i, M_j as polynomials.  sources: None (symbolic lam_1..lam_{n+1}) or a dict
    {w: Fraction} (missing w: 0) -- then the lam exponents stay zero."""
    S = Setup(n)
    zero = tuple([0] * S.nv)

    def lamterm(w):
        if sources is None:
            return {S.mono(l={w: 1}): Fr(1)}
        c = Fr(sources.get(w, 0))
        return {zero: c} if c else {}

    # doubled exponents for u, U, v, V
    def P(w):        # prod_{l<w} u_l
        return {l: 2 for l in range(1, w)}

    L = {}
    acc = {}
    for i in range(1, n + 1):
        invP = {l: -2 for l in range(1, i)}
        acc = padd(acc, pmul(lamterm(i), {S.mono(u=invP): Fr(1)}))
        L[i] = dict(acc)

    def Lw(w):       # lam_w / P_w
        return pmul(lamterm(w), {S.mono(u={l: -2 for l in range(1, w)}): Fr(1)})

    M = {}
    for j in range(1, n + 1):
        tot = {}
        for k in range(j + 1, n + 1):
            # P_k^2 c_k = P_k^2 u_k U_k R_k
            eu = {l: 4 for l in range(1, k)}
            eu[k] = eu.get(k, 0) + 2
            eU = {k: 2}
            eV = {}
            for m in range(k, n):
                eu[m + 1] = eu.get(m + 1, 0) + 2
                eV[m] = eV.get(m, 0) - 2
            pref = {S.mono(u=eu, U=eU, V=eV): Fr(1)}
            diff = {}
            for w in range(j + 1, k + 1):
                diff = padd(diff, Lw(w))
            tot = padd(tot, pmul(pref, diff))
        eu = {l: 2 for l in range(1, n + 1)}
        tot = padd(tot, pmul({S.mono(u=eu, v={n: 2}): Fr(1)}, lamterm(n + 1)))
        M[j] = tot
    return S, L, M


def hafnian(n, L, M):
    one = {tuple([0] * len(next(iter(L[1])))): Fr(1)} if L[1] else None

    @lru_cache(maxsize=None)
    def haf(idx):
        if not idx:
            return ((tuple([0] * len(next(iter(M[n])))), Fr(1)),)
        i0, rest = idx[0], idx[1:]
        out = {}
        for k in range(len(rest)):
            j = rest[k]
            sub = dict(haf(rest[:k] + rest[k + 1:]))
            out = padd(out, pmul(pmul(L[i0], M[j]), sub))
        return tuple(out.items())

    return dict(haf(tuple(range(1, n + 1))))


def integrate(n, S, poly):
    """int 2^{-n/2} Pi * poly; returns {lam exponent tuple: {k: q}} (value sum q sqrt(pi)^k)."""
    # Pi, doubled exponents
    off = [0] * S.nv
    off[S.iv(n)] -= 2
    for k in range(1, n):
        off[S.iv(k)] += 1
        off[S.iV(k)] += k - 2
    for k in range(1, n + 1):
        off[S.iu(k)] -= k + 1
    res = {}
    minexp = None
    for key, c in poly.items():
        e = [x + y for x, y in zip(key, off)]
        q, kpi = Fr(c), 0
        for k in range(1, n):
            for (a, b) in ((e[S.iu(k)], e[S.iU(k)]), (e[S.iv(k)], e[S.iV(k)])):
                m = min(a, b)
                minexp = m if minexp is None else min(minexp, m)
                dq, dk = beta2(a + 2, b + 2)
                q *= dq
                kpi += dk
        a, b, cc, d = e[S.iu(n)], e[S.iU(n)], e[S.iv(n)], e[S.iV(n)]
        assert d == 0
        minexp = min(minexp if minexp is not None else 0, a, b, cc)
        dq, dk = beta2(a + 2, b + cc + 4)
        q *= dq / Fr(cc + 2, 2)
        kpi += dk
        lk = tuple(e[S.il(1):])
        slot = res.setdefault(lk, {})
        slot[kpi] = slot.get(kpi, 0) + q / 2 ** (n // 2)
    return res, minexp


def kernel_general(n):
    S, L, M = build(n)
    H = hafnian(n, L, M)
    res, minexp = integrate(n, S, H)
    return res, minexp, len(H)


if __name__ == "__main__":
    import sympy as sp
    for n in [int(x) for x in sys.argv[1:]] or [2, 4]:
        res, minexp, nterms = kernel_general(n)
        lam = sp.symbols("l1:%d" % (n + 2))
        expr = 0
        for lk, d in res.items():
            for k, q in d.items():
                if q:
                    expr += sp.Rational(q.numerator, q.denominator) * sp.sqrt(sp.pi) ** k * sp.Mul(*[x ** p for x, p in zip(lam, lk)])
        print("n=%d: %d terms, min doubled exponent %s (> -2 required)" % (n, nterms, minexp))
        print("  K =", sp.factor(sp.expand(expr)))
