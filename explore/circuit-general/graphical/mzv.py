"""Exact multiple zeta values up to weight WMAX, reduced to a fixed algebra basis.

zeta(s1,...,sd) = sum_{n1>...>nd>=1} prod n_i^{-s_i},  s1 >= 2.
Iterated-integral word: x0^{s1-1} x1 ... x0^{sd-1} x1  (letters '0','1'); zeta(u) zeta(v) = zeta(u sh v).
Relations: double shuffle  sh(u,v) = st(s,t)  and Hoffman  sh('1', u) = st((1), s)  (divergent terms cancel).
Basis (irreducibles): z2, z3, z5, z7, z9, z11, z53 = zeta(5,3), z73 = zeta(7,3), z533 = zeta(5,3,3).
An MZV element is a dict {monomial: Fraction}, monomial = sorted tuple of irreducible names.
"""
import json
import os
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
WMAX = 11
IRRED = {'z2': (2,), 'z3': (3,), 'z5': (5,), 'z7': (7,), 'z9': (9,), 'z11': (11,),
         'z53': (5, 3), 'z73': (7, 3), 'z533': (5, 3, 3)}
IWEIGHT = {k: sum(v) for k, v in IRRED.items()}
DIMS = [1, 0, 1, 1, 1, 2, 2, 3, 4, 5, 7, 9]


def word(s):
    return ''.join('0' * (x - 1) + '1' for x in s)


def comp(w):
    out, c = [], 0
    for a in w:
        c += 1
        if a == '1':
            out.append(c)
            c = 0
    assert c == 0
    return tuple(out)


@lru_cache(maxsize=None)
def shuffle(u, v):
    if not u:
        return {v: 1}
    if not v:
        return {u: 1}
    out = {}
    for w, c in shuffle(u[1:], v).items():
        out[u[0] + w] = out.get(u[0] + w, 0) + c
    for w, c in shuffle(u, v[1:]).items():
        out[v[0] + w] = out.get(v[0] + w, 0) + c
    return out


@lru_cache(maxsize=None)
def stuffle(s, t):
    if not s:
        return {t: 1}
    if not t:
        return {s: 1}
    out = {}
    for r, c in stuffle(s[1:], t).items():
        k = (s[0],) + r
        out[k] = out.get(k, 0) + c
    for r, c in stuffle(s, t[1:]).items():
        k = (t[0],) + r
        out[k] = out.get(k, 0) + c
    for r, c in stuffle(s[1:], t[1:]).items():
        k = (s[0] + t[0],) + r
        out[k] = out.get(k, 0) + c
    return out


def compositions(w):
    """convergent compositions of weight w (s1 >= 2)"""
    out = []

    def rec(rem, pre):
        if rem == 0:
            if pre and pre[0] >= 2:
                out.append(tuple(pre))
            return
        for x in range(1, rem + 1):
            rec(rem - x, pre + [x])
    rec(w, [])
    return out


def monomials(w):
    names = sorted(IRRED, key=lambda k: (IWEIGHT[k], k))
    out = []

    def rec(i, rem, pre):
        if rem == 0:
            out.append(tuple(pre))
            return
        for j in range(i, len(names)):
            n = names[j]
            if IWEIGHT[n] <= rem:
                rec(j, rem - IWEIGHT[n], pre + [n])
    rec(0, w, [])
    return out


def _reduce_weight(w):
    """return {composition: {monomial: Fraction}} for all convergent compositions of weight w"""
    comps = compositions(w)
    if not comps:
        return {}
    # column order: eliminate long/complicated compositions first, keep short/simple ones free
    order = sorted(comps, key=lambda s: (-len(s), [-x for x in s]))
    col = {s: i for i, s in enumerate(order)}
    rows = []
    for a in range(2, w - 1):
        for s in compositions(a):
            for t in compositions(w - a):
                if (a, s) > (w - a, t):
                    continue
                r = {}
                for x, c in shuffle(word(s), word(t)).items():
                    k = comp(x)
                    r[k] = r.get(k, 0) + c
                for x, c in stuffle(s, t).items():
                    r[x] = r.get(x, 0) - c
                rows.append(r)
    for s in compositions(w - 1):
        r = {}
        for x, c in shuffle('1', word(s)).items():
            k = comp(x)
            r[k] = r.get(k, 0) + c
        for x, c in stuffle((1,), s).items():
            r[x] = r.get(x, 0) - c
        bad = [k for k, c in r.items() if c and k[0] == 1]
        assert not bad, bad
        rows.append(r)
    # sparse Gauss-Jordan
    piv = {}  # column index -> row dict (normalized, pivot coefficient 1)
    for r in rows:
        r = {col[k]: F(c) for k, c in r.items() if c}
        for j in sorted(list(r)):
            if j in r and j in piv:
                c = r.pop(j)
                for jj, cc in piv[j].items():
                    if jj != j:
                        r[jj] = r.get(jj, 0) - c * cc
                        if r[jj] == 0:
                            del r[jj]
        if not r:
            continue
        p = min(r)
        c0 = r[p]
        r = {j: c / c0 for j, c in r.items()}
        for j2, r2 in piv.items():
            if p in r2:
                c = r2.pop(p)
                for jj, cc in r.items():
                    if jj != p:
                        r2[jj] = r2.get(jj, 0) - c * cc
                        if r2[jj] == 0:
                            del r2[jj]
        piv[p] = r
    free = [j for j in range(len(order)) if j not in piv]
    assert len(free) == DIMS[w], (w, len(free), DIMS[w])

    def tofree(vec):  # vec over compositions -> over free columns
        out = {}
        for k, c in vec.items():
            j = col[k]
            if j in piv:
                for jj, cc in piv[j].items():
                    if jj != j:
                        out[jj] = out.get(jj, 0) - c * cc
            else:
                out[j] = out.get(j, 0) + c
        return {j: c for j, c in out.items() if c}
    mons = monomials(w)
    assert len(mons) == DIMS[w]
    # monomial -> composition vector (stuffle product of irreducible compositions)
    mat = []
    for m in mons:
        vec = {(): 1}
        for n in m:
            new = {}
            for s, c in vec.items():
                for x, cc in stuffle(s, IRRED[n]).items():
                    new[x] = new.get(x, 0) + c * cc
            vec = new
        mat.append(tofree(vec))
    # solve: free column j = sum_m A[j][m] * monomial m ; invert matrix M[m][j]
    n = len(mons)
    fidx = {j: i for i, j in enumerate(free)}
    M = [[F(0)] * n + [F(int(i == k)) for k in range(n)] for i in range(n)]
    for i, v in enumerate(mat):
        for j, c in v.items():
            M[i][fidx[j]] = F(c)
    # we have mon_i = sum_j M[i][j] free_j ; invert to free_j = sum_i Inv[j][i] mon_i
    A = [row[:] for row in M]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    inv = [row[n:] for row in A]  # inv[j][i]: free_j in terms of mon_i ... careful with transpose
    # A started as [M | I]; after reduction [I | M^{-1}] ; mon = M free  =>  free = M^{-1} mon
    # M here is indexed M[i=mon][j=free]: mon_i = sum_j M[i][j] free_j, i.e. mon = M free (matrix-vector)
    out = {}
    for s in comps:
        v = tofree({s: 1})
        res = {}
        for j, c in v.items():
            jj = fidx[j]
            for i in range(n):
                x = inv[jj][i]
                if x:
                    res[mons[i]] = res.get(mons[i], 0) + c * x
        out[s] = {m: x for m, x in res.items() if x}
    return out


_TABLE = {}


def table(w):
    if w in _TABLE:
        return _TABLE[w]
    path = os.path.join(HERE, 'data', 'mzv_w%d.json' % w)
    if os.path.exists(path):
        raw = json.load(open(path))
        tab = {tuple(json.loads(k)): {tuple(m.split('*')) if m else (): F(v) for m, v in val.items()}
               for k, val in raw.items()}
    else:
        tab = _reduce_weight(w)
        os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
        json.dump({json.dumps(list(k)): {'*'.join(m): str(v) for m, v in val.items()} for k, val in tab.items()},
                  open(path, 'w'))
    _TABLE[w] = tab
    return tab


def zeta(*s):
    """reduced value of the convergent MZV zeta(s)"""
    s = tuple(s)
    if not s:
        return {(): F(1)}
    return dict(table(sum(s))[s])


# ---- MZV element arithmetic (dict monomial -> Fraction) ----
def add(a, b, c=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + c * v
        if out[k] == 0:
            del out[k]
    return out


def mul(a, b):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = tuple(sorted(k1 + k2))
            out[k] = out.get(k, 0) + v1 * v2
    return {k: v for k, v in out.items() if v}


def scal(a, c):
    return {k: v * c for k, v in a.items() if v * c}


def weight(m):
    return sum(IWEIGHT[n] for n in m)


def tostr(a):
    if not a:
        return '0'
    parts = []
    for m in sorted(a, key=lambda m: (weight(m), m)):
        v = a[m]
        parts.append('%s%s' % (v, ('*' + '*'.join(m)) if m else ''))
    return ' + '.join(parts)


# ---- numerics (for checks) ----
def numeric(a, dps=30):
    import mpmath as mp
    mp.mp.dps = dps + 10
    vals = {'z2': mp.zeta(2), 'z3': mp.zeta(3), 'z5': mp.zeta(5), 'z7': mp.zeta(7), 'z9': mp.zeta(9),
            'z11': mp.zeta(11)}
    tot = mp.mpf(0)
    for m, v in a.items():
        t = mp.mpf(v.numerator) / v.denominator
        for n in m:
            if n not in vals:
                vals[n] = mzv_numeric(IRRED[n], dps)
            t *= vals[n]
        tot += t
    return tot


def mzv_numeric(s, dps=30):
    """zeta(s) by nested partial sums with Euler-Maclaurin-free tail: uses the
    integral representation via multiple polylog series at 1/2 (Hoelder convolution)."""
    import mpmath as mp
    mp.mp.dps = dps + 10
    w = word(s)  # integral from 0 to 1 of omega word read left (outer, near 1) to right (inner, near 0)
    # zeta(w) = sum_{k=0}^{n} I(1/2 -> 1; w[:k]) * I(0 -> 1/2; w[k:])
    # I(0->x; a1..an) with a_i in {0,1}: omega_0 = dt/t, omega_1 = dt/(1-t), outer letter first
    N = int((dps + 10) * 3.33) + 20

    def series(wd, x):
        # returns value of I(0->x; wd) using power series (wd must end with '1')
        # represent iterated integral as power series sum c_n x^n, built from innermost letter
        if not wd:
            return mp.mpf(1)
        c = [mp.mpf(0)] * (N + 1)  # coefficient of x^n
        # innermost: int_0^x dt/(1-t) = sum_{n>=1} x^n / n
        assert wd[-1] == '1'
        for n in range(1, N + 1):
            c[n] = mp.mpf(1) / n
        for a in reversed(wd[:-1]):
            if a == '0':
                c = [c[0]] + [c[n] / n for n in range(1, N + 1)]
                assert c[0] == 0
            else:
                # int_0^x f(t) dt/(1-t): coefficient: sum_{m<n} c_m / n
                new = [mp.mpf(0)] * (N + 1)
                acc = mp.mpf(0)
                for n in range(1, N + 1):
                    acc += c[n - 1]
                    new[n] = acc / n
                c = new
        return mp.fsum(c[n] * x ** n for n in range(1, N + 1))

    def dual(wd):
        # I(1/2 -> 1; wd) = I(0 -> 1/2; dual) with t -> 1-t: omega_0 <-> omega_1, order reversed
        return ''.join('1' if a == '0' else '0' for a in reversed(wd))
    tot = mp.mpf(0)
    half = mp.mpf(1) / 2
    for k in range(len(w) + 1):
        outer, inner = w[:k], w[k:]
        # outer part integrates from 1/2 to 1; inner from 0 to 1/2
        if inner and inner[-1] != '1':
            continue
        if outer and outer[0] == '1':
            pass
        vo = series(dual(outer), half) if outer else mp.mpf(1)
        vi = series(inner, half) if inner else mp.mpf(1)
        tot += vo * vi
    return tot


if __name__ == '__main__':
    import sys
    import time
    import mpmath as mp
    top = int(sys.argv[1]) if len(sys.argv) > 1 else WMAX
    for w in range(2, top + 1):
        t = time.time()
        tab = table(w)
        print('weight %d: %d MZVs reduced (%.1fs)' % (w, len(tab), time.time() - t), flush=True)
    # numerical checks of random reductions
    mp.mp.dps = 40
    for s in [(2, 1), (3, 1, 1), (4, 2), (3, 3, 2), (5, 2, 1), (2, 2, 2, 2), (4, 3, 1), (6, 2, 1), (3, 2, 2, 2),
              (5, 1, 1, 1, 1), (4, 4, 2), (2, 3, 3, 2), (7, 1, 1, 1, 1), (2, 2, 2, 2, 3)]:
        if sum(s) > top:
            continue
        a = mzv_numeric(s, 30)
        b = numeric(zeta(*s), 30)
        print('zeta%s = %s   |num - red| = %s' % (s, tostr(zeta(*s)), mp.nstr(abs(a - b), 3)))
