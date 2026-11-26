"""Sparse modular solver for the period systems, with rational reconstruction.

rows: list of (dict col -> Fraction, dict label -> Fraction) meaning sum coef*x_col = const.
Forward elimination mod p (pivot = most complex column of the reduced row), back
substitution in reverse pivot order, CRT over primes, rational reconstruction, and an
exact check of every reconstructed value against a sample of the original rows.
"""
import random
import time
from fractions import Fraction as F

PRIMES = [2305843009213693951, 2305843009213693921, 2305843009213693907, 2305843009213693867,
          2305843009213693803, 2305843009213693771, 2305843009213693687, 2305843009213693669,
          2305843009213693657, 2305843009213693597, 2305843009213693591, 2305843009213693573]


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


PRIMES = [p for p in PRIMES if _is_prime(p)]


def ratrec(a, m):
    """rational reconstruction of a mod m (Wang)"""
    a %= m
    if a == 0:
        return F(0)
    bound = int((m // 2) ** 0.5)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return F(r1, s1) if s1 > 0 else F(-r1, -s1)


def solve_mod(rows, ncols, labels, p, rank_key):
    inv = lambda x: pow(x, p - 2, p)

    def red(fr):
        return fr.numerator % p * inv(fr.denominator % p) % p
    L = len(labels)
    lidx = {l: i for i, l in enumerate(labels)}
    pivots = {}
    order = []
    for row, const in rows:
        r = {}
        for c, v in row.items():
            x = red(v)
            if x:
                r[c] = x
        k = [0] * L
        for l, v in const.items():
            k[lidx[l]] = red(v)
        while True:
            hit = [c for c in r if c in pivots]
            if not hit:
                break
            for c in hit:
                if c not in r:
                    continue
                f = r[c]
                pr, pk = pivots[c]
                for cc, vv in pr.items():
                    nv = (r.get(cc, 0) - f * vv) % p
                    if nv:
                        r[cc] = nv
                    else:
                        r.pop(cc, None)
                for i in range(L):
                    if pk[i]:
                        k[i] = (k[i] - f * pk[i]) % p
        if not r:
            if any(k):
                raise ArithmeticError('inconsistent system mod p')
            continue
        c = max(r, key=rank_key)
        f = inv(r[c])
        r = {cc: vv * f % p for cc, vv in r.items()}
        k = [x * f % p for x in k]
        pivots[c] = (r, k)
        order.append(c)
    sol = {}
    for c in reversed(order):
        r, k = pivots[c]
        val = list(k)
        ok = True
        for cc, vv in r.items():
            if cc == c:
                continue
            if cc not in sol:
                ok = False
                break
            s = sol[cc]
            for i in range(L):
                if s[i]:
                    val[i] = (val[i] - vv * s[i]) % p
        if ok:
            sol[c] = val
    return sol, len(pivots)


def solve(rows, ncols, rank_key=lambda c: c, verbose=True, check=200):
    t = time.time()
    labels = sorted({l for r, k in rows for l in k} | {'1'})
    M = 1
    acc = None
    prev = None
    for np_, p in enumerate(PRIMES):
        sol, rank = solve_mod(rows, ncols, labels, p, rank_key)
        if acc is None:
            acc = {c: list(v) for c, v in sol.items()}
            M = p
        else:
            # CRT combine
            minv = pow(M, -1, p)
            for c, v in sol.items():
                a = acc.get(c)
                if a is None:
                    continue
                for i in range(len(labels)):
                    a[i] = a[i] + M * ((v[i] - a[i]) * minv % p)
            M *= p
        rec = {}
        good = True
        for c, v in acc.items():
            out = {}
            for i, l in enumerate(labels):
                x = ratrec(v[i], M)
                if x is None:
                    good = False
                    break
                if x:
                    out[l] = x
            if not good:
                break
            rec[c] = out
        if verbose:
            print('  prime %d: rank %d, solved %d, reconstructed %s  %.1fs' % (np_ + 1, rank, len(sol), good,
                                                                               time.time() - t), flush=True)
        if good and prev == rec:
            break
        prev = rec if good else None
    else:
        raise ArithmeticError('rational reconstruction did not stabilize')
    # exact check on a sample of rows
    rng = random.Random(1)
    sample = rows if len(rows) <= check else rng.sample(rows, check)
    for row, const in sample:
        if not all(c in rec for c in row):
            continue
        lhs = {}
        for c, v in row.items():
            for l, x in rec[c].items():
                lhs[l] = lhs.get(l, 0) + v * x
        lhs = {l: x for l, x in lhs.items() if x}
        rhs = {l: x for l, x in const.items() if x}
        if lhs != rhs:
            raise ArithmeticError('reconstruction fails an original relation')
    return rec, rank
