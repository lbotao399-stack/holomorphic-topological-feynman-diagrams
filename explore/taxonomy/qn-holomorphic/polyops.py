"""Graded-commutative polynomials in jet letters (A, alpha), alpha in N^H, modulo a prime.

Letter (A, alpha) = d^alpha phi^A, parity pi[A].  Plane-wave convention: d^alpha phi -> u^alpha (no factorials).
Monomial: tuple of letters in increasing order (key (A, alpha)); odd letters at most once.
Poly: dict monomial -> residue mod P.
"""
from fractions import Fraction as Fr
import itertools
import math

P = 2147483647


def set_prime(p):
    global P
    P = p


def modq(x, p=None):
    p = p or P
    if isinstance(x, Fr):
        return x.numerator % p * pow(x.denominator % p, p - 2, p) % p
    return int(x) % p


def ratrec(a, p=None):
    p = p or P
    a %= p
    r0, r1, s0, s1 = p, a, 0, 1
    bound = int((p / 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def canon(word, pi):
    """sort a word of letters with Koszul sign; returns (sign, monomial) or (0, None)."""
    w = list(word)
    n = len(w)
    sign = 1
    # insertion sort tracking odd transpositions
    for i in range(1, n):
        j = i
        while j > 0 and w[j - 1] > w[j]:
            if pi[w[j - 1][0]] and pi[w[j][0]]:
                sign = -sign
            w[j - 1], w[j] = w[j], w[j - 1]
            j -= 1
    for i in range(1, n):
        if w[i] == w[i - 1] and pi[w[i][0]]:
            return 0, None
    return sign, tuple(w)


def padd(acc, mono, c):
    if c % P == 0:
        return
    v = (acc.get(mono, 0) + c) % P
    if v:
        acc[mono] = v
    else:
        acc.pop(mono, None)


def pscale(poly, c):
    c %= P
    return {m: v * c % P for m, v in poly.items() if v * c % P}


def psum(*polys):
    out = {}
    for p in polys:
        for m, v in p.items():
            padd(out, m, v)
    return out


def pmul_mono(m1, m2, pi):
    return canon(m1 + m2, pi)


def show(poly, names, H, maxterms=40, rat=True):
    items = sorted(poly.items())
    out = []
    for m, v in items[:maxterms]:
        c = ratrec(v) if rat else v
        lets = []
        for A, al in m:
            d = ''.join('d%d' % (a + 1) * k for a, k in enumerate(al) if k)
            lets.append((d + '.' if d else '') + names[A])
        out.append('%s %s' % (c, ' '.join(lets)))
    if len(items) > maxterms:
        out.append('... (%d terms)' % len(items))
    return '\n'.join(out)


def jets_upto(H, deg):
    """all multi-indices of total degree <= deg."""
    out = []
    for t in range(deg + 1):
        for c in itertools.combinations_with_replacement(range(H), t):
            al = [0] * H
            for a in c:
                al[a] += 1
            out.append(tuple(al))
    return out


def binom_split(alpha):
    """yield (beta, gamma, coefficient) with beta + gamma = alpha, coefficient = prod binom."""
    ranges = [range(a + 1) for a in alpha]
    for beta in itertools.product(*ranges):
        gamma = tuple(a - b for a, b in zip(alpha, beta))
        c = 1
        for a, b in zip(alpha, beta):
            c *= math.comb(a, b)
        yield beta, gamma, c
