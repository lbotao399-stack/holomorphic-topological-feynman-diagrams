"""Claim b3: number of jet monomials (9 distinct odd jets d^alpha mu^a, one jet d^beta gamma_b,
torus weight sum_mu (alpha - e_a) + beta + e_b = (3,3,3,3,3); then sum|alpha| + |beta| = 15 + 9 - 1 = 23).

Own generating function, organized by n_a = #{mu jets with index a}:
   count = sum_b sum_{n_1+..+n_5 = N} [t^{c + n - e_b}] prod_a E_{n_a}(t) * prod_i (1 - t_i)^{-1},
   E_n(t) = e_n(all monomials t^alpha, alpha in N^5) = sum_{lambda |- n} eps_lambda / z_lambda prod_{k in lambda} prod_i (1 - t_i^k)^{-1},
so every coefficient factorizes over the five coordinates.  Checked against explicit enumeration.
"""
import sys, os, time
from fractions import Fraction as Fr
from itertools import product, combinations
from functools import lru_cache
from math import factorial
from collections import Counter
from sympy.utilities.iterables import partitions as sym_partitions

D5 = 5


def parts_list(n):
    out = []
    for p in sym_partitions(n):
        p = dict(p)
        lst = tuple(sorted(k for k, m in p.items() for _ in range(m)))
        z = 1
        for k, m in p.items():
            z *= k ** m * factorial(m)
        eps = (-1) ** (n - len(lst))
        out.append((lst, Fr(eps, z)))
    return out


@lru_cache(maxsize=None)
def coef1(parts, c, with_gamma):
    """[t^c] prod_{k in parts} (1 - t^k)^{-1} * (1 - t)^{-with_gamma}."""
    if c < 0:
        return 0
    a = [0] * (c + 1); a[0] = 1
    ks = list(parts) + ([1] if with_gamma else [])
    for k in ks:
        for n in range(k, c + 1):
            a[n] += a[n - k]
    return a[c]


def compositions(N, k):
    if k == 1:
        yield (N,)
        return
    for i in range(N + 1):
        for rest in compositions(N - i, k - 1):
            yield (i,) + rest


def count_gf(nmu, target, ngamma=1):
    PL = {n: parts_list(n) for n in range(nmu + 1)}
    PL[0] = [((), Fr(1))]
    tot = Fr(0)
    bs = range(D5) if ngamma else [None]
    for n in compositions(nmu, D5):
        for b in bs:
            c = [target[i] + n[i] - (1 if b == i else 0) for i in range(D5)]
            if min(c) < 0:
                continue
            for choice in product(*[PL[n[a]] for a in range(D5)]):
                parts = tuple(sorted(sum((list(pl) for pl, _ in choice), [])))
                w = Fr(1)
                for _, cz in choice:
                    w *= cz
                val = 1
                for i in range(D5):
                    val *= coef1(parts, c[i], bool(ngamma))
                    if val == 0:
                        break
                tot += w * val
    assert tot.denominator == 1
    return int(tot)


def jets_upto(maxdeg):
    J = []
    for d in range(maxdeg + 1):
        for alpha in product(range(d + 1), repeat=D5):
            if sum(alpha) == d:
                for a in range(D5):
                    J.append((d, tuple(alpha[i] - (1 if i == a else 0) for i in range(D5))))
    return J   # sorted by degree


def count_brute(nmu, nder, target):
    """Explicit enumeration: increasing index tuples of distinct mu jets with total degree <= nder, then the
    gamma jets d^beta gamma_b with beta + e_b = target - sum_mu weights."""
    J = jets_upto(nder)
    tot = 0
    def rec(start, k, deg, wsum):
        nonlocal tot
        if k == nmu:
            c = [target[i] - wsum[i] for i in range(D5)]
            if sum(c) != (nder - deg) + 1:
                return
            for b in range(D5):
                if all((c[i] - (1 if i == b else 0)) >= 0 for i in range(D5)):
                    tot += 1
            return
        for j in range(start, len(J)):
            d, w = J[j]
            if deg + d * (nmu - k) > nder and d > 0:
                # jets are sorted by degree: remaining picks have degree >= d
                break
            rec(j + 1, k + 1, deg + d, [wsum[i] + w[i] for i in range(D5)])
    rec(0, 0, 0, [0] * D5)
    return tot


if __name__ == '__main__':
    t0 = time.time()
    cases = [(1, 5, (1,) * 5), (2, 6, (1,) * 5), (4, 3, (0,) * 5), (3, 2, (0,) * 5),
             (3, 5, (1, 1, 1, 0, 0)), (5, 6, (1, 1, 0, 0, 0)), (5, 4, (0, 0, 0, 0, 0)),
             (9, 6, (0, 0, 0, -1, -1)), (9, 7, (1, 0, 0, -1, -1)), (6, 6, (2, 0, -1, 0, 0))]
    for nmu, nder, tg in cases:
        assert sum(tg) == nder + 1 - nmu
        b = count_brute(nmu, nder, tg)
        # gf with negative target components: shift handled since c = target + n - e_b must be >= 0
        g = count_gf(nmu, tg)
        print('n_mu=%d n_der=%d target=%s: brute %d, generating function %d, equal %s (%.1f s)'
              % (nmu, nder, tg, b, g, b == g, time.time() - t0))
    N = count_gf(9, (3,) * 5)
    print('mu^9 gamma, weight (3,3,3,3,3) [23 derivatives]:', N, '; explorer 63414968175297950:', N == 63414968175297950)
    print('sec', round(time.time() - t0, 1))
