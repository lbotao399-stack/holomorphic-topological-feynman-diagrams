"""Dimension of the jet-monomial space of the sector mu^9 gamma (g^1, K=0) of {I^12}_0:
9 distinct odd jets d^alpha mu^a, one even jet d^beta gamma_b, 23 derivatives, torus weight
#d_i + [b=i] - #[a=i] = 3 for every i (SL(5)-invariant densities times eps^3).
Generating function: e_9 over mu-jets (Newton: e_9 = sum_{lambda |- 9} eps_lambda p_lambda / z_lambda),
p_k = (sum_a t_a^{-k}) prod_i (1-(t_i q)^k)^{-1}; gamma series (sum_b t_b) prod_i (1-t_i q)^{-1}."""
from fractions import Fraction
from itertools import product
from math import factorial
from collections import Counter
from sympy.utilities.iterables import partitions

def series_coeffs(parts, nmax):
    c = [0] * (nmax + 1); c[0] = 1
    for k in list(parts) + [1]:
        for n in range(k, nmax + 1):
            c[n] += c[n - k]
    return c

def zlam(lam):
    z = 1
    for k, m in lam.items():
        z *= k ** m * factorial(m)
    return z

def count(nmu=9, nder=23, target=3):
    tot = Fraction(0)
    for lam in partitions(nmu):
        lam = dict(lam)
        parts = [k for k, m in lam.items() for _ in range(m)]
        eps = (-1) ** (nmu - len(parts))
        f = series_coeffs(parts, 60)
        # prefactor prod_{parts k} (sum_a t_a^{-k}) : multinomial over assignment of parts to directions
        pref = Counter({(0,) * 5: 1})
        for k in parts:
            new = Counter()
            for c, n in pref.items():
                for a in range(5):
                    cc = list(c); cc[a] += k; new[tuple(cc)] += n
            pref = new
        s = 0
        for c, n in pref.items():
            for b in range(5):
                ns = [target + c[i] - (1 if i == b else 0) for i in range(5)]
                if min(ns) < 0: continue
                prod_ = 1
                for x in ns: prod_ *= f[x]
                s += n * prod_
        tot += Fraction(eps * s, zlam(lam))
    return tot

if __name__ == '__main__':
    print('jet monomials mu^9 gamma, 23 derivatives, weight (3,3,3,3,3):', count())
    print('check small case mu^1 gamma^0 ... e_1 with 1 derivative weight (0,..): skipped')
