"""Own exact Dirichlet integration of a ratio density  g = prod_b b^{exps_b} * P(xi)  over a product of
simplices:  int_{Delta^m} prod_{j=0}^m b_j^{a_j} db_1..db_m = prod_j Gamma(a_j+1) / Gamma(m+1+sum_j a_j).
Values are kept as (Fraction q, k) meaning q * sqrt(pi)^k."""
from fractions import Fraction as Fr
import math


def gamma_half(x):
    """Gamma(x) for x > 0 integer or half-integer: (q, k) with Gamma = q sqrt(pi)^k"""
    x = Fr(x)
    if x <= 0:
        raise ValueError("divergent Dirichlet integral: Gamma(%s)" % x)
    if x.denominator == 1:
        return Fr(math.factorial(int(x) - 1)), 0
    assert x.denominator == 2
    j = int(x - Fr(1, 2))
    return Fr(math.factorial(2 * j), 4 ** j * math.factorial(j)), 1


def dirichlet(exps):
    q, k = Fr(1), 0
    for a in exps:
        g, gk = gamma_half(Fr(a) + 1)
        q *= g
        k += gk
    m = len(exps) - 1
    g, gk = gamma_half(m + 1 + sum(Fr(a) for a in exps))
    return q / g, k - gk


def integrate_engine_density(red, dens):
    """int g d xi of an htk RatioDensity with numeric sources, using only (num.terms(), exps, blocks).
    Returns {power of pi: Fraction}."""
    bidx = {b: i for i, b in enumerate(red.bary)}
    cidx = {c: i for i, c in enumerate(red.coords)}
    acc = {}
    for mon, c in dens.num.terms():
        q = Fr(int(c.numerator), int(c.denominator)) if hasattr(c, "numerator") else Fr(c)
        k = 0
        for blk in red.blocks:
            ex = [dens.exps[bidx[blk.bary[0]]]] + [dens.exps[bidx[b]] + mon[cidx[b]] for b in blk.bary[1:]]
            dq, dk = dirichlet(ex)
            q *= dq
            k += dk
        assert k % 2 == 0
        acc[k // 2] = acc.get(k // 2, Fr(0)) + q
    return {p: v for p, v in acc.items() if v != 0}
