"""Mixed order-two brackets in the manuscript's unnormalized alternating convention.

The bracket and the CE differential use the project's formal trace algebra with N
indeterminate. The Cartan restriction is evaluated independently on diagonal waves.
Run: python3 notes/interpretation/mixed_brackets.py
"""
from fractions import Fraction as F
from itertools import product
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'deepen' / 'deform' / 'code'))
import ce


def generator(component, colors):
    # Taylor expansion of B(partial_i A,A):  (p_0^i-p_1^i) B(X_0,X_1).
    monomial = ce.monomial_full(2, {(0, component): 1})
    return ce.alt_raw(ce.raw_from({monomial: 1}, colors), 2)


def cubic(colors):
    raw = ce.raw_from(ce.bracket_poly(0, 1, 3), colors)
    return ce.scale(ce.alt_raw(raw, 3), F(1, 6))


def mixed_expected(component):
    # On sum p=0, alt(p_0^i [p_0 p_1] tr0 tr12) =
    # 2 [p_0 p_1] sum_cyclic p_j^i tr(j)tr(kl).
    p0 = {ce.monomial_full(3, {(0, component): 1}): 1}
    polynomial = ce.poly_mul(p0, ce.bracket_poly(0, 1, 3))
    return ce.alt_raw(ce.raw_from(polynomial, {((0,), (1, 2)): {0: 1}}), 3)


def diagonal_value(cochain, matrices, momenta):
    n = len(matrices[0])
    independent_momenta = [x for p in momenta[:-1] for x in p]
    result = 0
    for (colors, powers), coefficient in cochain.items():
        color = prod(sum(prod(matrices[i][j] for i in cycle) for j in range(n))
                     for cycle in colors)
        momentum = prod(x**d for x, d in zip(independent_momenta, powers))
        result += color*momentum*ce.np_eval(coefficient, n)
    return result


def main():
    # The invariant forms corresponding to the identity and the central projection.
    identity = [generator(i, {((0, 1),): {0: 1}}) for i in range(2)]
    central = [generator(i, {((0,), (1,)): {-1: 1}}) for i in range(2)]
    linear = {(((0,),), ()): {0: 1}}
    d = cubic({((0, 1, 2),): {0: 1}, ((0, 2, 1),): {0: 1}})
    s1 = cubic({((0,), (1, 2)): {0: 1}, ((0, 2), (1,)): {0: 1},
                ((0, 1), (2,)): {0: 1}})
    s0 = cubic({((0,), (1,), (2,)): {0: 1}})
    m1 = ce.lin((1, s1), ({-1: -3}, s0))
    traceless = ce.lin((1, d), ({-1: -2}, s1), ({-2: 4}, s0))

    for arity, cochain in [(1, linear)] + [(2, g) for g in identity+central] + [
            (3, j) for j in (d, s1, s0, traceless, m1)]:
        assert not ce.d_ce(cochain, arity)
        if arity > 1:
            assert not ce.bracket(linear, 1, cochain, arity)
    # Both T vertices carry dbar z_1 dbar z_2, whose exterior square vanishes.

    for left, right in product(identity+central, repeat=2):
        assert not ce.bracket(left, 2, right, 2)
    for i in range(2):
        for j in (traceless, m1, s0, d, s1):
            assert not ce.bracket(identity[i], 2, j, 3)
        for j in (traceless, s0):
            assert not ce.bracket(central[i], 2, j, 3)
        mixed = ce.bracket(central[i], 2, m1, 3)
        assert mixed == mixed_expected(i)
        assert ce.bracket(m1, 3, central[i], 2) == ce.scale(mixed, -1)
        assert not ce.d_ce(mixed, 3)
        # On gl_2, D = S1-S0. The coefficient of M1 in x0*D+x1*S1 is x0+x1.
        assert ce.at_N(ce.bracket(central[i], 2, d, 3), 2) == ce.at_N(mixed, 2)
        assert ce.bracket(central[i], 2, s1, 3) == mixed
        for n in (1, 2, 3, 5):
            momenta = ((1, 1), (-1, 0), (0, -1))
            if n == 1:
                assert diagonal_value(mixed, ((2,), (-3,), (5,)), momenta) == 0
                continue
            identity_matrix = [1]*n
            h = [1, -1]+[0]*(n-2)
            arguments = (identity_matrix, h, h)
            assert diagonal_value(m1, arguments, momenta) == 2*n
            assert diagonal_value(mixed, arguments, momenta) == 4*n
    print('Both components closed. Cartan evaluation 4N for N=2,3,5; gl_1 vanishes.')
    print('Mixed class vanishes iff chi*(b_c-b_s)=0 for N>=2.')
    print('gl_2: chi=x0+x1; top condition x1*(x0+x1)=0.')


if __name__ == '__main__':
    main()
