"""Exact finite marked-field changes on the flat fivefold quadratic profile."""

from itertools import product
from math import comb, factorial, prod

import sympy as sp


def transverse_derivative(genus, metric, operator):
    """Apply the diagonal second-order operator g times to B^g at zero."""
    total = 0
    for exponents in product(range(genus + 1), repeat=3):
        last = genus - sum(exponents)
        if last < 0:
            continue
        total += prod(
            comb(2 * k, k) * (h * weight) ** k
            for k, h, weight in zip((*exponents, last), metric, operator)
        )
    return sp.Rational(factorial(genus) ** 2, 2 ** genus) * total


def check():
    # The correlated profile at level j has coefficient a_k (-B)^(k-j)/(k-j)!.
    # Its classical implicit scalar u equals B+s to first order in the wave s.
    background, wave, aux = sp.symbols("B s u")
    strengths = {2: 3, 3: -2, 4: 5}
    tree = background + wave + sum(
        a * (-background) ** (k-j) * aux**j / (factorial(k-j) * factorial(j))
        for k, a in strengths.items() for j in range(k+1)
    )
    assert sp.Poly((tree-aux).subs(aux, background+wave).expand(), wave).nth(1) == 0
    assert sp.expand(tree.subs({aux: background, wave: 0})-background) == 0

    # Change affine volume coordinates before constructing the Cartesian Hessian.
    frame = sp.eye(5)
    for i, j, value in ((0, 1, 2), (1, 3, -1), (3, 0, 3), (2, 4, 2), (4, 0, -1)):
        frame[i, :] += value * frame[j, :]
    assert frame.det() == 1
    tangent = frame[0, :]
    for metric in ((1, 1, 1, 1), (1, -2, 1, 3)):
        original = frame.T * sp.diag(0, *metric) * frame
        for momentum in (-2, 3):
            variation = momentum**3 * tangent.T * tangent
            for genus in range(1, 7):
                # A rank-one wave makes each 4-by-4 cofactor affine. The exact
                # test polynomial has degree at most g+1 in the wave parameter.
                values = []
                for sample in range(genus+2):
                    hessian = original + sample * variation
                    if genus == 1:
                        values.append(hessian.det())
                        continue
                    operator = frame * hessian.adjugate() * frame.T
                    assert operator.is_diagonal()
                    value = sample * momentum**(2*genus+1) * operator[0, 0]**genus
                    if genus in strengths:
                        # Only (-B)^g has total transverse jet order 2g.
                        normal = [operator[j, j] for j in range(1, 5)]
                        value += (sp.Rational(strengths[genus] * (-1)**genus, factorial(genus))
                                  * transverse_derivative(genus, metric, normal))
                    values.append(value)
                # Newton forward formula differentiated at zero.
                first = 0
                for degree in range(1, len(values)):
                    values = [right-left for left, right in zip(values, values[1:])]
                    first += sp.Rational((-1)**(degree-1), degree) * values[0]
                assert first == prod(metric)**genus * momentum**(2*genus+1)
                for m2, m3, m4 in product(range(2*genus+1), repeat=3):
                    if m2 + 2*m3 + 3*m4 == 2*genus:
                        slots = 4*genus+1 + m2+m3+m4
                        level = 2*m2 + 3*m3 + 4*m4
                        assert 2*slots - 2*level - 4*genus == 2
    print("Tree profile, marked degree and polarized Hessian changes pass through genus six.")


if __name__ == "__main__":
    check()
