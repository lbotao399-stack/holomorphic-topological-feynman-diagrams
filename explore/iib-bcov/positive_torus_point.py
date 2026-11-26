"""The top bare primary loop on the positive-degree compactification colors of T^4."""
from fractions import Fraction
from itertools import permutations
from math import factorial

from compactified import extremal, torus_product
from scalar_fan import independent_integral, witness


def polarized_root(colors, external):
    """Unnormalized product of odd point letters, in the displayed order."""
    total = 0
    for order in permutations(range(len(colors))):
        # Koszul sign of the odd evaluation letters in this Wick assignment.
        sign = (-1)**sum(i > j for at, i in enumerate(order) for j in order[at+1:])
        value = 0  # Unit of the exterior algebra, in bit-mask coordinates.
        for i in order:
            factor, value = torus_product(value, colors[i])
            sign *= factor
        factor, value = torus_product(value, external)
        if value == 15:
            total += sign*factor
    return total


def classical_quartic():
    """Absolute classical variation, with three spectator three-forms and one colored vector."""
    one_forms = (1,2,4,8)
    dual = [(torus_product(a, 15^a)[0], 15^a) for a in one_forms]
    scalar_colors = (*dual[:3], (1,0))
    scalar_momenta = ((-1,-1,0),(0,0,0),(0,0,0),(1,0,0))
    vector_color = dual[3]
    assert tuple(map(sum, zip(*scalar_momenta,(0,1,0)))) == (0,0,0)
    terms = []
    for differentiated, (scalar_sign, color) in enumerate(scalar_colors):
        sign, varied = torus_product(vector_color[1], color)
        sign *= scalar_sign*vector_color[0]*scalar_momenta[differentiated][0]
        if not sign:
            continue
        colors = list(scalar_colors)
        colors[differentiated] = (sign, varied)
        for order in permutations(range(4)):
            value = 1
            for a,j in zip(one_forms, order):
                s, result = torus_product(a, colors[j][1])
                value *= s*colors[j][0]*(result == 15)
            if value:
                terms.append(value)
    assert len(terms) == 1
    return abs(terms[0])


if __name__ == '__main__':
    assert classical_quartic() == 1
    for n in (1,2,3):
        assert witness(n, 1) > 0
        geometric = independent_integral(n, 1)
        assert geometric == extremal(n)
        colors = [1 << i for i in range(n+1)]
        complement = 15 ^ sum(colors)
        assert polarized_root(colors, complement) == factorial(n+1)
        assert polarized_root([colors[0], *colors[:-1]], complement) == 0
        print('Genus, normalized fan, polarized product:', n, geometric,
              geometric*polarized_root(colors, complement))
    assert geometric*polarized_root([1,2,4,8], 0) == Fraction(1,10080)
