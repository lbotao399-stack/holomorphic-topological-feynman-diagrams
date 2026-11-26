"""Ring contraction and the rooted box on scalar components paired with H^2(Y)."""
from fractions import Fraction
from itertools import permutations, product
from math import comb

from compactification import cubic, k3_product, torus_product


def integrated_green(i, j):
    """Integral of min(T_i,T_j)-T_i*T_j on the standard four-edge simplex."""
    linear = Fraction(min(i, j), 24)
    quadratic = sum(Fraction(1 + (a == b), 120)
                    for a in range(i) for b in range(j))
    return linear - quadratic


def rooted_box(momenta, v):
    """Sum all orders around the rooted box, with vector label 2."""
    determinant = sum(
        (-1)**sum(order[i] > order[j] for i in range(3) for j in range(i+1, 3))
        * momenta[0][order[0]] * momenta[1][order[1]] * momenta[2][order[2]]
        for order in permutations(range(3))
    )
    charges = [sum(x*y for x, y in zip(p, v)) for p in momenta]
    assert charges[2] == 0
    coefficient = sum(integrated_green(i, order.index(2)+1)*charges[label]
                      for order in permutations(range(3))
                      for i, label in enumerate(order, 1))
    result = determinant**2 * coefficient / 2
    assert result == Fraction(determinant**2 * sum(charges), 24)
    return result


def check_middle_ring(multiply, top, middle, parity):
    """Trace all tree contractions for two or three cut middle labels."""
    labels = range(top+1)
    dual = {}
    for i in labels:
        j, = [j for j in labels if cubic(multiply, top, i, j, 0)]
        dual[i] = (j, cubic(multiply, top, i, j, 0))
    handle = [0]*(top+1)
    for i, (j, metric) in dual.items():
        sign, output = multiply(i, j)
        handle[output] += (-1)**parity(i)*metric*sign
    assert not any(handle[:-1])
    for a in middle:
        partners = []
        for b in middle:
            sign, channel = multiply(a, b)
            if sign:
                assert channel == top
                partners.append(b)
            for external in labels:
                tree = sum(cubic(multiply, top, a, b, i)*metric
                           * cubic(multiply, top, j, external, 0)
                           for i, (j, metric) in dual.items())
                assert tree == sign * (external == 0)
            for c in middle:
                assert sign*multiply(channel, c)[0] == 0
        assert len(partners) == 1
        assert handle[top]*multiply(top, a)[0] == 0
    return len(middle)


if __name__ == '__main__':
    basis = [(1,0,0), (0,1,0), (0,0,1)]
    assert rooted_box(basis, (1,1,0)) == Fraction(1,12)
    for v in ((2,-1,0), (0,3,0), (1,-2,0)):
        for a, b in product((-2,1,3), repeat=2):
            rooted_box(((a,2,1), (1,b,-1), basis[2]), v)
    for multiply, top, middle, parity, expected in (
        (k3_product, 23, range(1,23), lambda i: 0, 252),
        (torus_product, 15, [i for i in range(16) if i.bit_count() == 2],
         lambda i: i.bit_count() % 2, 20),
    ):
        dimension = check_middle_ring(multiply, top, middle, parity)
        assert comb(dimension+1, 2)-1 == expected
        print('Middle rank, one-loop quadratic kernel, unit-wave coefficient:',
              dimension, expected, Fraction(1,12))
