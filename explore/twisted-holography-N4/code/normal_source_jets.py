"""Normal derivatives of the D3 source, tested on regular and Laurent polyvectors."""
from itertools import product
from math import comb, factorial

from single_particle import div_A0, div_A2, rank


def exponents(degree):
    for first in range(degree + 1):
        for second in range(degree + 1 - first):
            yield first, second, degree - first - second


def power(linear, degree):
    """Expand an ordinary power by multiplying linear forms in the normal frame."""
    terms = {(0, 0, 0): 1}
    for _ in range(degree):
        result = {}
        for monomial, coefficient in terms.items():
            for axis, amplitude in enumerate(linear):
                if amplitude:
                    raised = list(monomial)
                    raised[axis] += 1
                    key = tuple(raised)
                    result[key] = result.get(key, 0) + coefficient * amplitude
        terms = result
    return terms


def cutoff_germ(terms):
    """Divergence of the normal trivector near the deleted plane, in the volume frame."""
    germ = {}
    for monomial, coefficient in terms.items():
        for (_, index, slots), value in div_A0(((0, 0), monomial, (2, 3, 4))).items():
            axis, = set(range(3)) - {slot - 2 for slot in slots}
            # Contract the remaining bivector into dw1 dw2 dw3.
            germ[axis, index] = (-1)**axis * value * coefficient
    return germ


def laurent_measure(terms, index, axis):
    primitive = ((0, 0), index, (0, 1, axis + 2))
    return sum(value * terms.get(monomial, 0)
               for (_, monomial, _), value in div_A2(primitive).items())


def check():
    assert not cutoff_germ({(0, 0, 0): 1})
    for degree in (2, 4, 6, 8):
        rows = []
        for first, second, _ in exponents(degree):
            linear = (1, first, second)
            terms = power(linear, degree)
            germ = cutoff_germ(terms)
            lower = power(linear, degree - 1)
            for axis, index in product(range(3), exponents(degree - 1)):
                value = germ.get((axis, index), 0)
                assert value == degree * linear[axis] * lower.get(index, 0)
                assert laurent_measure(terms, index, axis) == -value
            rows.append(germ)
        assert rank(rows) == len(rows) == comb(degree + 2, 2)
        axial = {(0, 0, degree): 1}
        index = (0, 0, degree - 1)
        assert cutoff_germ(axial) == {(2, index): degree}
        assert -laurent_measure(axial, index, 2) * factorial(degree - 1) == factorial(degree)
        print(f'normal degree {degree}: {len(rows)} independent gradients, first jet order {degree-1}')


if __name__ == '__main__':
    check()
