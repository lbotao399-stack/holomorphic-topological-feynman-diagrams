"""Leading Laurent coefficients of the five-line K3,3 projector."""
from fractions import Fraction as F
from functools import cache
from math import comb, prod
from flint import fmpq

from quartic_five_cut import FREE, RETAINED, cut_images, entries, gram_jet, gram_parts
from quartic_planar_double_pole import multiply, pochhammer


@cache
def beta(a, b, i, j):
    # Each argument has the form (integer at ep=0, coefficient of ep).
    return multiply([(pochhammer(*a, i), 1), (pochhammer(*b, j), 1),
                     (pochhammer(a[0]+b[0], a[1]+b[1], i+j), -1)])


def coefficients(factor, lowest):
    valuation, lead, logarithmic = factor
    if valuation < lowest:
        raise ValueError('The pole exceeds its analytic bound.')
    return (lead if valuation == lowest else F(0),
            lead*logarithmic if valuation == lowest else lead if valuation == lowest+1 else F(0))


def scale(coeffs, regular):
    return coeffs[0]*regular[0], coeffs[1]*regular[0]+coeffs[0]*regular[1]


@cache
def open_poles(i, j, k, n, h):
    """Coefficients of ep^-3 and ep^-2 of the positive-variable moment."""
    i, j, k = sorted((i, j, k), reverse=True)
    if i < 0:
        raise ValueError('A free-line exponent must be nonnegative.')
    e = (1, -1)
    triple, double = F(0), F(0)
    for m in range(i+1):
        # The apparent angular poles cancel between the two beta normalizations.
        angular = multiply([(pochhammer(F(1, 2), -1, h), 1),
                            (pochhammer(1+h+i-m, -1, m), 1),
                            (pochhammer(*e, h+m), -1)])
        assert angular[0] >= 0
        radial = beta((2, -2), e, 2*h+j+k+n, h+i-m)
        common = F(-1, 4)**h * comb(i, m)
        x_terms = [beta(e, e, h+j+y+m, h+i-m) for y in range(n+1)]
        z_terms = [beta(e, e, h+k+z+m, h+i-m) for z in range(n+1)]
        for y, first in enumerate(x_terms):
            for z, second in enumerate(z_terms[:n-y+1]):
                moment = multiply([(v, 1) for v in (angular, radial, first, second)])
                a, b = coefficients(moment, -3)
                multinomial = (-1)**(y+z)*comb(n, y)*comb(n-y, z)
                triple += common*multinomial*a
                double += common*multinomial*b
    return triple, double


@cache
def lower_poles(terms, h):
    triple, double = F(0), F(0)
    for (i, j, k, n), coefficient in terms:
        a, b = open_poles(i-1, j-1, k-1, n, h)
        sign = (-1)**((j+k) % 2)
        triple += sign*coefficient*a
        double += sign*coefficient*b
    return triple, double


@cache
def hermite_poles(n, h):
    """Regular maximal-cut coefficient and ep^-4, ep^-3 lower coefficients."""
    top = multiply([(pochhammer(2*h, -2, n), 1), (pochhammer(4*h, -4, n), -1)])
    base, _, _ = gram_parts()
    x = base.context().gens()[3]
    off_cut = (4*x-2)*base-x*(x-1)*base.derivative(3)
    fourth, third = F(0), F(0)
    for j in range(n):
        # The descending Hermite recurrence gives this product for q_j.
        primitive = multiply([(pochhammer(j+1+2*h, -2, n-j-1), 1),
                              (pochhammer(j+4*h, -4, n-j), -1),
                              (pochhammer(F(2*h-1, 2), -1, 1), 1)])
        pole, finite = coefficients(primitive, -1)
        if not pole and not finite:
            continue
        a, b = lower_poles(entries(x**j*off_cut), h-1)
        fourth += pole*a
        third += finite*a+pole*b
    return coefficients(top, 0), (fourth, third)


@cache
def unit_poles(h):
    if h == 0:
        return F(0), F(0)
    source = h-1 if h > 0 else h
    base, _, top_gram = gram_parts()
    # B_0=x²(x-1)²/16. Its normalized maximal-cut moment is a beta ratio.
    top = multiply([(pochhammer(2*source, -2, 2), 1)]*2 +
                   [(pochhammer(4*source, -4, 4), -1)])
    v0, v1 = (v/16 for v in coefficients(top, 0))
    remainder = [F(0), lower_poles(entries(base-top_gram), source)[0]]
    for n, c in ((4, 1), (3, -2), (2, 1)):
        for j, value in enumerate(hermite_poles(n, source)[1]):
            remainder[j] += c*value/16
    if h > 0:
        return tuple(a+b for a, b in zip(remainder, scale(unit_poles(h-1), (v0, v1))))
    assert v0
    residual = tuple(a-b for a, b in zip(unit_poles(h+1), remainder))
    return scale(residual, (1/v0, -v1/v0**2))


@cache
def polar_poles(n, h):
    top, lower = hermite_poles(n, h)
    return tuple(a+b for a, b in zip(lower, scale(unit_poles(h), top)))


def one_five_pole(powers, shift):
    target = tuple(powers[i]-1 for i in RETAINED)
    if min(target) < 0:
        return F(0), F(0)
    base, _, _ = gram_parts()
    dots = sum(target)
    h = shift-dots
    value, slope = (base.context().constant(0) for _ in range(2))
    for k, polynomial in gram_jet(target):
        c0, c1 = F(1), F(0)
        for i in range(k):
            factor = F(2*shift-1-2*i, 2*(i+1))
            c0, c1 = c0*factor, c1*factor-c0/(i+1)
        term = polynomial*base**(dots-k)
        value += fmpq(str(c0))*term
        slope += fmpq(str(c1))*term
    offsets = []
    for axis, idx in enumerate(FREE[:3]):
        derivative = base.derivative(axis)
        for order in range(1, powers[idx]):
            multiplier = fmpq(2*h-1, 2)*derivative
            slope = (base*slope.derivative(axis)+multiplier*slope-derivative*value)/order
            value = (base*value.derivative(axis)+multiplier*value)/order
            h -= 1
        offsets.append(-min(powers[idx], 1))
    fourth, third = F(0), F(0)
    for (i, j, k, n), coefficient in entries(value):
        i += offsets[0]
        j += offsets[1]
        k += offsets[2]
        if max(i, j, k) < 0:
            a, b = polar_poles(n, h)
            fourth += coefficient*a
            third += coefficient*b
        else:
            sign = (-1)**((j+k) % 2)
            third += sign*coefficient*open_poles(i, j, k, n, h)[0]
    if offsets == [-1, -1, -1]:
        for (i, j, k, n), coefficient in entries(slope):
            if not (i or j or k):
                third += coefficient*polar_poles(n, h)[0]
    norm = F(8**shift, prod((2*j+1)*(2*j+2)*(2*j+3) for j in range(shift)))
    log = sum(F(2, 2*j+k) for j in range(shift) for k in (1, 2, 3))
    return fourth*norm, (third+log*fourth)*norm


def physical_poles(powers, shift):
    """Return the ep^-3,ep^-2 coefficients of ordinary M2."""
    projections = [one_five_pole(tuple(powers[i-1] for i in g), shift) for g in cut_images()]
    fourth, third = map(sum, zip(*projections))
    return -fourth/2, -third/2-3*fourth/4


if __name__ == '__main__':
    from pathlib import Path
    import json
    rows = json.loads(Path(__file__).with_name('quartic_planar_terminal.json').read_text())
    for row in rows:
        actual = physical_poles(row['uncut'], 1)
        expected = F(row['top,planar_double,elementary_double'][2])
        assert actual == (0, expected), (row, actual)
    print('All', len(rows), 'physical M2 double poles agree, with zero triple poles.')
