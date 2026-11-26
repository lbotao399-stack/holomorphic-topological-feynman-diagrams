"""Bracket of a tangent Laurent vector and scalar states in the three-shear complex."""
from fractions import Fraction
from itertools import product
from math import factorial, prod

import projective_flux as form


CYCLIC = ((0, 1, 2), (1, 2, 0), (2, 0, 1))


def oriented_pair(i, j):
    return (1 << i) | (1 << j), 1 if i < j else -1


def differentiate(value, orders):
    for order, coordinate in zip(orders, form.w):
        for _ in range(order):
            value = form.derivative(value, coordinate)
    return value


def differential(terms, active=(0, 1, 2)):
    result = [{} for _ in range(8)]
    for mask, value in enumerate(terms):
        signed = (-1)**mask.bit_count()
        result[mask] = form.add(result[mask], form.scale(form.differential(value), signed),
                                form.scale(form.divergence(value), signed))
        for i in active:
            bit = 1 << i
            if bit & mask:
                continue
            sign = (-1)**((bit-1) & mask).bit_count()
            shear = form.scale(form.xi[0], form.w[i])
            action = form.scale(form.bracket(shear, value), sign)
            result[bit | mask] = form.add(result[bit | mask], action)
    return result


def vector_lift():
    states = [{} for _ in range(8)]
    direction = form.xi[1]
    states[0] = form.scale(form.wedge(form.U, direction), form.z1**2/2)
    for i, j, k in CYCLIC:
        states[1 << i] = form.scale(form.wedge(form.b[i], direction), form.z1)
        mask, sign = oriented_pair(i, j)
        states[mask] = form.scale(direction, sign * form.v[k] * form.u)
    return states


def check_angular(orders):
    value = differentiate(form.U, orders)
    scalar = form.scale(value, form.z2)
    vector = vector_lift()
    for coefficient in [*differential(vector), *differential([scalar]+[{}]*7)]:
        form.check_zero(coefficient)
    # The right input has even form degree and no exterior parameter.
    bracket = [form.bracket(part, scalar) for part in vector]
    for coefficient in differential(bracket):
        form.check_zero(coefficient)
    for mask in (0, 1, 2, 4, 7):
        form.check_zero(bracket[mask])
    weighted = {}
    for i, j, k in CYCLIC:
        mask, sign = oriented_pair(i, j)
        c = form.scale(bracket[mask], sign)
        form.check_zero(form.add(c, form.scale(value, -form.v[k]*form.u)))
        weighted = form.add(weighted, form.scale(c, form.w[k]))
        higher = list(orders)
        higher[k] += 1
        weight = sum(orders) + 3
        form.check_zero(form.add(form.scale(c, weight), differentiate(form.U, higher)))
        # A primitive on this oriented pair, with the derivative taken before pullback.
        divisor = weight * (higher[i]+1)
        higher[i] += 1
        primitive = [{} for _ in range(8)]
        primitive[1 << j] = form.scale(differentiate(form.U, higher), form.z1/divisor)
        primitive[mask] = form.scale(differentiate(form.b[i], higher), -form.R.one*sign/divisor)
        for index, image in enumerate(differential(primitive, (i, j))):
            target = bracket[mask] if index == mask else {}
            form.check_zero(form.add(image, form.scale(target, -1)))
    form.check_zero(form.add(weighted, form.scale(value, -1)))


def check_laurent(powers, poles):
    # The Čech pullback of d_y^orders (1/prod y_i) uses no volume Jacobian.
    orders = [(h-1)//p for p, h in zip(powers, poles)]
    prefactor = [p*(a+1)-h for a, p, h in zip(orders, powers, poles)]
    normal = (-1)**sum(orders)*prod(factorial(a) for a in orders)
    total = Fraction(0)
    for k in range(3):
        higher = [a + int(i == k) for i, a in enumerate(orders)]
        derivative = (-1)**sum(higher)*prod(factorial(a) for a in higher)
        coefficient = Fraction(-derivative, normal*(sum(orders)+3))
        assert coefficient == Fraction(orders[k]+1, sum(orders)+3)
        exponents = [b-p*(a+1) for b,p,a in zip(prefactor,powers,higher)]
        exponents[k] += powers[k]
        assert exponents == [-h for h in poles]
        total += coefficient
    # Pairing with w^(poles-1) extracts the three simple poles.
    assert total == 1


def check_linear_vector(orders):
    degree = sum(orders)
    value = differentiate(form.U, orders)
    normal_euler = form.add(*(form.scale(direction, coord) for direction, coord in
                              zip(form.xi[2:], form.w)))
    euler = form.add(normal_euler, form.scale(form.xi[0], form.z1),
                     form.scale(form.xi[1], (degree-1)*form.z2))
    vector = form.wedge(value, euler)
    scalar = [{} for _ in range(8)]
    scalar[0] = form.scale(form.U, form.z1**2/2)
    for i, j, k in CYCLIC:
        scalar[1 << i] = form.scale(form.b[i], form.z1)
        mask, sign = oriented_pair(i, j)
        scalar[mask] = {0: sign*form.v[k]*form.u}
    form.check_zero(form.divergence(vector))
    for coefficient in [*differential([vector]+[{}]*7), *differential(scalar)]:
        form.check_zero(coefficient)
    # An even Lie action (a two-form times a vector) passes exterior parameters.
    output = [form.bracket(vector, part) for part in scalar]
    for coefficient in differential(output):
        form.check_zero(coefficient)
    for mask in (0, 1, 2, 4, 7):
        form.check_zero(output[mask])
    weighted = {}
    for i, j, k in CYCLIC:
        mask, sign = oriented_pair(i, j)
        c = form.scale(output[mask], sign)
        form.check_zero(form.add(c, form.scale(value, form.v[k]*form.u)))
        weighted = form.add(weighted, form.scale(c, form.w[k]))
    form.check_zero(form.add(weighted, value))
    # Multiplication by w_s on a Laurent derivative records the ordinary normal
    # vector projection.  Distinct multi-indices are distinguished at fixed s.
    for s, order in enumerate(orders):
        if order:
            coefficient = (-1)**degree * prod(factorial(a) for a in orders)
            lower = list(orders)
            lower[s] -= 1
            assert coefficient == -order*(-1)**(degree-1)*prod(factorial(a) for a in lower)


def check():
    for orders in ((0,0,0),(0,0,1),(2,0,1),(1,1,1)):
        check_angular(orders)
    for orders in product(range(4), repeat=3):
        if 1 <= sum(orders) <= 3:
            check_linear_vector(orders)
    for powers in ((1,1,1),(1,1,2),(1,1,7),(2,3,4)):
        for poles in product(range(1,8), repeat=3):
            check_laurent(powers, poles)
    print('Three-shear bracket closed; all pair restrictions exact; 1372 normalized Laurent residues equal 1.')
    print('All 19 transverse derivatives of orders one through three have closed linear-shear inputs and the predicted bracket.')


if __name__ == '__main__':
    check()
