"""Flux-position Gaussian integration for the separated one-loop pentagons."""
from itertools import permutations
from math import factorial, prod

import sympy as sp
from sympy.integrals.rationaltools import ratint_ratpart
from puncture_pentagon import POINTS, UNIT, coefficient


def angular(points, order, derivative):
    path = (0, *order, 4)
    displacements = [points[j] - points[i] for i, j in zip(path, path[1:])]

    def absolute_row(vertex, covector):
        tail = [covector.dot(d) if k < path.index(vertex) else 0
                for k, d in enumerate(displacements)]
        return sp.Matrix([*covector, *tail]).T

    def relative(covector):
        return absolute_row(0, covector) - absolute_row(4, covector)

    x, y = points[0][2:, :], points[4][2:, :]
    anti = [relative(UNIT[:, k]) for k in (0, 1)]
    anti += [absolute_row(0, sp.Matrix([0, 0, *y.cross(x)])),
             absolute_row(4, sp.Matrix([0, 0, *y]))]
    normal = points[1][2:, :]
    # Expand i_normal d^3(wbar) directly, valid also on every coordinate hyperplane.
    test = [absolute_row(2, UNIT[:, k]) for k in (0, 1, 2)]
    flux_anti = sum((-1)**k * normal[k] * sp.Matrix.vstack(
        *anti, *[absolute_row(1, UNIT[:, 2+j]) for j in range(3) if j != k], *test
    ).det(method='domain-ge') for k in range(3))
    odd = [relative(UNIT[:, k]) for k in range(5) if k != derivative]
    odd += [absolute_row(1, UNIT[:, k]) for k in (0, 1)]
    odd += [absolute_row(2, UNIT[:, k]) for k in (2, 3, 4)]
    return flux_anti * (-1)**derivative * sp.Matrix.vstack(*odd).det(method='domain-ge')


def flux_polynomials(order, derivative):
    t, h = sp.symbols('t h')
    path = (0, *order, 4)
    k = path.index(1)
    left, right = POINTS[path[k-1]], POINTS[path[k+1]]
    mean = (1-t)*left + t*right
    points = list(POINTS)
    points[1] = mean[:2, :].col_join((1-h)*mean[2:, :])
    numerator = sp.factor(angular(points, order, derivative))
    variance = sp.expand(t*(1-t)*(left-right).dot(left-right) + h*mean[2:, :].dot(mean[2:, :]))
    fixed = prod((POINTS[i]-POINTS[j]).dot(POINTS[i]-POINTS[j])
                 for i, j in zip(path, path[1:]) if 1 not in (i, j))
    return numerator, variance, fixed


def radial_integral(numerator, variance, fixed):
    """Integrate the relative inverse width of the flux Gaussian."""
    t, h = sp.symbols('t h')
    a, b = variance.subs(h, 0), variance.diff(h)
    coefficients = sp.Poly(sp.cancel(h**2*numerator/(1-h)), h)
    total = 0
    for (k,), coefficient in coefficients.terms():
        assert 0 <= k <= 8
        radial = sum(sp.binomial(8-k, j)*(-b)**j/((k+j+1)*(a+b)**(k+j+1))
                     for j in range(9-k))
        total += coefficient * radial / a**(9-k)
    return sp.factor(t**5*(1-t)**5*total/fixed**6)


def hermite_integral(rational):
    t = sp.Symbol('t')
    p, q = sp.fraction(rational)
    primitive, logarithmic = ratint_ratpart(sp.Poly(p, t), sp.Poly(q, t), t)
    assert sp.cancel(primitive.diff(t) + logarithmic - rational) == 0
    constant = primitive.subs(t, 1) - primitive.subs(t, 0)
    p, q = sp.fraction(logarithmic)
    poly = sp.Poly(q, t)
    linear, quadratic = sp.Poly(p/poly.LC(), t), poly.monic()
    a, b = quadratic.nth(1), quadratic.nth(0)
    discriminant = a*a-4*b
    logarithm = linear.nth(1)/2
    radical_logarithm = linear.nth(0) - a*logarithm
    return constant, a, b, discriminant, logarithm, radical_logarithm


def certificate():
    groups = {}
    for order in permutations((1, 2, 3)):
        for derivative in range(2):
            anti, odd, _ = coefficient(order, derivative)
            assert angular(POINTS, order, derivative) == anti*odd
            n, variance, fixed = flux_polynomials(order, derivative)
            key = derivative, variance, fixed
            groups[key] = groups.get(key, 0) + n
    rows = {}
    for (direction, variance, fixed), n in groups.items():
        if n == 0:
            continue
        rational = radial_integral(n, variance, fixed)
        const, a, b, disc, log, simple = hermite_integral(rational)
        assert log == 0
        rows[direction, fixed] = const, simple, a, b
        assert disc == a**2-4*b

    # Rational bounds on integral_0^1 dt/(t^2+a*t+b) from the atanh series.
    lower = upper = 0
    for fixed in (12, 6, 8):
        const, simple, a, b = rows[0, fixed]
        assert rows[1, fixed] == (-const, -simple, a, b)
        c, disc = -2*b-a, a**2-4*b
        ratio = disc/c**2
        assert 0 < ratio < 1 and c > 0
        partial = sum(ratio**j/(2*j+1) for j in range(24))
        tail = ratio**24/(49*(1-ratio))
        signed_tail = -2*simple*tail/c
        middle = const - 2*simple*partial/c
        lower += middle + min(0, signed_tail)
        upper += middle + max(0, signed_tail)
        print(f"quadratic ({a}, {b}), rational term {const}, simple-pole coefficient {simple}")
    factor = sp.Rational(factorial(9), 2*factorial(5)**2)
    assert factor == sp.Rational(63, 5)
    assert 1 < 10**12*factor*lower < 10**12*factor*upper < sp.Rational(11, 10)
    print("After flux integration the two coefficients are opposite, with 1 < 10^12 I_1/pi^5 < 11/10.")


if __name__ == '__main__':
    certificate()
