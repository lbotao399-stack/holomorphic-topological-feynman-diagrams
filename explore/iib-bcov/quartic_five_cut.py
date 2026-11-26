"""Five-residue projector onto the first gamma master of the cut K3,3 family."""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import comb, prod

from flint import fmpq, fmpq_mpoly_ctx

from quartic_nonplanar_cut import AUTOMORPHISMS, endpoint_gram
from quartic_planar_subcut import rising


# Original zero-based slots left after cutting p2²,p3²,p4²,p5²,p7².
FREE = (0, 5, 7, 8)
RETAINED = (1, 2, 3, 4, 6)


def entries(poly):
    return tuple((tuple(map(int, exponent)), Fraction(str(value))) for exponent, value in poly.terms())


@cache
def gram_parts():
    ctx = fmpq_mpoly_ctx.get(('u', 'v', 'w', 'x'))
    blocks = defaultdict(dict)
    for exponent, value in endpoint_gram(0).terms():
        key = tuple(int(exponent[i]) for i in RETAINED)
        blocks[key][tuple(int(exponent[i]) for i in FREE)] = value / 16
    blocks = {key: ctx.from_dict(terms) for key, terms in blocks.items()}
    base = blocks.pop((0,) * 5)
    x = ctx.gens()[3]
    return base, blocks, x**2 * (x-1)**2 / 16


@cache
def gram_jet(target):
    base, blocks, _ = gram_parts()
    ctx = base.context()
    parts = [(e, poly) for e, poly in blocks.items() if all(x <= y for x, y in zip(e, target))]

    @cache
    def coefficient(degree, order):
        if order == 0:
            return ctx.constant(int(not any(degree)))
        if not order <= sum(degree) <= 4 * order:
            return ctx.constant(0)
        return sum((p * coefficient(tuple(y-x for x, y in zip(e, degree)), order-1)
                    for e, p in parts if all(x <= y for x, y in zip(e, degree))), ctx.constant(0))

    return [(i, p) for i in range(sum(target) + 1) if (p := coefficient(target, i))]


def differentiated(powers, shift, s):
    """Move the three free-line powers down to simple poles before integration."""
    target = tuple(powers[i]-1 for i in RETAINED)
    base, _, _ = gram_parts()
    dots = sum(target)
    h = shift - dots
    result = base.context().constant(0)
    for order, polynomial in gram_jet(target):
        factor = prod((s+shift-j)/(j+1) for j in range(order))
        result += fmpq(str(factor)) * polynomial * base**(dots-order)
    offsets = []
    for axis, index in enumerate(FREE[:3]):
        power = powers[index]
        for order in range(1, power):
            result = (base * result.derivative(axis)
                      + fmpq(str(s+h)) * base.derivative(axis) * result) / order
            h -= 1
        offsets.append(-min(power, 1))
    return entries(result), h, offsets


def beta_ratio(a, b, i, j):
    return rising(a, i) * rising(b, j) / rising(a+b, i+j)


@cache
def open_moment(i, j, k, n, h, s):
    """Moment of a^i b^j c^k x^n B^(s+h), divided by the five-line unit cut at s.

    The positive variables a=p1², b=−p6², c=−p8² are permuted by the cut.
    For l=x+b+c, y=b/l, z=c/l, the two roots in a are
    (sqrt((1-l)(1-y)(1-z)) +/- sqrt(yz))². Integrating a first gives the
    finite binomial sum below. The other three variables give beta ratios.
    """
    i, j, k = sorted((i, j, k), reverse=True)
    if i < 0:
        raise ValueError('At least one free-line exponent must be nonnegative.')
    e = s + Fraction(3, 2)
    angular = Fraction(-1, 4)**h * beta_ratio(Fraction(1, 2), s+1, 0, h)
    common = beta_ratio(2*e, e, 2*h+j+k+n, h) * beta_ratio(e, e, h, h)**2
    total = 0
    for m in range(i+1):
        angular_polynomial = comb(i, m) * rising(e+h+i-m, m) / rising(e+h, m)
        radial = beta_ratio(2*e+2*h+j+k+n, e+h, 0, i-m)
        for y in range(n+1):
            for z in range(n-y+1):
                # x^n = l^n (1-y-z)^n in the normalized integration coordinates.
                multinomial = (-1)**(y+z) * comb(n, y) * comb(n-y, z)
                left = beta_ratio(e+h, e+h, j+y+m, i-m)
                right = beta_ratio(e+h, e+h, k+z+m, i-m)
                total += angular_polynomial * radial * multinomial * left * right
    return angular * common * total


def signed_moment(i, j, k, n, h, s):
    sign = -1 if (j+k) % 2 else 1
    return sign * open_moment(i, j, k, n, h, s)


@cache
def lower_integral(terms, h, s):
    return sum(c * signed_moment(i-1, j-1, k-1, n, h, s)
               for (i, j, k, n), c in terms)


@cache
def polynomial_reduction(n, h, s):
    """Reduce x^n/(uvw) to the same sector with no ISP and lower sectors."""
    v = s+h
    top = rising(2*v+1, n)/rising(4*v+2, n)
    if n == 0:
        return top, 0
    residual = [0]*(n+1)
    residual[n], residual[0] = 1, -top
    base, _, _ = gram_parts()
    x = base.context().gens()[3]
    off_cut = (4*x-2)*base - x*(x-1)*base.derivative(3)
    lower = 0
    for j in reversed(range(n)):
        primitive = residual[j+1]/(j+2+4*v)
        residual[j] += (j+1+2*v)*primitive
        if primitive:
            lower += primitive * lower_integral(entries(x**j * off_cut), h-1, s)
    assert not residual[0]
    return top, v*lower


@cache
def dimension_step(h, s):
    """Unit fully polar moment at h+1 in terms of the unit at h."""
    base, _, top_gram = gram_parts()
    lower = lower_integral(entries(base-top_gram), h, s)
    quartic, cubic, quadratic = (polynomial_reduction(n, h, s) for n in (4, 3, 2))
    return ((quartic[0]-2*cubic[0]+quadratic[0])/16,
            lower+(quartic[1]-2*cubic[1]+quadratic[1])/16)


@cache
def polar_moment(n, h, s):
    if h > 0:
        factor, part = dimension_step(h-1, s)
        upper, lower = polar_moment(0, h-1, s)
        upper, lower = factor*upper, part+factor*lower
    elif h < 0:
        factor, part = dimension_step(h, s)
        upper, lower = polar_moment(0, h+1, s)
        upper, lower = upper/factor, (lower-part)/factor
    else:
        upper, lower = 1, 0
    factor, part = polynomial_reduction(n, h, s)
    return factor*upper, factor*lower+part


def integrate(terms, h, offsets, shift, s):
    upper, lower = 0, 0
    for (i, j, k, n), coefficient in terms:
        i += offsets[0]
        j += offsets[1]
        k += offsets[2]
        if max(i, j, k) < 0:
            a, b = polar_moment(n, h, s)
            upper += coefficient*a
            lower += coefficient*b
        else:
            lower += coefficient*signed_moment(i, j, k, n, h, s)
    d = 2*s+5
    scale = Fraction(8**shift)/prod((d+2*j-1)*(d+2*j-2)*(d+2*j-3) for j in range(shift))
    return upper*scale, lower*scale


def one_cut(powers, shift, s):
    if any(powers[i] < 1 for i in RETAINED):
        return 0, 0
    terms, h, offsets = differentiated(powers, shift, s)
    return integrate(terms, h, offsets, shift, s)


@cache
def cut_images():
    images = {}
    for g in AUTOMORPHISMS:
        if g[0] == 0:
            images[tuple(sorted(g[i] for i in (1, 6, 8)))] = g[1:]
    assert set(images) == {(1,6,8), (2,5,7), (3,7,8), (4,5,6)}
    return tuple(images.values())


def projector(powers, shift, s, apply=one_cut):
    values = [apply(tuple(powers[i-1] for i in g), shift, s) for g in cut_images()]
    assert all(a == values[0][0] for a, _ in values)
    ep = -s-Fraction(1, 2)
    # Each of the four five-line units equals ep/(3ep-2) times ordinary M2.
    return values[0][0], sum(b for _, b in values)*ep/(3*ep-2)


def symbolic_cut(powers, shift, s):
    """Exact rational function using interpolation only for a bounded polynomial numerator."""
    if any(powers[i] < 1 for i in RETAINED):
        return 0, 0
    degree = sum(max(0, power-1) for power in powers)
    terms = defaultdict(int)
    # Every moved derivative multiplies by at most one linear factor in s.
    for j in range(degree+1):
        weight = prod((s-k)/(j-k) for k in range(degree+1) if k != j)
        samples, h, offsets = differentiated(powers, shift, Fraction(j))
        for exponents, coefficient in samples:
            terms[exponents] += weight*coefficient
    return integrate(((p, c) for p, c in terms.items() if c), h, offsets, shift, s)


def laurent_coefficient(expression, power, field):
    expression = field(expression)
    if not expression:
        return Fraction(0)
    numerator = {k[0]: Fraction(str(v)) for k, v in expression.numer.to_dict().items()}
    denominator = {k[0]: Fraction(str(v)) for k, v in expression.denom.to_dict().items()}
    lead = min(denominator)
    coefficients = {}
    for j in range(min(numerator)-lead, power+1):
        correction = sum(v*coefficients.get(j+lead-i, 0) for i, v in denominator.items() if i != lead)
        coefficients[j] = (numerator.get(j+lead, 0)-correction)/denominator[lead]
    return coefficients.get(power, Fraction(0))


def check_terminal():
    import json
    from pathlib import Path
    from sympy.polys.domains import QQ
    from sympy.polys.fields import field
    numbers, ep = field('ep', QQ)
    rows = json.loads(Path(__file__).with_name('quartic_planar_terminal.json').read_text())
    for row in rows:
        _, lower = projector(row['uncut'], 1, -ep-Fraction(1, 2), symbolic_cut)
        assert laurent_coefficient(lower, -2, numbers) == Fraction(row['top,planar_double,elementary_double'][2])
    print('All', len(rows), 'physical five-line double poles match the independent reduction.')


def check_form(folder, executable):
    from four_loop import form, read_vectors, tagged
    from shift import FIELD, K33, gram, k33, ordinary_coordinates, vector

    def no_lines(omitted):
        return tuple(int(j not in omitted) for j in range(1, 9))

    def expression(powers):
        return '*'.join(f'1/(p{i}.p{i})^{p}' if p > 0 else f'(p{i}.p{i})^{-p}'
                        for i, p in enumerate(powers, 1) if p)

    probes = [no_lines(missing) for missing in k33.MISSING]
    probes.extend(no_lines(missing) for missing in ((1,6,8), (1,), (5,), (1,6), (6,8), (1,8), (2,5,7)))
    probes.extend(((2,1,1,1,1,1,1,1), (1,2,1,1,1,2,1,1), (1,2,1,2,0,1,1,2), (1,1,1,1,-1,2,1,1)))
    for shift in (0, 1):
        inputs = [(i, expression(p)+('*'+gram(K33) if shift else '')) for i, p in enumerate(probes)]
        rows = read_vectors(form(folder, executable, 'five_controls_'+str(shift), tagged(inputs), edges=K33))
        if shift == 0:
            basis = [vector(rows[i], 0, K33) for i in range(5)]
        for i, p in enumerate(probes):
            coords = ordinary_coordinates(vector(rows[i], 0, K33), basis)
            for s in (Fraction(-8, 11), Fraction(2, 7)):
                d = 2*s+5
                scale = Fraction(8, (d-1)*(d-2)*(d-3)) if shift else 1
                expected = tuple(scale*Fraction(FIELD.to_sympy(coords[j]).subs(k33.EP, -s-Fraction(1, 2)))
                                 for j in (0, 2))
                assert projector(p, shift, s) == expected
    print('64 exact simultaneous comparisons of eight-line and five-line coordinates.')


if __name__ == '__main__':
    from argparse import ArgumentParser
    from pathlib import Path
    from tempfile import TemporaryDirectory
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--form', type=Path)
    parser.add_argument('--forcer', type=Path)
    args = parser.parse_args()
    if args.form and args.forcer:
        with TemporaryDirectory(prefix='bcov-five-cut-') as temp:
            folder = Path(temp)
            for name in ('forcer', 'forcer.h'):
                (folder/name).symlink_to((args.forcer/name).resolve(strict=True))
            check_form(folder, args.form.resolve(strict=True))
    elif args.form or args.forcer:
        parser.error('--form and --forcer are used together')
    else:
        check_terminal()
