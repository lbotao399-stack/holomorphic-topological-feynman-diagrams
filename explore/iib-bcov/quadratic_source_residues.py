"""Local BCOV source coefficients on a quadratic and normal D3 linking moments."""
from fractions import Fraction
from math import comb, factorial, prod

from flint import fmpq_mpoly_ctx

from elliptic_descendants import vertex as elliptic_vertex
from protected_source_tower import ward_sum
from two_loop_single_descendant import single_vertex


def normal_moment(genus, order):
    """Coefficient on the unit monomial of degree 2*order with unit determinant."""
    count = genus - order
    integral = sum(Fraction(comb(count, k) * 3**k, 2*order+2*k+1)
                   for k in range(count+1))
    return integral / (4**order * 24**count * factorial(count))


def triples(total):
    for first in range(total+1):
        for second in range(total-first+1):
            yield first, second, total-first-second


def laurent_checks():
    ctx = fmpq_mpoly_ctx.get(('w1', 'w2', 'w3'))
    w = ctx.gens()
    cases = [
        ((1,0,0), (0,1,0), ((0,0,0),(0,0,0),(0,0,1)), (0,0)),
        ((1,2,-1), (0,1,2), ((1,-2,3),(-2,0,1),(3,1,-1)), (2,-3)),
        ((1,-2,1), (2,-4,2), ((1,0,0),(0,1,0),(0,0,1)), (-1,2)),
    ]
    checked = 0
    for left, right, hessian, z in cases:
        v = tuple(left[(i+1)%3]*right[(i+2)%3] - left[(i+2)%3]*right[(i+1)%3]
                  for i in range(3))
        delta = sum(v[i]*hessian[i][j]*v[j] for i in range(3) for j in range(3))
        quadratic = sum((z[0]*a+z[1]*b)*x for a,b,x in zip(left,right,w))
        quadratic += sum(hessian[i][j]*w[i]*w[j] for i in range(3) for j in range(3))/2
        powers = [{e: Fraction(str(c)) for e,c in (quadratic**m).to_dict().items()}
                  for m in range(7)]
        for genus in range(1,5):
            for source in range(2*genus-2, 3*genus+2):
                # Explicit Laurent differentiation of 1/(w1*w2*w3), followed
                # by multiplication by each ordinary quadratic power.
                residues = {}
                for j in range(genus+1):
                    m = source+1-2*genus-j
                    if m < 0:
                        continue
                    k = genus-j
                    factor = Fraction(delta**j, 24**j*factorial(j)*4**k*(2*k+1)*factorial(m))
                    for pole in triples(2*k):
                        derivative = prod(a**i for a,i in zip(v,pole))
                        for exponent, coefficient in powers[m].items():
                            test = tuple(i-l for i,l in zip(pole,exponent))
                            if min(test) >= 0:
                                residues[test] = residues.get(test,0) + factor*derivative*coefficient
                d = 3*genus-source-1
                for total in range(2*genus+1):
                    for test in triples(total):
                        expected = Fraction()
                        if 0 <= d <= genus and total == 2*d:
                            expected = normal_moment(genus,d)*delta**(genus-d)*prod(
                                a**i for a,i in zip(v,test))
                        assert residues.get(test,0) == expected, (genus,source,test)
                        checked += 1
    print('Laurent monomial tests on three brane quadratics:', checked, flush=True)


def first_jet_checks():
    unit = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    extra = ((1,1,-1,1,0),(-1,2,1,0,1))
    soft = (1,-1,2,0,1)
    for genus in (1,2):
        for hard in ((*unit,extra[0]), (*unit,*extra)):
            def value(waves):
                root = tuple(-sum(row) for row in zip(*waves))
                inputs = (root,*waves)
                if genus == 2:
                    return single_vertex(inputs, full_source=True)
                levels = (len(inputs)-5,) + (0,)*len(waves)
                return elliptic_vertex(inputs, levels)
            measured = {t: value((*hard,tuple(t*v for v in soft))) for t in (-2,-1,0,1,2)}
            assert measured[0] == value(hard)
            assert 8*(measured[1]-measured[-1]) == measured[2]-measured[-2]
            assert measured[1]+measured[-1] != 2*measured[0]
    print('Four fivefold soft-primary families have zero linear jet and nonzero curvature.', flush=True)


def ward_checks():
    # Reconstruct the two-loop derivative splittings from the five-dimensional
    # Ward diagrams, then evaluate the full-rank quadratic and the Laurent flux.
    a = ward_sum(2,(1,0))[0]
    plus = ward_sum(2,(1,1))[0]
    minus = ward_sum(2,(1,-1))[0]
    odd = (plus-minus)/4
    mixed = (plus+minus)/2-2*a
    assert (a,odd,mixed)==(Fraction(1,1920),0,Fraction(1,576))
    # At the normal origin the extra axial quadratic may be retained once,
    # twice (identical-field 2!), or left as its value w3^2/2 by string.
    for d,expected in enumerate((Fraction(1,240),Fraction(7,720),Fraction(1,80))):
        # The coefficient of w3^(2d) is a finite Laurent residue.
        source = 5-d
        direct = Fraction()
        for j,base in enumerate((a,mixed,mixed/2)):
            m = source-3-j
            if m >= 0:
                direct += base * factorial(4-2*j) / (2**m*factorial(m))
        assert direct == expected == normal_moment(2,d)
    print('Two-loop constant, quadratic and quartic unit linking tests:',
          *(normal_moment(2,d) for d in range(3)), flush=True)
    elliptic = ward_sum(1, (1,))[0]
    assert elliptic == Fraction(1, 24)
    for genus, local, expected in [
            (1, (elliptic, elliptic), (Fraction(1, 6), Fraction(1, 4))),
            (2, (a, mixed, mixed/2), (Fraction(1, 72), Fraction(1, 24), Fraction(1, 16)))]:
        rows = cubic_laurent_rows(genus, 2, local)
        assert tuple(rows) == expected
        print('Fivefold graph weights on the axial double pole:', genus, *rows, flush=True)


def root_checks():
    stable = 0
    for genus in range(1,15):
        source = 3*genus-1
        # Directly combine the individual flat Laurent moments.
        measured = sum(Fraction((-1)**d, 2**d*factorial(d))*normal_moment(genus,d)
                       for d in range(genus+1))
        assert measured == Fraction(1, 24**genus*factorial(genus))
        leaf = sum(Fraction((-1)**j, factorial(j)*factorial(source+1-j))
                   for j in range(source+1))
        assert leaf == Fraction((-1)**source, factorial(source+1))
        assert source+1 > genus
        for loops in range(1,genus+1):
            for marked_genus in range(genus-loops+1):
                for exceptional_at_root in (0,1):
                    for degree in range(1,2*loops+2-exceptional_at_root):
                        hard = degree + exceptional_at_root
                        if marked_genus == 0 and hard < 2:
                            continue
                        donors = 2*marked_genus+2+source-hard
                        assert 2*donors > 10*marked_genus
                        stable += 1
    print('Marked-vertex exclusions through genus fourteen:', stable, flush=True)


def cubic_laurent_rows(genus, axis, local_weights):
    """Pair a scalar with a double pole on axis i against w_i w3^(2j)."""
    pole_order = 1 + int(axis == 2)
    rows = []
    for shift in range(genus + 1):
        moment = Fraction()
        for v, coefficient in enumerate(local_weights):
            derivatives = 2*(genus-v)
            powers = genus-v-shift
            if powers < 0:
                continue
            differentiated = prod(range(pole_order, pole_order+derivatives))
            moment += coefficient * differentiated / (2**powers*factorial(powers))
        rows.append(moment)
    return rows


def cubic_laurent_checks():
    for genus in range(1, 15):
        local = [Fraction(1, 24**v*factorial(v)*4**(genus-v)*factorial(2*(genus-v)+1))
                 for v in range(genus+1)]
        for axis in range(3):
            rows = cubic_laurent_rows(genus, axis, local)
            if genus == 1 and axis == 2:
                assert rows == [Fraction(1, 6), Fraction(1, 4)]
            if genus == 2 and axis == 2:
                assert rows == [Fraction(1, 72), Fraction(1, 24), Fraction(1, 16)]
            measured = sum(value*Fraction((-1)**j, 2**j*factorial(j))
                           for j, value in enumerate(rows))
            assert measured == Fraction(1, 24**genus*factorial(genus))
            for tested_axis in range(3):
                # Differentiation and multiplication restore every original
                # pole. Multiplication by the test w_J then selects I=J.
                residual_powers = tuple(-1-int(i==axis)+int(i==tested_axis) for i in range(3))
                assert (residual_powers == (-1,-1,-1)) == (axis == tested_axis)
    print('All cubic Laurent scalar pairings through genus fourteen agree.', flush=True)


def regular_ward_checks():
    """Move derivatives onto the full regular source before filtering its axial degree."""
    ctx = fmpq_mpoly_ctx.get(('z1', 'z2', 'w1', 'w2', 'w3', 'c'))
    z1, z2, w1, w2, w3, c = ctx.gens()
    mixed = z1*w1 + z2*w2
    axial = c*w3**2/2
    for genus in range(1, 5):
        for closed_genus in range(genus+1):
            residual = 3*(genus-closed_genus)
            for index, normal in enumerate((w1,w2,w3)):
                projected = ctx.constant(0)
                prediction = ctx.constant(0)
                extra = int(index < 2)
                for k in range(closed_genus+1):
                    det_power = closed_genus-k
                    denom = 24**det_power*factorial(det_power)*4**k*factorial(2*k+1)
                    for j in range(min(genus,residual+k)+1):
                        exponent = residual+k-j
                        term = normal*(-axial)**j*(mixed+axial)**exponent
                        for _ in range(k):
                            term = 2*c*(term.derivative(0).derivative(2)
                                        + term.derivative(1).derivative(3)) + term.derivative(4).derivative(4)
                        projected += c**det_power*term/(denom*factorial(j)*factorial(exponent))
                    numerator = factorial(residual+k+1+extra)
                    prediction += (c**closed_genus*normal*mixed**residual*numerator
                                   / (factorial(residual)*factorial(residual+1+extra)*(denom//2**k)))
                assert all(power[-1]>genus for power in (projected-prediction).to_dict())
    print('Regular separating Ward sources at all closed genera through genus four agree.', flush=True)


def annular_potential_checks():
    # Killing the two tangent coordinates on the delta potential leaves
    # a single normal variable. Differentiate before taking its divergence.
    ctx = fmpq_mpoly_ctx.get(('w3', 'c'))
    w3, c = ctx.gens()
    for measured in range(1,10):
        level = 3*measured-1
        for loop in range(0,level//2+2):
            coefficient = ctx.constant(0)
            for k in range(loop+1):
                j = loop-k
                m = level+1-2*loop-j
                if m < 0:
                    continue
                term = (c*w3**2/2)**m
                for _ in range(2*k):
                    term = term.derivative(0)
                coefficient += (term*c**j
                                / (factorial(m)*24**j*factorial(j)*4**k*factorial(2*k+1)))
            if loop == measured:
                expected = sum(Fraction(factorial(2*k),
                                        2**k*factorial(k)*24**(loop-k)*factorial(loop-k)
                                        *4**k*factorial(2*k+1)) for k in range(loop+1))
                assert coefficient == c**measured*expected.numerator/expected.denominator
            if loop > measured:
                assert coefficient == 0
            assert all(power[1]>=measured+2 for power in coefficient.derivative(0).to_dict())
    print('Annular primary germs through measured genus nine have the predicted axial order.', flush=True)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--ward', action='store_true', help='reconstruct all two-loop local linking moments from Ward graphs')
    args = parser.parse_args()
    laurent_checks()
    first_jet_checks()
    root_checks()
    cubic_laurent_checks()
    annular_potential_checks()
    if args.ward:
        regular_ward_checks()
        ward_checks()
