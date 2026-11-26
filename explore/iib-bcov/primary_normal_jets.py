"""Elementary primary scalar about a two-normal germ; exact fivefold checks."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from bulk_three_hexagon import potentials
from bulk_eight import BULK
from elementary_scalar import elementary_three_loop
from fourfold_bulk import fixture
from longitudinal_theta import primary_coefficient
from marked_boundaries import quartics
from planar_bulk import response
from quartic_scalar_response import BIPARTITE, total_channels, interval
from worldsheet_primary import fourfold_primary

BASIS = [tuple(int(i == j) for i in range(5)) for j in range(5)]
U, W, X, Q, Y = BASIS
S = (1, -1, 0, 0, 0)
V = (1, 1, 0, 0, 0)


def scaled(factor, source):
    return tuple(factor * entry for entry in source)


def add(first, second):
    return tuple(a + b for a, b in zip(first, second))


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def projected(genus, normal, exe):
    """Rooted reduced graph units; conversion to heat units supplies (-1)**genus."""
    if genus == 2:
        return primary_coefficient([U, W, *normal], Q, U), F(0)
    return elementary_three_loop([U, W, *normal], Q, U, exe=exe)


def quartic(genus, exe):
    # All g y-marks have exactly two derivatives by homogeneity and double softness.
    # The x sign mean represents cosh(x)-1; only its second and fourth jets survive.
    count = genus - 1
    result = [F(0), F(0)]
    for positive in range(count + 1):
        normal = [Y] * genus + [X] * positive + [scaled(-1, X)] * (count - positive)
        weight = F(comb(count, positive), 2**count * factorial(count) * factorial(genus))
        for i, value in enumerate(projected(genus, normal, exe)):
            result[i] += weight * value
    return tuple(result)


def cubic_pair(genus, u, v, exe):
    """Two separately labeled cubic normal jets and a unit quadratic Hessian."""
    count = 2 * genus - 3
    result = [F(0), F(0)]
    for x_marks in range(count + 1):
        quadratic = [X] * x_marks + [Y] * (count - x_marks)
        for t, s in product((-1, 1), repeat=2):
            # The double zero at every mark leaves only two spare derivatives.
            # Odd parity in both labeled slots therefore selects degree (3,3,2,...).
            normal = [*quadratic, scaled(t, u), scaled(s, v)]
            weight = F(t * s * comb(count, x_marks), 4 * factorial(count))
            for i, value in enumerate(projected(genus, normal, exe)):
                result[i] += weight * value
    return tuple(result)


def nested_fourfold():
    # In the order (A vector,C vector), polarization along B.
    # Four scalar waves remain after removing three normal Gaussian insertions.
    waves = [p[1:] for p in fixture(3)[3:-1]]
    vectors = (waves[0], waves[4])
    scalars = [waves[i] for i in (1, 2, 3, 5)]
    polarization = (0, 1, 0, 0)
    total = sum(dot(p, polarization) for p in scalars)
    result = [F(0), F(0)]
    for sign, first, second in ((-1, *vectors), (1, *reversed(vectors))):
        def accumulate(inputs, coefficient):
            root = tuple(-sum(col) for col in zip(*inputs))
            for i, value in enumerate(fourfold_primary(3, [*inputs, root])):
                result[i] += sign * total * coefficient * value

        accumulate([*scalars, first], total)
        for k, p in enumerate(scalars):
            if dot(p, polarization):
                fused = scalars.copy()
                fused[k] = add(p, first)
                accumulate([*fused, second], dot(p, polarization))
    return tuple(result)


def quartic_hexagon():
    # A two-slot Fourier sign mean of cosh(x)-1 selects its quartic response.
    # Project the original ten-point insertion; keep the root and both Taylor factorials.
    result = F(0)
    for positive in range(3):
        momenta = [U, W, Q, S, *([Y] * 3), *([X] * positive), *([scaled(-1, X)] * (2-positive))]
        root = tuple(-sum(col) for col in zip(*momenta))
        phi = potentials([*momenta, root])
        result += (F(comb(2, positive), 4 * factorial(3) * factorial(2))
                   * -dot(root, V) * (dot(phi[2], V) - dot(phi[3], V)))
    return result


def four_tangent_ward(genus):
    # Keep the unabsorbed axial covectors in the fourfold. Compute the full
    # vertex after each tangent fusion, before averaging the two normal jets.
    tangent = lambda p: (p[0], p[1], p[3], 0)
    count = genus - 1
    sums = [[F(0), F(0)] for _ in range(2)]
    for signs in product((-1, 1), repeat=count):
        axial = [(0, 0, 0, sign) for sign in signs]
        for sign, first, second in ((-1, Q, S), (1, S, Q)):
            def insert(slot, fields, contraction):
                waves = [*axial, *map(tangent,fields)]
                root = tuple(-sum(col) for col in zip(*waves))
                weight = F(sign * contraction, 2**count * factorial(count))
                for i, val in enumerate(fourfold_primary(genus,[*waves,root])):
                    sums[slot][i] += weight * val
            insert(0, (U,W,first), 4)
            insert(1, (add(U,first),W,second), 2)
            insert(1, (U,add(W,first),second), 2)
    return sums


def check_normal_coordinates():
    import sympy as sp
    u, y, x = sp.symbols('u y x')
    inverse = 3-sp.sqrt(9-6*x)
    assert sp.series(inverse**2/2,x,0,5).removeO() == x*x/2+x**3/6+5*x**4/72
    rotation = sp.Matrix([u*sp.cos(y/3)+y*sp.sin(y/3),
                          -u*sp.sin(y/3)+y*sp.cos(y/3)])
    assert sp.trigsimp(rotation.dot(rotation)-u*u-y*y) == 0
    assert sp.trigsimp(rotation.jacobian((u,y)).det()) == 1-u/3
    shear = sp.Matrix([x,y+x*x])
    assert shear.jacobian((x,y)).det() == 1


def polynomial_shear(genus, period, exe):
    # The cubic in B=(x²+(y+x²)²)/2 is x²y, and the fourth x derivative is 12.
    # 6x²y=(x+y)^3-(x-y)^3-2y^3; polarize the incoming cubics separately.
    terms = ((1,add(X,Y)), (-1,add(X,scaled(-1,Y))), (-2,Y))
    value = [12*x for x in period]
    for i,j in combinations_with_replacement(range(3),2):
        a,u = terms[i]
        b,v = terms[j]
        for k,entry in enumerate(cubic_pair(genus,u,v,exe)):
            value[k] += F(a*b,2 if i==j else 1)*entry
    return tuple(value)


def shear_tree_channels():
    """Individual classical fusions need literal jet derivatives before summation."""
    channels = {}
    powers = (2,3,2,2)  # shear generator, cubic variation and two Hessian slots
    cubics = ((1,add(X,Y)), (-1,add(X,scaled(-1,Y))), (-2,Y))
    for n in range(3):
        for coeff,v in cubics:
            axes = [X,v]+[X]*n+[Y]*(2-n)
            for values in product(*(range(d+1) for d in powers)):
                p = [U,W,Q]+[scaled(m,z) for m,z in zip(values,axes)]
                root = tuple(-sum(col) for col in zip(*p))
                p.append(root)
                # Homogeneous total normal degree nine makes this finite difference exact.
                # The 2! for the generator x² cancels the repeated Hessian Taylor divisor.
                weight = coeff*comb(2,n)*prod(F((-1)**(d-m)*comb(d,m),factorial(d))
                                             for d,m in zip(powers,values))
                for i,j in combinations(range(8),2):
                    def insertion(slot, vector):
                        if i == slot:
                            return dot(p[j],vector)
                        return dot(p[i],vector) if j == slot else 0
                    pairing = (-dot(root,Y)*(insertion(2,U)-insertion(7,U))
                               +dot(root,U)*(insertion(3,Y)-insertion(7,Y)))
                    if pairing:
                        fused = [v for h,v in enumerate(p) if h not in (i,j)]+[add(p[i],p[j])]
                        amplitude = sum(x*y for x,y in zip(BULK,quartics(fused)))
                        channels[i,j] = channels.get((i,j),F(0))+weight*pairing*amplitude
    nonzero = {pair:value for pair,value in channels.items() if value}
    assert nonzero == {(3,4):F(1,2592), (3,5):-F(1,19440), (3,6):-F(1,19440)}
    assert sum(nonzero.values()) == F(11,38880)


def shear_hexagon():
    result = F(0)
    for t,s in product((-1,1),repeat=2):
        # Exactly one x and three y Hessian slots. Odd y and even x select
        # the cubic x²y; the other normal insertions are doubly soft.
        p = [U,W,Q,X,add(scaled(t,X),scaled(s,Y)),X,Y,Y,Y]
        root = tuple(-sum(col) for col in zip(*p))
        phi = potentials([*p,root])
        first = dot(phi[2],U)-dot(phi[-1],U)
        second = dot(phi[3],Y)-dot(phi[-1],Y)
        assert first == second == 0
        result += F(s,6)*(-dot(root,Y)*first+dot(root,U)*second)
    return result


def general_shear():
    # Extract the cubic mark of the original bare graph matching polynomial.
    # Test K and J independently, including oblique quadratic Hamiltonian vectors.
    for g in range(2,7):
        count = 2*g-2
        for u,l in (((1,2),(-1,3)),((2,1),(1,0)),((1,0),(0,1))):
            v = (-l[1],l[0])
            for k,j in ((F(1),F(0)),(F(0),F(1))):
                actual = F(0)
                for n in range(count+1):
                    donors = ((1,0),)*n+((0,1),)*(count-n)
                    for sign in (-1,1):
                        plane = (l,scaled(sign,u),*donors)
                        phi = response(plane,k,j)
                        actual += F(sign*comb(count,n),2*factorial(count))*dot(phi[0][:2],v)
                expected = dot(u,v)*(k*(g*dot(l,l)*dot(u,u)-dot(l,u)**2)/factorial(g)**2
                                     +j*dot(l,u)**2)
                assert actual == expected


def main():
    five = (F(740466719, 423686495232000), -F(469153, 326918592000))
    six = (F(841479311, 211843247616000), -F(53359, 16345929600))
    with TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'nine'
        source = Path(__file__).with_name('operator') / 'nine_sextics.cpp'
        run(['c++', '-std=c++17', '-O3', str(source), '-o', str(exe)], check=True)
        for genus in (2, 3):
            period = quartic(genus, exe)
            assert period == ((F(1, 62208), F(0)) if genus == 2 else five)
            aligned = cubic_pair(genus, X, X, exe)
            assert aligned == ((F(0), F(0)) if genus == 2 else tuple(2*a-b for a,b in zip(five,six)))
            mixed = cubic_pair(genus, X, add(X, Y), exe)
            j1 = tuple(m-a for m,a in zip(mixed,aligned))
            j2 = tuple(a-b for a,b in zip(aligned,j1))
            expected_cubics = {
                2: ((F(1,77760),F(0)),(-F(1,77760),F(0))),
                3: ((F(28153087,21727512576000),-F(2971,2794176000)),
                    (-F(500673587,282457663488000),F(158827,108972864000)))}
            assert (j1,j2) == expected_cubics[genus]
            ward = four_tangent_ward(genus)
            assert ward == [[factor*x for x in period]
                            for factor in (-4,2*(4**genus-2))]
            tuned = tuple(F(5,3)*a+b/2 for a,b in zip(period,aligned))
            expected = ((F(5,186624),F(0)) if genus == 2 else
                        (F(3399295819,1271059485696000),-F(1076227,490377888000)))
            assert tuned == expected
            resistance = ((-F(7,77760),F(0)) if genus == 2 else
                          (-F(4059837269,70614415872000),F(1285913,27243216000)))
            assert j1 == tuple(-2*x-F(genus,factorial(genus)**2)*k
                               for x,k in zip(period,resistance))
            shear = polynomial_shear(genus,period,exe)
            assert tuple(a-12*b for a,b in zip(shear,tuned)) == tuple(
                F(4*genus,factorial(genus)**2)*k for k in resistance)
            assert genus != 2 or shear == (F(11,77760),F(0))
            print('Genus', genus, 'quartic', period, 'cubic trace and full contractions', j1, j2, flush=True)
            # Independently rotate both cubic axes, including a nonunit mutual pairing.
            u, v = add(X, Y), add(X, scaled(2, Y))
            expected = tuple(dot(u,u)*dot(v,v)*dot(u,v)*p + dot(u,v)**3*q
                             for p,q in zip(j1,j2))
            assert cubic_pair(genus, u, v, exe) == expected
    nested = nested_fourfold()
    graphs = ((-F(9261411221, 56491532697600), F(459733, 3405402000)),
              (-F(829523, 1975228416000), F(0)), (-F(1, 1451520), F(0)))
    projected_graphs = tuple(sum(column)/3 for column in zip(*graphs))
    assert nested == projected_graphs
    print('Six-primary nested worldsheet terms:', nested, 'graphs:', projected_graphs, flush=True)
    hexagon = quartic_hexagon()
    assert hexagon == F(1,2177280)
    assert (120*five[0]-hexagon,120*five[1]) == (F(738845099,3530720793600),-F(469153,2724321600))
    print('Quartic background hexagon projection:', hexagon, flush=True)
    four = tuple(-(F(11,576)*k+j)/6 for k,j in total_channels(BIPARTITE))
    assert tuple(6*x for x in four) == (F(100224698200645422317757433,3414152210624483328000000000),
                                       F(31025870232372083393461,395156505859315200000000),
                                       -F(93558855827167151,784040686228800000))
    assert interval(four)[0] > 0
    shear = tuple(-k/96-2*j for k,j in total_channels(BIPARTITE))
    assert interval(shear)[0] > 0
    assert tuple(a-12*b for a,b in zip(shear,four)) == tuple(k/36 for k,j in total_channels(BIPARTITE))
    print('Volume-preserving normal coordinate change, genus four:',four, interval(four), flush=True)
    check_normal_coordinates()
    shear_tree_channels()
    general_shear()
    print('Normal shear, genus-three hexagon:',shear_hexagon(),flush=True)


if __name__ == '__main__':
    main()
