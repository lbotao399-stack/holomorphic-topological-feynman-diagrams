"""Maximal D3-flux pole of the saturated BCOV trace on closed scalar inputs."""
from fractions import Fraction
import json
from itertools import combinations
from math import comb, factorial, prod
from pathlib import Path
import sys
from subprocess import run
from tempfile import TemporaryDirectory

IIB = Path(__file__).resolve().parents[2] / 'iib-bcov'
sys.path[:0] = [str(IIB / 'operator'), str(IIB)]
from flint import fmpz_mat, fmpz_mpoly_ctx
from four_loop_planar import restricted_symbols
from quartic_scalar_response import BIPARTITE, interval, total_channels


def matching_coefficient(n, tangent):
    """Coefficient of s_0^2 s_1 ... s_(2n-2) in the three planar symbols."""
    length = 2*n-2
    first, *rest = tangent
    assert len(rest) == length
    scale = (-2)**(n-2)*prod(range(1, length-2, 2))*prod(rest)
    total = sum(rest)
    squares = sum(x*x for x in rest)
    a = 2*(length-1)*(length*first**2-squares)
    d = (length-1)*first*total-total**2+squares
    return tuple(scale*x for x in (a, 2*a+2*(2*length-3)*d, d))


def projected_period(n):
    """K_n and J_n multipliers before differentiating the unit normal residue."""
    a, b, d = matching_coefficient(n, [2]+[1]*(2*n-2))
    assert 2*(3*n-1)*a-5*b+2*(2*n+1)*d == 0
    assert a == 6*d and 3*b == (4*n-1)*a
    cross = (-1)**n * comb(2*n,n)
    k = Fraction(cross*b, 2*n*(n-1)*factorial(n))
    j = Fraction(cross*factorial(n)*(2*n*a-b), 2*n*(n-1)*(2*n+1))
    common = Fraction(comb(2*n,n)*factorial(2*n-2), factorial(n))
    assert (k,j) == (common*Fraction(4*n-1,factorial(n)), common*factorial(n))
    return k,j


def verify_five_row_minors():
    # Two complementary mixed waves live in (dw_2,dz_2); the flux is
    # allowed arbitrary normal momentum. Test every minor, including recoil.
    plane = [(2,2),(1,1),(-3,1),(4,1),(0,1),(2,1),(-1,1)]
    for left,right,normal in ((2,-3,(4,-2,5)),(3,3,(1,4,-1)),(-1,2,(2,-2,0))):
        base = [(x,y,0,0,0) for x,y in plane]
        base.extend([(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)])
        actual = [(x,y,0,0,0) for x,y in plane]
        actual.extend([(0,0,left,1,0),(0,0,right,1,0),
                       (normal[0],0,normal[1],0,normal[2])])
        for rows in (base,actual):
            rows.append(tuple(-sum(xs) for xs in zip(*rows)))
        for subset in combinations(range(11),5):
            original, transformed = (fmpz_mat([rows[i] for i in subset]).det()
                                     for rows in (base,actual))
            assert transformed == original*(left-right)*normal[2]


def check_integrated_vertices():
    from bulk_eight import BULK, quartics
    from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers

    for n in (2,3):
        points, weights = [], []
        slots = 2*n-2
        # Difference once in every identical slope, twice in the first.
        # Their total degree is exactly 2n, so this extracts one monomial.
        for first in range(3):
            for active in range(slots+1):
                plane = [(first,2)]+[(1,1)]*active+[(0,1)]*(slots-active)
                rows = [(x,y,0,0,0) for x,y in plane]
                rows.extend([(0,0,1,1,0),(0,0,0,1,0),(0,0,0,0,1)])
                rows.append(tuple(-sum(xs) for xs in zip(*rows)))
                points.append(rows)
                weights.append(Fraction((-1)**(first+slots-active)*comb(2,first)*comb(slots,active),2))
        if n == 2:
            measured = [(sum(v*w for v,w in zip(BULK,quartics(rows))),) for rows in points]
            periods = [(Fraction(-7,77760), Fraction(-1,20736))]
        else:
            with TemporaryDirectory() as tmp:
                native = Path(tmp)/'nine'
                run(['c++','-std=c++17','-O3',str(IIB/'operator'/'nine_sextics.cpp'),
                     '-o',str(native)],check=True)
                values = orbit_integers(native,points)
            data = json.loads((IIB/'operator'/'three_loop_coefficients.json').read_text())
            assert tuple(map(tuple,data['basis_words'])) == DOUBLE_SOFT
            columns = [tuple(map(Fraction,col)) for col in zip(*data['rational_and_zeta'])]
            measured = [tuple(sum(v*w for v,w in zip(row,col)) for col in columns) for row in values]
            periods = [(Fraction(-4059837269,70614415872000), Fraction(-2078212919,635529742848000)),
                       (Fraction(1285913,27243216000), Fraction(131371,49037788800))]
        k,j = projected_period(n)
        for column,(first_period,second_period) in zip(zip(*measured),periods,strict=True):
            original = sum(x*y for x,y in zip(weights,column)) * (-1)**n*comb(2*n,n)
            assert original == k*first_period+j*second_period


def check():
    for n in range(2,6):
        ctx = fmpz_mpoly_ctx.get(tuple('s'+str(i) for i in range(2*n-1)))
        slopes = ctx.gens()
        for tangent in ([2]+[1]*(2*n-2), [(-1)**i*(i+1) for i in range(2*n-1)]):
            symbols = restricted_symbols(list(zip(slopes, tangent)))
            actual = tuple(int(symbol.to_dict().get((2,)+(1,)*(2*n-2),0)) for symbol in symbols)
            assert actual == matching_coefficient(n,tangent)
        left,right = slopes[:2]
        exponent = (n,n)+(0,)*(2*n-3)
        assert ((left-right)**(2*n)).to_dict()[exponent] == (-1)**n*comb(2*n,n)
    for n in range(2,31):
        projected_period(n)
    verify_five_row_minors()
    check_integrated_vertices()
    assert matching_coefficient(4, [2,1,1,1,1,1,1]) == (2160,10800,360)
    k,j = projected_period(4)
    assert (k,j) == (Fraction(2625,2),50400)
    symbol = tuple(k*a+j*b for a,b in total_channels(BIPARTITE))
    assert interval(symbol)[1] < 0
    print('Genus four, scalar trace in the 1, zeta(3), zeta(5) basis:', *symbol)
    print('Ninth-order Laurent pole per unit flux:', *(factorial(8)*x for x in symbol))
    print('Every five-row minor and both planar coefficient extractions passed.')


if __name__ == '__main__':
    check()
