"""Test the single-flux primary variation against the seven- and nine-field vertices."""
from fractions import Fraction
from itertools import product
import json
from math import factorial
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

IIB = Path(__file__).resolve().parents[2] / 'iib-bcov'
sys.path.insert(0, str(IIB))
sys.path.insert(0, str(IIB / 'operator'))
from bulk_eight import BULK, quartics
from elementary_scalar import elementary_three_loop
from planar_trace_all_genera import elementary_plane
from three_loop_soft_jet import orbit_integers


def trace(plane, strengths, normal, shears, exe):
    # Coordinates dw1, dw2, dz2, dz1, dw3; the scalar-vector probe measures
    # the scalar trace with unit pairing against the balancing momentum.
    c, d = strengths
    s, t = shears
    sources = [(a,b,0,0,0) for a,b in plane] + [(0,0,c,0,s),(0,0,0,d,t)]
    flux = (normal[0],normal[1],0,0,normal[2])
    root = tuple(-sum(xs) for xs in zip(*sources,flux))
    if len(sources) == 5:
        return (sum(x*y for x,y in zip(BULK, quartics([*sources,flux,root]))),)
    return elementary_three_loop(sources,flux,(0,0,0,-Fraction(1,d),0),exe=exe)


def highest_pole(exe, normal, radii):
    # Each minor is affine along these rays. Its two tangent rows need at
    # most one of the three sloping mixed inputs to supply a normal row.
    points, weights = [], []
    for bits in product(range(2), repeat=3):
        if not any(bits):
            continue
        for h in radii:
            i, j = (r for r in radii if r != h)
            weight = -Fraction(i*i+j*j, 2*h*(h*h-i*i)*(h*h-j*j))
            for sign in (-1, 1):
                s, t, u = (sign*h*b for b in bits)
                sources = [(1,0,0,0,0)]*2 + [(0,1,0,0,0)]*2
                sources += [(s,0,1,0,0),(0,t,0,1,0),(u,u,1,1,0),
                            (normal[0],normal[1],0,0,normal[2])]
                points.append([*sources,tuple(-sum(xs) for xs in zip(*sources))])
                weights.append((-1)**(3-sum(bits)) * sign * weight)
    rows = orbit_integers(exe, points)
    mixed = [sum(a*b for a,b in zip(weights,col)) for col in zip(*rows)]
    data = json.loads((IIB/'operator'/'three_loop_coefficients.json').read_text())
    return tuple(sum(a*Fraction(b) for a,b in zip(mixed,col))
                 for col in zip(*data['rational_and_zeta']))


def check():
    elementary = (-Fraction(49188779,2715939072000), Fraction(7783,523908000))
    j3 = (-Fraction(2078212919,635529742848000), Fraction(131371,49037788800))
    with TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'nine'
        run(['c++','-std=c++17','-O3',str(IIB/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        for q, shears in [((0,0,1),(0,0)), ((3,-1,-3),(-1,2)), ((-2,1,0),(3,-1))]:
            # The two-loop polynomial-input symbol also allows unequal strengths.
            plane = [(2,0),(-3,0),(0,-2)]
            value, = trace(plane,(2,-1),q,shears,exe)
            assert value == -Fraction(2**2 * 3**2 * 2**4 * 2**4 * q[2]**4, 31104)
            axis = trace([(1,0)]*3+[(0,1)]*2,(1,1),q,shears,exe)
            assert axis == tuple(q[2]**6 * x for x in elementary)
            for parameter in (0,1):
                planar, = trace(elementary_plane(2,parameter),(1,1),q,shears,exe)
                assert planar == -Fraction(q[2]**4,17280)
            points = [trace(elementary_plane(3,r),(1,1),q,shears,exe) for r in (0,1)]
            assert tuple((16*a-7*b)/9 for a,b in zip(*points)) == tuple(30*q[2]**6*x for x in j3)
        maximal = (Fraction(20401023829,105921623808000), -Fraction(6463627,40864824000))
        assert highest_pole(exe,(0,0,1),(1,2,3)) == maximal
        assert highest_pole(exe,(3,-4,-2),(1,2,4)) == tuple(64*x for x in maximal)
    assert Fraction(factorial(4),31104) == Fraction(1,1296)
    assert Fraction(factorial(4),17280) == Fraction(1,720)
    assert 30*factorial(6) == 21600
    print("Both full vertices reproduce the scalar-trace and normal-pole coefficients.")


if __name__ == "__main__":
    check()
