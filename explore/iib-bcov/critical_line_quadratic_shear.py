"""Elementary bulk crossing residues for B=((z+ax²+bxy+dy²)²+c x²y²)/2.

For c != 0 set X=b²/c and Y=(b²-4ad)/c. The coefficients below multiply
c^g in the project's reduced genus-g scalar. The transverse Gaussian on
the branch (x,0,-ax²) is dx/(x sqrt(c)), so its residue is c^(g-1/2)
times the same polynomial. The volume-preserving shear sends this curve
and the other critical branch to the two straight axes of the unsheared germ.

The genus-two zero locus is X=-(17Y²-22Y+30)/66. On it the physical
three-loop combination has one simple real root in (-1/2,-4/9).
Choosing a=1, b²=cX, d=c(X-Y)/4 kills both elementary residues. For c>0
and a,b,d real the two-loop numerator is instead strictly positive since
17Y²-22Y+30 has discriminant -1556. These are fixed-flat-coordinate
elementary operations. Their complete relative quantum transport is open.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import sympy as sp

from critical_line_massive_shear import packets, stencils
from planar_cubic_primary import CORNERS, ROOT, labeled_vertex


def monomial_stencil(degree):
    """Fourier stencils with first term x^degree y^(4-degree)/2."""
    odd = ((-1,-F(1,2)), (1,F(1,2)))
    cube = ((-2,-F(1,2)), (-1,F(1)), (1,-F(1)), (2,F(1,2)))
    pure = (((-2,0,0),F(1,2)), ((-1,0,0),-F(2)), ((1,0,0),-F(2)), ((2,0,0),F(1,2)))
    if degree in (0,4):
        rows = pure
    elif degree in (1,3):
        rows = tuple(((x,y,0),wx*wy/2) for x,wx in cube for y,wy in odd)
    else:
        return stencils()[2]
    if degree < 2:
        rows = tuple(((y,x,z),w) for (x,y,z),w in rows)
    return rows


def projected(genus, degrees, exe):
    k = len(degrees)
    quadratic,cubic,_ = stencils()
    selections = [(quadratic,k),(cubic,2*(genus-k))]
    selections.extend((monomial_stencil(j),n) for j,n in Counter(degrees).items())
    return signed_vertex(genus,selections,exe)


def signed_vertex(genus,selections,exe):
    letters = {p for stencil,_ in selections for p,w in stencil}
    actions = [{p: tuple(s*p[i] for s,i in zip(sign,order)) for p in letters}
               for sign in CORNERS for order in permutations(range(3))]
    weights = defaultdict(F)
    for factors in product(*(packets(stencil,n) for stencil,n in selections)):
        waves = sum((wave for wave,weight in factors),())
        weight = prod(w for _,w in factors)
        key = min(tuple(sorted(perm[p] for p in waves)) for perm in actions)
        weights[key] += weight
    weights = {k:w for k,w in weights.items() if w}
    vertices = labeled_vertex(genus,weights,exe)
    return tuple(sum(w*t for w,t in zip(weights.values(),col)) for col in zip(*vertices))


def quartic_channels(genus,exe):
    channels = {}
    for k in range(genus+1):
        for powers in combinations_with_replacement(range(5),k):
            if sum(powers) != 2*k:
                continue
            reflected = tuple(4-i for i in reversed(powers))
            if reflected in channels:
                channels[powers] = channels[reflected]
            else:
                channels[powers] = projected(genus,powers,exe)
            print('crossing',genus,powers,channels[powers],flush=True)
    return channels


def invariant_columns(genus,channels):
    # Set c=1. In an area-preserving frame the shear is uxy, u²=Y.
    # The original xy becomes x²+mxy+n y² by a further diagonal scaling.
    # Its discriminant is one, and u*m=b by the polarized discriminant.
    X,Y,M = sp.symbols('X Y M')
    n = (M-1)/4
    # M=m². The extra factor m for degrees 1 and 3 is supplied in pairs.
    coeff = (n*n,2*n,Y+(3*M-1)/2,2,1)
    columns = [0,0]
    for powers,values in channels.items():
        odd = sum(j%2 for j in powers)
        assert not odd%2
        multiplicity = prod(coeff[j] for j in powers)*M**(odd//2)*Y**(genus-len(powers))
        columns = [a+multiplicity*v for a,v in zip(columns,values)]
    return [sp.Poly(sp.cancel(p.subs(M,X/Y)) if p else 0,X,Y) for p in columns]


def simultaneous_zeros(two,three):
    X,Y = two[0].gens
    first = sp.Poly(66*X+17*Y**2-22*Y+30,X,Y)
    assert two[0].as_expr()*155520 == first.as_expr() and two[1].is_zero
    values = [sp.rem(p.as_expr(),first.as_expr(),X).expand() for p in three]
    for value in values:
        print('three-loop column on the two-loop zero locus',sp.factor(value),flush=True)
    assert all(sp.degree(v,Y)==3 for v in values)
    assert sp.gcd(*values) == 1
    # The two independent period columns have no common root on this conic.
    # Their physical linear combination has exactly one real root.
    z = sp.Symbol('z')
    measured = sp.Poly(values[0]+z*values[1],Y)
    partial = sum(F(1,i**3) for i in range(1,65))
    ends = (partial+F(1,2*65**2), partial+F(1,2*64**2))
    discr = sp.discriminant(measured.diff().as_expr(),Y)
    assert sp.Poly(discr,z).LC() > 0
    for endpoint in ends:
        assert discr.subs(z,endpoint) < 0 < measured.LC().subs(z,endpoint)
        assert measured.eval(-F(1,2)).subs(z,endpoint) < 0 < measured.eval(-F(4,9)).subs(z,endpoint)
    print('unique real root between -1/2 and -4/9; three distinct complex roots',flush=True)
    x,y,z,a,b,d,c = sp.symbols('x y z a b d c')
    W = ((z+a*x**2+b*x*y+d*y**2)**2+c*x**2*y**2)/2
    transverse = sp.hessian(W,(y,z)).det().subs({y:0,z:-a*x**2})
    assert sp.expand(transverse) == c*x**2
    return values


def check_stencils():
    for degree in range(5):
        rows=monomial_stencil(degree)
        for i in range(5):
            for j in range(5-i):
                if i+j<2:
                    continue
                jet=sum(w*x**i*y**j for (x,y,z),w in rows)
                target=F(factorial(i)*factorial(j),2) if (i,j)==(degree,4-degree) else 0
                assert jet==target, (degree,i,j,jet,target)


def holdout(genus,channels,exe):
    # Use the *original* x,y coordinates of the shear x²-y² at c=1.
    # For k=g-1 there are precisely two z(x²-y²) donors. The determinant
    # discriminant here is Y=4 and b=0; neither quadratic is proportional.
    quadratic = stencils()[0]
    cubic = tuple(((i,0,j),F(j,2)) for i,j in product((-1,1),repeat=2))
    cubic += tuple(((0,i,j),-F(j,2)) for i,j in product((-1,1),repeat=2))
    quartic = defaultdict(F)
    for degree,scale in ((0,1),(2,-1),(4,1)):
        for p,w in monomial_stencil(degree):
            quartic[p] += scale*w
    rows = ((quadratic,genus-1),(cubic,2),(tuple((p,w) for p,w in quartic.items() if w),genus-1))
    direct = signed_vertex(genus,rows,exe)
    quartic_coeffs = (F(1,16),0,F(7,2),0,1)
    expected = tuple(sum(4*v[j]*prod(quartic_coeffs[i] for i in powers)
                         for powers,v in channels.items() if len(powers)==genus-1) for j in range(2))
    assert direct == expected
    print('original-coordinate mixed shear holdout',genus,direct,flush=True)


if __name__ == '__main__':
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        check_stencils()
        polynomials = []
        for genus in (2,3):
            channels = quartic_channels(genus,exe)
            holdout(genus,channels,exe)
            coeffs = invariant_columns(genus,channels)
            print('invariant columns',genus,[sp.factor(p.as_expr()) for p in coeffs],flush=True)
            polynomials.append(coeffs)
        simultaneous_zeros(*polynomials)
