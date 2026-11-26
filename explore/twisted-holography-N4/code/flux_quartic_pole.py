"""Fourth normal derivative of the flux in the polarized seven-primary BCOV contact."""
from math import factorial
from pathlib import Path
from random import Random
import sys

import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'iib-bcov'))
from bulk_eight import BULK, quartics


def check():
    a,b,d,e,h,s,t,u = sp.symbols('a b d e h s t u')
    rows = sp.Matrix([[a,0,0,0],[0,b,0,0],[s,0,d,0],
                      [0,t,0,e],[u,u,h,h]])
    # Each determinant has one flux row q. With all other rows in the
    # four-plane it contributes its omitted four-row minor times q_3.
    cofactors = [(-1)**i*rows.extract([j for j in range(5) if j!=i],range(4)).det()
                 for i in range(5)]+[0,0]
    centered = [7*v-sum(cofactors) for v in cofactors]
    polynomial = BULK[0]*sum(v**2 for v in centered)**2 + BULK[1]*sum(v**4 for v in centered)
    mixed = sp.factor(sp.diff(polynomial,s,t,u).subs({s:0,t:0,u:0}))
    assert mixed == -a**2*b**2*d**2*e**2*h**2*(a+b)*(d*e-d*h-e*h)/38880
    assert mixed.subs({a:1,b:1,d:1,e:1,h:1})==sp.Rational(1,19440)
    # The tangent divergence adds (d+h) in the z2 channel and -(e+h) in z1.
    exponents = (3,2,4,3,2)
    for factor,expected in ((d+h,sp.Rational(-1,1620)),(-e-h,0)):
        jet = sp.Poly(24*mixed*factor,a,b,d,e,h)
        assert jet.coeff_monomial(exponents)==expected
    # Compare with the integral five-dimensional code on full conserved waves.
    rng = Random(420910)
    for _ in range(12):
        values = {v:rng.randint(-3,3) for v in (a,b,d,e,h,s,t,u)}
        five = [tuple(map(int,row))+(0,) for row in rows.subs(values).tolist()]
        q = (rng.randint(-3,3),rng.randint(-3,3),0,0,rng.choice((-2,-1,1,2)))
        six = [*five,q]
        seven = [*six,tuple(-sum(col) for col in zip(*six))]
        assert sum(c*v for c,v in zip(BULK,quartics(seven)))==polynomial.subs(values)*q[4]**4
    # Fourth powers of fifteen normal linear forms span the symmetric fourth power.
    powers = [(i,j) for i in range(5) for j in range(5-i)]
    coordinates = [[factorial(4)*x**i*y**j//(factorial(i)*factorial(j)*factorial(4-i-j))
                    for i,j in powers] for x,y in powers]
    assert sp.Matrix(coordinates).rank()==15
    print('Maximal fourth normal derivative, polynomial z2 component -1/1620, fifteen quartic modes.')


if __name__=='__main__':
    check()
