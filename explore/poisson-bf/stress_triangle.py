"""Bomans–Wu central stress tensor triangle in the conserved source convention."""
from itertools import product
from math import factorial
import sympy as sp

from triangle import bracket


def integral(poly, u, v):
    return sum(c * sp.Rational(factorial(i)*factorial(j), factorial(i+j+2))
               for (i,j),c in sp.Poly(poly,u,v).terms())


def verify():
    u,v,zeta = sp.symbols('u v zeta')
    p = sp.Matrix(sp.symbols('p:2'))
    q = sp.Matrix(sp.symbols('q:2'))
    k = -p-q
    det = bracket(p,q)
    y12 = u*p-v*q
    incoming = [y12-p,y12,y12+q]
    momenta = [p,q,k]
    # Coefficients of the Schwinger polynomial are kept factored in the sources.
    def directed(slots):
        factors = [sp.Poly(edge[i]+zeta*external[i],u,v)
                   for edge,external,i in zip(incoming,momenta,slots)]
        return det*integral(factors[0]*factors[1]*factors[2],u,v)
    # Multiply both published tensors by pi**2 for an exact rational comparison.
    vector_weight = -(1-2*zeta)
    extended_weight = sp.Rational(3,2)*zeta*(1-zeta)*(1-2*zeta)
    for i,j,l in product(range(2),repeat=3):
        one = directed((i,j,l))
        both = one.subs(zeta,1-zeta)-one
        vector = -det**2*(p[i]*(l-j)+q[j]*(i-l)+k[l]*(j-i))/24
        extended = -det*p[i]*q[j]*k[l]/3
        assert sp.expand(both-vector_weight*vector-extended_weight*extended) == 0
    fermion_charge = 2*zeta-1
    a = sp.Rational(3,32)*(3*fermion_charge**3-fermion_charge)
    c = sp.Rational(1,32)*(9*fermion_charge**3-5*fermion_charge)
    assert sp.expand(vector_weight-16*(a-c)) == 0
    assert sp.expand(extended_weight-sp.Rational(4,3)*(3*c-2*a)) == 0
    free = sp.Rational(1,3)
    assert 1+3*vector_weight.subs(zeta,free) == 0
    assert 3*extended_weight.subs(zeta,free) == free
    print('two Wick orientations; all eight components at symbolic momenta; vector and chiral anomalies passed')


if __name__ == '__main__':
    verify()
