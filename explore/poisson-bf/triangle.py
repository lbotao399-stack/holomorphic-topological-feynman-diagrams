"""Directed Poisson triangle in the paper's simplex coordinates."""
from itertools import permutations
from math import factorial
import sympy as sp

r1,r2=sp.symbols('r1 r2')
t12,t13,t23=1-r1-r2,r1,r2
p=sp.Matrix(sp.symbols('p0:2'))
q=sp.Matrix(sp.symbols('q0:2'))


def bracket(p,q):
    return p[0]*q[1]-p[1]*q[0]


def labeled(p,q):
    y12=t13*p-t23*q
    y13=p-y12
    y23=q+y12
    assert all(sp.expand(x)==0 for x in t12*y12+t23*y23-t13*y13)
    factors=[bracket(p,-y13),bracket(q,y12),bracket(-p-q,y23)]
    assert all(sp.expand(actual+t*bracket(p,q))==0 for actual,t in zip(factors,(t23,t13,t12)))
    return -bracket(p,q)**4*sp.Rational(1,120)


def verify():
    momenta=(p,q,-p-q)
    total=0
    for s in permutations(range(3)):
        sign=(-1)**sum(a>b for j,a in enumerate(s) for b in s[j+1:])
        total+=sign*labeled(momenta[s[0]],momenta[s[1]])
    assert sp.expand(total)==0
    for i in range(5):
        for j in range(5):
            value=sp.integrate(r1**(1+i)*r2**(1+j)*(1-r1-r2),(r2,0,1-r1),(r1,0,1))
            assert value/(factorial(i)*factorial(j))==sp.Rational((i+1)*(j+1),factorial(i+j+5))
    print('three Poisson vertices; simplex moments and the full alternating cochain passed')


if __name__=='__main__':
    verify()
