"""Finite part of the three bubble products in the physical measure at dimension 11."""
from math import factorial
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]/'explore/nonSP/verify/k33_forcer'))
import replay

HALF = sp.Rational(11,2)
EULER, LOG2 = sp.symbols('Euler log2')


def gamma_expansion(z, slope):
    """Power, leading term, and first two regular log coefficients of Gamma(z+slope*t)."""
    z = sp.Rational(z)
    if z.is_integer and z<=0:
        m = -int(z)
        leading = (-1)**m / (sp.Integer(factorial(m))*slope)
        return (-1,leading,slope*(sp.harmonic(m)-EULER),
                slope**2*(sp.zeta(2)+sp.harmonic(m,2))/2)
    half = not z.is_integer
    base = sp.Rational(1,2) if half else sp.S.One
    step = int(z-base)
    linear = slope*(-EULER-2*LOG2 if half else -EULER)
    quadratic = slope**2*sp.zeta(2)*(3 if half else 1)/2
    for i in range(min(0,step),max(0,step)):
        sign = 1 if step>0 else -1
        linear += sign*slope/(base+i)
        quadratic -= sign*slope**2/(2*(base+i)**2)
    return 0,sp.gamma(z),linear,quadratic


def bubble(j):
    # GschemeConstants(j,0)=(1+j)G(1+j*ep,1), ep=-7/2+t, D/2=11/2-t.
    a = 1-sp.Rational(7*j,2)
    factors = ((a+1-HALF,j+1,1),(HALF-a,-1-j,1),(HALF-1,-1,1),
               (a,j,-1),(2*HALF-a-1,-2-j,-1))
    order,leading,linear,quadratic = 0,sp.Integer(j+1),0,0
    for argument,slope,power in factors:
        v,f,g,h = gamma_expansion(argument,slope)
        order += power*v
        leading *= f**power
        linear += power*g
        quadratic += power*h
    return order,sp.simplify(leading),sp.expand(linear),sp.expand(quadratic)


def beta_product_limit(coefficients):
    """Combine rational coefficients in t=ep+7/2, including the projective 1/Gamma(11/2)."""
    values = {i:0 for i in (-2,-1,0)}
    for exponents,coefficient in zip(((2,1,0),(2,0,1),(1,1,1)),coefficients,strict=True):
        if not coefficient:
            continue
        order,leading,linear,quadratic = 0,1,0,0
        for j,power in enumerate(exponents):
            v,f,g,h = bubble(j)
            order += power*v
            leading *= f**power
            linear += power*g
            quadratic += power*h
        terms = (leading,leading*linear,leading*(quadratic+linear**2/2))
        series = replay.laurent(coefficient,-order)
        assert not series or min(series)+order>=-2
        for n,value in series.items():
            for k,term in enumerate(terms):
                if n+order+k<=0:
                    values[n+order+k] += sp.Rational(value)*term
    assert sp.expand(values[-2])==sp.expand(values[-1])==0
    result = sp.expand(values[0]/sp.gamma(HALF))
    assert not result.free_symbols
    return result
