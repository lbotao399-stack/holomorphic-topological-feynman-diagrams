"""Integrate the doubled-edge core by cutting a lead of its bubble and reducing the kite."""
from collections import defaultdict
from functools import cache
from math import factorial, prod
from pathlib import Path
import json

import sympy as sp
from sympy.polys.fields import field

K, D = field('D', sp.QQ)
HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'reconstruct_cyclic' / 'core03-10.json'
EDGES = [[0,1],[0,4],[0,5],[1,2],[1,5],[2,3],[2,5],[3,4],[3,4]]


@cache
def shifted_bubble(a, n, shifted=False):
    """G(a,b) divided by its base at (1,b0), b=b0+n with n nonnegative."""
    b0 = 2-D/2 if shifted else K.one
    b = b0+n
    if a > 1:
        i = a-1
        return shifted_bubble(i,n,shifted)*(i+b-D/2)*(D-i-b-1)/(i*(D/2-i-1))
    if n:
        j = b-1
        return shifted_bubble(1,n-1,shifted)*(1+j-D/2)*(D-j-2)/(j*(D/2-j-1))
    return K.one


@cache
def kite(a,b,c,d,e):
    """Five lines k, Q-k, l, Q-l, k-l. Return the two bubble-product coefficients."""
    if e == 0:
        if min(a,b,c,d)==0:
            return K.zero,K.zero  # scaleless one-line factor
        return shifted_bubble(a,b-1)*shifted_bubble(c,d-1),K.zero
    if c == 0:
        return K.zero,shifted_bubble(d,e-1)*shifted_bubble(a,b+d+e-2,True)
    if d == 0:
        return K.zero,shifted_bubble(c,e-1)*shifted_bubble(b,a+c+e-2,True)
    lower_e1 = kite(a+1,b,c,d,e-1)
    lower_c = kite(a+1,b,c-1,d,e)
    lower_e2 = kite(a,b+1,c,d,e-1)
    lower_d = kite(a,b+1,c,d-1,e)
    return tuple((a*(u-v)+b*(w-x))/(D-2*e-a-b)
                 for u,v,w,x in zip(lower_e1,lower_c,lower_e2,lower_d))


def coefficients(record, swap):
    by_kite = defaultdict(lambda: K.zero)
    for exponents, coefficient in record['terms']:
        a,b,c,d,e = [exponents[i]+1 for i in (0,2,3,6,4)]
        key = (c,d,a,b,e) if swap else (a,b,c,d,e)
        weight = sp.QQ(coefficient*prod(factorial(i) for i in exponents),record['denominator'])
        by_kite[key] += weight*shifted_bubble(exponents[7]+1,exponents[8])
    first = second = K.zero
    for key, weight in by_kite.items():
        u,v = kite(*key)
        first += weight*u
        second += weight*v
    return first,second


def finite_value(first,second):
    # X=G(1,1) is regular. Y=G(1,2-D/2) has a simple pole at D=11.
    ds = sp.Symbol('D')
    half = sp.Rational(11,2)
    x = sp.gamma(2-half)*sp.gamma(half-1)**2/sp.gamma(9)
    # Residue of Gamma(3-D) at D=11 is -1/8!.
    residue_y = -sp.gamma(half-1)*sp.gamma(9)/(factorial(8)*sp.gamma(2-half)*sp.gamma(3*half-3))
    A = first.as_expr().subs(ds,11)
    B = ((second/(D-11)).as_expr()).subs(ds,11)
    assert A.is_Rational and B.is_Rational
    return sp.expand(sp.simplify((A*x**3+B*x**2*residue_y)/sp.gamma(half)))


def doubled_coefficients():
    result = {}
    for record in json.loads(SOURCE.read_text()):
        assert record['edges']==EDGES and record['orbit_representatives'] and not any(record['h'])
        first,second = coefficients(record,False)
        assert (first,second)==coefficients(record,True)
        result[record['N_power']] = finite_value(first,second)
    return result


if __name__=='__main__':
    print(doubled_coefficients())
