"""Integer-polynomial two-loop kite reduction in z = D/2."""
from functools import cache
from flint import fmpz_poly

Z = fmpz_poly([0,1])
ONE = fmpz_poly([1])
ZERO = fmpz_poly()


def reduce_pair(x,y):
    g = x.gcd(y)
    if y.leading_coefficient()*g.leading_coefficient() < 0:
        g = -g
    return x//g,y//g


@cache
def bubble(a,n,shifted):
    b = (2-Z if shifted else ONE)+n
    if a>1:
        i = a-1
        x,y = bubble(i,n,shifted)
        return reduce_pair(x*(i+b-Z)*(2*Z-i-b-1), y*i*(Z-i-1))
    if n:
        j = b-1
        x,y = bubble(1,n-1,shifted)
        return reduce_pair(x*(1+j-Z)*(2*Z-j-2), y*j*(Z-j-1))
    return ONE,ONE


def product(a,b):
    x,y = a
    w,t = b
    return reduce_pair(x*w,y*t)


def native_kite(a,b,c,d,e):
    # Exchange loop momenta and reverse both propagator pairs. Among the
    # equivalent orientations minimize the sum reduced by the recurrence.
    orient = ((a,b,c,d),(b,a,d,c),(c,d,a,b),(d,c,b,a))
    a,b,c,d = min(orient,key=lambda v:(min(v[2:])!=0,sum(v[2:]),v))
    return coefficients(a,b,c,d,e)


@cache
def coefficients(a,b,c,d,e):
    if e==0:
        if min(a,b,c,d)==0:
            return ZERO,ZERO,ONE
        x,y = product(bubble(a,b-1,False),bubble(c,d-1,False))
        return x,ZERO,y
    if c==0 or d==0:
        if d==0:
            a,b,c,d = b,a,d,c
        x,y = product(bubble(d,e-1,False),bubble(a,b+d+e-2,True))
        return ZERO,x,y
    left = right = ZERO
    denominator = ONE
    for weight,child in ((a,(a+1,b,c,d,e-1)),(-a,(a+1,b,c-1,d,e)),
                         (b,(a,b+1,c,d,e-1)),(-b,(a,b+1,c,d-1,e))):
        x,y,t = native_kite(*child)
        g = denominator.gcd(t)
        new,old = t//g,denominator//g
        left = left*new+weight*x*old
        right = right*new+weight*y*old
        denominator *= new
    denominator *= 2*Z-2*e-a-b
    g = denominator.gcd(left).gcd(right)
    if denominator.leading_coefficient()*g.leading_coefficient()<0:
        g = -g
    return left//g,right//g,denominator//g
