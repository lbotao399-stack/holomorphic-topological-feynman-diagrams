"""Finite projected normal-coordinate Ward transport, including a rank-one Hessian."""
from itertools import product
from math import factorial

import sympy as sp


def jets(B, w):
    zero = dict.fromkeys(w,0)
    def entries(order):
        return {index:sp.diff(B,*(w[i] for i in index)).subs(zero)
                for index in product(range(2),repeat=order)}
    return entries(2),entries(3),entries(4)


def contractions(B, w):
    h,t,u=jets(B,w)
    H=sp.Matrix(2,2,lambda i,j:h[i,j])
    inv=H.inv()
    traced=[sum(inv[i,j]*t[i,j,k] for i,j in product(range(2),repeat=2)) for k in range(2)]
    trace_square=(sp.Matrix(traced).T*inv*sp.Matrix(traced))[0]
    full=sum(t[i,j,k]*t[l,m,n]*inv[i,l]*inv[j,m]*inv[k,n]
             for i,j,k,l,m,n in product(range(2),repeat=6))
    fourth=sum(u[i,j,k,l]*inv[i,j]*inv[k,l] for i,j,k,l in product(range(2),repeat=4))
    return H.det(),trace_square,full,fourth


def check_composition():
    x,y=sp.symbols('x y')
    # Noncommuting volume-preserving shears; drop all degrees above four.
    def pull(B,f):
        value=B.subs({x:f[0],y:f[1]},simultaneous=True)
        return sum(c*x**i*y**j for (i,j),c in sp.Poly(value,x,y).terms() if i+j<=4)
    B=(x*x+2*x*y+3*y*y)/2+(2*x**3-3*x*y*y+y**3)/6+x**4/8
    f=(x+y*y,y)
    h=(x,y+2*x*x)
    second=pull(B,f)
    third=pull(second,h)
    reverse=pull(pull(B,h),f)
    assert third != reverse
    k,j=sp.symbols('K J')
    for g in (2,3,4,6):
        def split(field):
            d,t,s,u=contractions(field,(x,y))
            intrinsic=d**g*(u-t-2*s/3)
            primitive=d**g*(-g*k*t/(2*factorial(g)**2)+(k/factorial(g)**2-j)*s/6)
            return intrinsic,primitive
        a,b,c,r=map(split,(B,second,third,reverse))
        assert a[0]==b[0]==c[0]==r[0]
        assert sp.expand((a[1]-b[1])+(b[1]-c[1])-(a[1]-c[1])) == 0
        assert sp.expand(c[1]-r[1]) == 0  # their commutator starts at coordinate degree three


def rank_one():
    a,b,c,d,U,r,s=sp.symbols('a b c d U r s')
    I=sp.Rational(1,62208)
    def two(B,C,D,fourth):
        return I*a*a*fourth+a*(B*D-C*C)/77760
    invariant=a*a*U+3*a*(b*d-c*c)
    transformed=(b+2*a*s,c+a*r,d,U+6*(c*r-d*s)+3*a*r*r)
    bp,cp,dp,up=transformed
    assert sp.expand(a*a*up+3*a*(bp*dp-cp*cp)-invariant) == 0
    assert sp.expand(two(*transformed)-two(b,c,d,U)
                     -sp.Rational(11,155520)*a*a*(c*r-d*s+a*r*r/2)) == 0
    x,y,t=sp.symbols('x y t',nonzero=True)
    # On the nonsingular dense stratum the normalized invariant part is unique.
    H,tr,full,fourth=contractions((x*x+t*y*y)/2+y**3/6,(x,y))
    assert sp.simplify(H**2*(fourth-tr-2*full/3)) == -5/(3*t)
    assert sp.simplify(H**3*(fourth-tr-2*full/3)) == -sp.Rational(5,3)


if __name__=='__main__':
    check_composition()
    rank_one()
    print('Noncommuting normal changes, intrinsic projection and rank-one polynomial transport passed.')
