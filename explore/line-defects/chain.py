"""Fixed endpoints on a holomorphic line, two integrated bulk vertices in C^3."""
from functools import cache
from fractions import Fraction as F
from math import comb, factorial, prod
import sympy as sp

u,v,s = sp.symbols('u v s')
w,z = sp.symbols('w z')
variables = (u,v,s)
xx,xy,yy = u*(v+s),u*s,s*(u+v)


def grounded(incidence,conductances,fixed,positions):
    """Return harmonic extension, edge currents, covariance and boundary quadratic form."""
    B = sp.Matrix(incidence)
    weighted = sp.diag(*conductances)
    Bf = B[:,[i for i in range(B.cols) if i not in fixed]]
    Bd = B[:,fixed]
    laplace = Bf.T*weighted*Bf
    inverse = laplace.inv()
    transfer = -inverse*Bf.T*weighted*Bd
    boundary = Bd.T*weighted*(Bd+Bf*transfer)
    return transfer*positions,weighted*(Bd+Bf*transfer)*positions,inverse,boundary


@cache
def wick(hx,hy,bx,by):
    """One transverse complex coordinate; h and b powers denote holomorphic and barred factors."""
    if hx+hy!=bx+by:
        return 0
    terms = []
    for j in range(max(0,bx-hy),min(hx,bx)+1):
        # j X-X contractions, bx-j Y-X contractions.
        count = comb(hx,j)*comb(hy,bx-j)*factorial(bx)*factorial(by)
        terms.append(count*xx**j*xy**(hx+bx-2*j)*yy**(hy-bx+j))
    return sum(terms)


def simplex(polynomial):
    total = 0
    for powers,c in sp.Poly(polynomial,*variables).terms():
        total += c*sp.Rational(sp.prod(factorial(n) for n in powers),factorial(sum(powers)+2))
    return total


def integrated(hol,bar,denominators,linear_x=False):
    """Integral divided by (-2 pi i)^6, times the appropriate power of |w-z|^2.

    Input powers have order (x1,x2,y1,y2). Include the transverse barred determinant.
    linear_x inserts x^0. Return its coefficient and the exponent of |w-z|^{-2}.
    """
    r,t,m,n = hol
    k,h,l,j = bar
    a,b,c = denominators
    left = wick(r,m,k+1,l)*wick(t,n,h,j+1)
    right = wick(r,m,k,l+1)*wick(t,n,h+1,j)
    if left==right:
        return 0, a+b+c-6-sum(hol)
    radial = a+b+c-6-sum(hol)
    if radial<=0:
        raise ValueError('the radial integral diverges for this numerator')
    # xx*yy - xy^2 = u*v*s*(u+v+s). Keep the simplex variables homogeneous.
    density=sp.cancel((left-right)*u**(2-a)*v**(2-b)*s**(2-c))
    if not density.is_polynomial(*variables):
        raise ValueError('a separate boundary convergence check is needed')
    if linear_x:
        density *= (v+s)*w+u*z
    constant = sp.Rational(factorial(radial-1),factorial(a-1)*factorial(b-1)*factorial(c-1))
    return sp.expand(constant*simplex(density)),radial


@cache
def pairings(hx,hy,bx,by):
    """Multiplicities of xx,xy,yy contractions in a single coordinate."""
    if min(hx,hy,bx,by)<0 or hx+hy!=bx+by:
        return ()
    return tuple((j,hx+bx-2*j,hy-bx+j,
                  comb(hx,j)*comb(hy,bx-j)*factorial(bx)*factorial(by))
                 for j in range(max(0,bx-hy),min(hx,bx)+1))


@cache
def bridge_moment(i,j,k,bar_x,bar_y,p,q):
    # After det(cov) has removed one u,v,s factor, expand the covariance
    # and the longitudinal means of (x^0)^p (y^0)^q on u+v+s=1.
    a=i+j-bar_x
    c=k+j-bar_y
    assert min(a,c)>=0
    values=[F(0)]*(p+q+1)
    for m in range(p+1):   # x^0 = w*(v+s)+z*u
        for n in range(q+1):  # y^0 = w*s+z*(u+v)
            left=i+p-m
            right=k+n
            factor=comb(p,m)*comb(q,n)
            for h in range(left+1):
                for l in range(right+1):
                    powers=(a+m+l,left-h+right-l,c+q-n+h)
                    values[m+n]+=F(factor*comb(left,h)*comb(right,l)*prod(map(factorial,powers)),
                                    factorial(sum(powers)+2))
    return tuple(values)


def gauge_box(hol,bar,p=1,q=0):
    """Coefficients of w^(p+q-j) z^j in pi^2 (w-z) M/dbar(z0), Fernandez C.2."""
    if min(*hol,*bar,p,q)<0:
        raise ValueError('monomial powers must be nonnegative')
    r,t,m,n=hol
    k,h,l,j=bar
    values=[F(0)]*(p+q+1)
    for sign,a,b,c,d in ((r*n,r-1,t,m,n-1),(-t*m,r,t-1,m-1,n)):
        if not sign:
            continue
        for i0,j0,k0,v0 in pairings(a,c,k,l):
            for i1,j1,k1,v1 in pairings(b,d,h,j):
                row=bridge_moment(i0+i1,j0+j1,k0+k1,k+h,l+j,p,q)
                for e,coefficient in enumerate(row):
                    values[e]+=sign*v0*v1*coefficient
    # (2+k+h)! and (2+l+j)! from the vertex factors cancel the two
    # endpoint Schwinger factorials, leaving the central 2! and the form prefactor.
    return tuple(v/F(-16*prod(map(factorial,bar))) for v in values)


def check():
    b=[[-1,1,0,0],[0,-1,1,0],[0,0,-1,1]]
    mean,current,cov,boundary=grounded(b,[1/u,1/v,1/s],[0,3],sp.Matrix([w,z]))
    total=u+v+s
    assert all(sp.cancel(x)==0 for x in (current-sp.ones(3,1)*(z-w)/total))
    assert sp.simplify(cov-sp.Matrix([[xx,xy],[xy,yy]])/total)==sp.zeros(2)
    assert sp.simplify(mean-sp.Matrix([((v+s)*w+u*z)/total,(s*w+(u+v)*z)/total]))==sp.zeros(2,1)
    assert sp.simplify(total*boundary-sp.Matrix([[1,-1],[-1,1]]))==sp.zeros(2)

    value,radial=integrated((0,1,1,0),(0,0,0,0),(3,3,3),True)
    assert radial==1 and sp.expand(value+(2*w+z)/48)==0
    # Fernandez (7.10), with the differential dbar(z0) stripped off.
    pref=-sp.Rational(1,4)*(2*sp.pi*sp.I)**-2*(4*sp.pi**2)**-3*8
    assert sp.simplify(pref*(-2*sp.pi*sp.I)**6*value-(2*w+z)/(96*sp.pi**2))==0

    tests=[((3,0,0,3),(1,1,1,1),(5,3,5),36,-(22*w+13*z)/35,True),
           ((3,0,0,3),(0,2,2,0),(5,3,5),36,-2*(4*w+3*z)/35,True),
           ((2,0,0,2),(0,1,1,0),(4,3,4),sp.Rational(81,64),-sp.Rational(3,16),False)]
    for hol,bar,den,pref,target,linear in tests:
        value,radial=integrated(hol,bar,den,linear)
        assert radial==1 and sp.simplify(pref*(-2*sp.I)**6*value-target)==0,(value,target)
    assert gauge_box((0,1,1,0),(0,0,0,0))==(F(1,48),F(1,96))
    assert gauge_box((0,1,1,0),(0,0,0,0),0,0)==(F(1,32),)
    assert gauge_box((0,1,1,0),(0,0,0,0),0,1)==(F(1,96),F(1,48))
    print('grounded covariance, gauge box, three gravity boxes passed')


if __name__=='__main__':
    check()
