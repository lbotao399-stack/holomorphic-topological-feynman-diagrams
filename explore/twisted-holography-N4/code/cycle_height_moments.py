"""Full edge-support Newton moments of a cycle with arbitrary conserved vector charges."""
from itertools import combinations
from functools import lru_cache
from fractions import Fraction
from math import factorial,prod

import sympy as sp


def simplex(polynomial, coordinates):
    """Integrate over t_i >= 0, sum(t_i) <= 1."""
    terms = sp.Poly(polynomial, *coordinates).terms()
    return sum(coeff*sp.prod(factorial(k) for k in powers)
               /factorial(len(coordinates)+sum(powers)) for powers,coeff in terms)


def merge_parallel(powers):
    """Effective exponent and raw moment factor for parallel Euclidean capacities."""
    if not powers or any(type(k) is not int or k<0 for k in powers):
        raise ValueError('Supply the nonnegative powers on a nonempty parallel family.')
    effective = sum(powers)+len(powers)-1
    factor = sp.Rational(2**(len(powers)-1)*sp.prod(factorial(2*k+1) for k in powers),
                         factorial(2*effective+1))
    return effective, factor


def cycle_moment(gram, powers, root=0):
    """Coefficient moment before the exponential Taylor signs and factorials.

    gram is the Gram matrix of the cumulative edge charge offsets in the cycle.
    The output is the signed full-support integral against product(a_i**powers[i]).
    """
    gram = sp.Matrix(gram)
    count = len(powers)
    if gram.shape!=(count,count) or gram!=gram.T or not 0<=root<count:
        raise ValueError('Supply a symmetric Gram matrix on the edges and a valid root.')
    if any(type(k) is not int or k<0 for k in powers):
        raise ValueError('All powers must be nonnegative integers.')
    others = tuple(i for i in range(count) if i!=root)
    # Relative Gram entries are invariant under a common translation of all offsets.
    relative = sp.Matrix(count,count,lambda i,j: gram[i,j]-gram[i,root]-gram[root,j]+gram[root,root])
    total = 0
    for size in range(1,count):
        for active in combinations(others,size):
            jacobian = relative.extract(active,active).det()
            if jacobian==0:
                continue
            t = sp.symbols('t:'+str(size))
            norm = sum(t[j]*t[k]*relative[v,w] for j,v in enumerate(active) for k,w in enumerate(active))
            height = [(norm-2*sum(z*relative[v,i] for z,v in zip(t,active))+relative[i,i])/2
                      for i in range(count)]
            integrand = (1-sum(t))*height[root]**(powers[root]+1)/(powers[root]+1)
            for i in others:
                if i in active:
                    integrand *= height[i]**powers[i]
                else:
                    integrand *= -height[i]**(powers[i]+1)/(powers[i]+1)
            total -= jacobian*simplex(integrand,t)
    return sp.factor(total)


@lru_cache(maxsize=32)
def cycle_polynomial(powers):
    """Compile one universal polynomial in the Gram entries relative to the first offset."""
    size = len(powers)
    gram = sp.zeros(size)
    entries = []
    for i in range(1,size):
        for j in range(1,i+1):
            symbol = sp.Symbol('g%d%d'%(i,j))
            entries.append(symbol)
            gram[i,j] = gram[j,i] = symbol
    polynomial = sp.Poly(cycle_moment(gram,powers), *entries)
    denominator, integer = polynomial.clear_denoms()
    return int(denominator),tuple((int(c),tuple((i,e) for i,e in enumerate(m) if e))
                                  for m,c in integer.terms())


@lru_cache(maxsize=8192)
def numeric_cycle(gram,powers):
    """Evaluate on actual Gram entries, using only simplex variables in the polynomial ring."""
    from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx
    n=len(powers)
    g=tuple(tuple(fmpq(str(x)) for x in row) for row in gram)
    total=fmpq(0)
    for size in range(1,n):
        ctx=fmpq_mpoly_ctx.get(tuple('u'+str(i) for i in range(size)))
        t=ctx.gens()
        for active in combinations(range(1,n),size):
            jac=fmpq_mat([[g[i][j] for j in active] for i in active]).det()
            if not jac: continue
            norm=sum(g[i][j]*x*y for i,x in zip(active,t) for j,y in zip(active,t))
            h=[(norm+g[i][i])/2-sum(x*g[i][j] for j,x in zip(active,t)) for i in range(n)]
            integrand=(1-sum(t))*h[0]**(powers[0]+1)/(powers[0]+1)
            for i in range(1,n):
                if i in active: integrand*=h[i]**powers[i]
                else: integrand*=-h[i]**(powers[i]+1)/(powers[i]+1)
            total-=jac*sum(c*prod(factorial(k) for k in m)/factorial(size+sum(m))
                           for m,c in integrand.terms())
    return Fraction(str(total))


def cycle_value(offsets, powers):
    vectors = [tuple(y-x for x,y in zip(offsets[0],p)) for p in offsets[1:]]
    if len(powers)+sum(powers)>=5:
        full=[(0,)*len(offsets[0]),*vectors]
        gram=tuple(tuple(sum(a*b for a,b in zip(v,w)) for w in full) for v in full)
        return numeric_cycle(gram,powers)
    values = [sum(a*b for a,b in zip(v,w)) for i,v in enumerate(vectors) for w in vectors[:i+1]]
    denominator,terms = cycle_polynomial(powers)
    from math import prod
    return Fraction(sum(c*prod(values[i]**e for i,e in m) for c,m in terms),denominator)


def check():
    a,b,c = sp.symbols('a b c')
    d = (b+c-a)/2
    triangle = sp.Matrix([[0,0,0],[0,c,d],[0,d,b]])
    expected = -(a+b+c)**3/960+(a+b+c)*(a*b+b*c+c*a)/120-7*a*b*c/240
    assert all(sp.expand(cycle_moment(triangle,(0,0,0),root)-expected)==0 for root in range(3))
    assert expected.subs({a:1,b:1,c:1})==sp.Rational(17,960)
    # One-dimensional cut boxes on two paths in series give an independent
    # formula: the minimum edge energy on either path is exponentially distributed.
    for first,second in ((1,1),(1,2),(1,3),(2,2),(2,3)):
        offsets = sp.Matrix([0]*first+[1]*second)
        count = first+second
        direct = (-1)**(count-1)*sp.Rational(first*second*factorial(2*first-1)*factorial(2*second-1),
                                            2**(count-2)*factorial(2*count))
        assert cycle_moment(offsets*offsets.T,(0,)*count)==direct
    parallel = sp.diag(0,1)
    for j,k in ((0,1),(2,0),(1,2),(2,2)):
        # Integrate the failed cut rho_1+rho_2<1 directly.
        direct = -sp.Rational(factorial(2*j+1)*factorial(2*k+1),
                              2**(j+k)*factorial(2*j+2*k+4))
        assert cycle_moment(parallel,(j,k))==direct
    for powers in ((0,0),(0,0,0),(2,1,0),(0,1,1,2)):
        effective,factor = merge_parallel(powers)
        direct = -sp.Rational(sp.prod(factorial(2*k+1) for k in powers),
                              2**sum(powers)*factorial(2*sum(powers)+2*len(powers)))
        assert -factor/(2**(effective+1)*(effective+1))==direct
    # Two direct paths parallel to one two-edge series path with unit charge.
    # The two cut failures share both direct paths. The full-support moment is
    # the positive integral over their intersection, in rho coordinates.
    r = sp.Symbol('r')
    direct = sp.integrate(r**3*(1-r)**4/24,(r,0,1))
    effective,factor = merge_parallel((0,0))
    offsets = sp.Matrix([0,1,1])
    assert factor*cycle_moment(offsets*offsets.T,(effective,0,0))==direct==sp.Rational(1,6720)
    for points,powers in (([(1,2),(0,1),(2,2),(1,0)],(1,0,1,0)),
                          ([(0,0,1),(1,0,1),(0,1,2),(1,-1,0)],(0,1,0,0))):
        coords = sp.Matrix(points)
        gram = coords*coords.T
        assert len({cycle_moment(gram,powers,root) for root in range(4)})==1
    print('Mixed triangle:',sp.factor(expected))
    print('All cycle root exchanges, weighted parallel cuts and charge-free series substitutions passed.')


if __name__=='__main__':
    check()
