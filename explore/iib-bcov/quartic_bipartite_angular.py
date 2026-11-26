"""Positive Gegenbauer expansion of a bipartite four-loop period in even dimension."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import comb,factorial,prod

from quartic_projective_ibp import K33


@cache
def radial_weights(n,b,q):
    """Connection coefficients of the Riesz kernel on the (2q-1)-sphere."""
    lam=q-1
    if not 1<=b<=lam:raise ValueError('The two angular propagator powers must lie in [1,q-1].')
    value=F(comb(b+n-1,n),comb(lam+n-1,n))
    terms=[value]
    for k in range(lam-b):
        value*=F((b+n+k)*(b-lam+k),(k+1)*(n+lam+1+k))
        terms.append(value)
    return tuple(terms)


@cache
def triangle_mode(n,b0,b1,b2,q):
    """Below unit radius return the powers of r and log(r) and their rational coefficients."""
    dimension=2*(q-b0)
    result=defaultdict(F)
    measure=F(2*(q-1),factorial(q-1)*(n+q-1))
    for k1,c1 in enumerate(radial_weights(n,b1,q)):
        u=n+2*k1
        for k2,c2 in enumerate(radial_weights(n,b2,q)):
            v=n+2*k2
            c=measure*c1*c2
            low=dimension+u-2*b2
            delta=low-v
            result[low,0]+=c/F(dimension+u+v)
            result[v,0]+=c/F(u+v+2*(b1+b2)-dimension)
            if delta:
                result[v,0]+=c/F(delta)
                result[low,0]-=c/F(delta)
            else:
                result[v,1]-=c
    return tuple((r,h,c) for (r,h),c in sorted(result.items()) if c)


@cache
def angular_weight(i,j,k,q):
    """Mean of three Gegenbauer harmonics of degrees j+k, k+i, i+j."""
    lam=q-1
    s=i+j+k
    return (F(lam*factorial(lam-1),factorial(2*lam-1))
            *prod(s+d for d in range(lam+1,2*lam))
            *prod(comb(lam-1+t,t) for t in (i,j,k)))


def shell_term(i,j,k,columns,q):
    angular=angular_weight(i,j,k,q)*F(2,factorial(q-1))
    summand=F(0)
    for exterior in (False,True):
        radial_dimension=2*(sum(b[2] for b in columns)-q) if exterior else 2*q
        modes=[triangle_mode(n,2*q-sum(b) if exterior else b[0],b[1],b[2],q)
               for n,b in zip((j+k,k+i,i+j),columns)]
        for t in product(*modes):
            power=radial_dimension+sum(v[0] for v in t)
            logs=sum(v[1] for v in t)
            assert power>0
            summand+=F((-1)**logs*factorial(logs),power**(1+logs))*prod(v[2] for v in t)
    assert summand>0
    return angular*summand


def columns_and_scale(p,q,partition):
    a=(0,2,4);b=(1,3,5)
    if partition:a,b=b,a
    idx={edge:i for i,edge in enumerate(K33)}
    columns=tuple(tuple(q-p[idx[tuple(sorted((u,v)))]] for u in a) for v in b)
    scale=prod(F(factorial(q-1-a),factorial(a-1)) for a in p)
    return columns,scale


def through(p,q,last,partition=0):
    columns,scale=columns_and_scale(p,q,partition)
    shells=[sum(shell_term(i,j,s-i-j,columns,q) for i in range(s+1) for j in range(s-i+1))
            for s in range(last+1)]
    return scale*sum(shells)


def controls():
    import json
    from pathlib import Path
    q=3
    data=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())['K33']
    # Euler's integral gives strict rational bounds for both zeta constants.
    def zeta_bounds(power,cutoff=128):
        partial=sum(F(1,j**power) for j in range(1,cutoff+1))
        return (partial+F(1,(power-1)*(cutoff+1)**(power-1)),
                partial+F(1,(power-1)*cutoff**(power-1)))
    z3=zeta_bounds(3);z5=zeta_bounds(5)
    for row in data:
        p=row['powers']
        r,b,c=map(F,row['1,zeta3,zeta5'])
        lower=r+min(b*z for z in z3)+min(c*z for z in z5)
        norm=F(prod(factorial(a-1) for a in p),factorial(q-1))
        for order in (0,1):
            values=[norm*through(p,q,n,order) for n in (0,2,6)]
            assert 0<values[0]<values[1]<values[2]<lower
        print('positive six-shell projective lower sum',p,values[-1],flush=True)
    # With all angular propagators at q-1 only the constant Gegenbauer
    # connection coefficient remains. This checks the normalization of both
    # the convolution and the outer angular average against Gegenbauer norms.
    for q in range(2,8):
        for n in range(7):
            assert radial_weights(n,q-1,q)==(1,)
            assert angular_weight(n,0,0,q)==F(q-1,n+q-1)*comb(2*q-3+n,n)


if __name__=='__main__':
    controls()
