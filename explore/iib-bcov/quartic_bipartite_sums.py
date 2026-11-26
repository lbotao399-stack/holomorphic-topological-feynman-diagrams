"""Partial fractions and single-shell harmonic sums for the three triangle modes."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import comb,factorial,prod

from quartic_bipartite_angular import angular_weight,columns_and_scale,triangle_mode


@cache
def connection(q,b,k):
    lam=q-1
    roots=tuple(range(b+k,lam))+tuple(range(lam+1,lam+k+1))
    scale=F((-1)**k*comb(lam-b,k)*factorial(lam-1),factorial(b-1))
    return scale,roots


@cache
def partial_fractions(roots):
    """Separate reciprocal powers of n+a without ever forming a high-degree polynomial."""
    profile=Counter(roots)
    out=[]
    for a,p in profile.items():
        coeff=[F(1)]+[F(0)]*(p-1)
        for b,h in profile.items():
            if b==a:continue
            jets=[F((-1)**j*comb(h+j-1,j),(b-a)**(h+j)) for j in range(p)]
            coeff=[sum(coeff[i]*jets[j-i] for i in range(j+1)) for j in range(p)]
        out.extend((a,p-j,c) for j,c in enumerate(coeff) if c)
    return tuple(out)


@cache
def rational_triangle(b0,b1,b2,q):
    """Terms c*r^(n+2l)*log(r)^h/(n+a)^p on the inner unit interval."""
    lam=q-1
    terms=defaultdict(F)
    def put(l,h,roots,weight):
        for a,p,c in partial_fractions(tuple(sorted(roots))):terms[l,h,a,p]+=weight*c
    for k1,k2 in product(range(q-b1),range(q-b2)):
        c1,t1=connection(q,b1,k1)
        c2,t2=connection(q,b2,k2)
        roots=(*t1,*t2,lam)
        factor=F(2*lam,factorial(lam))*c1*c2
        low=q-b0+k1-b2
        put(low,0,(*roots,q-b0+k1+k2),factor/2)
        put(k2,0,(*roots,b1+b2-q+b0+k1+k2),factor/2)
        delta=2*(low-k2)
        if delta:
            put(k2,0,roots,factor/delta)
            put(low,0,roots,-factor/delta)
        else:put(k2,1,roots,-factor)
    return tuple((*key,c) for key,c in sorted(terms.items()) if c)


def polynomial_binomial(s,a,lam):
    # binom(s-n+lam-1,lam-1) in powers of x=n+a.
    coeff=[1]
    for j in range(1,lam):
        t=s+a+j
        coeff=[(coeff[k]*t if k<len(coeff) else 0)-(coeff[k-1] if k else 0)
               for k in range(len(coeff)+1)]
    return [F(c,factorial(lam-1)) for c in coeff]


@cache
def harmonic(p,a,s):
    """Sum (u+a)^(-p) for 0<=u<=s including polynomial exponents."""
    if s<0:return F(0)
    return harmonic(p,a,s-1)+F(s+a)**(-p)


@cache
def excluded_diagonal(p,q,a,b,t):
    """Sum (u+a)^(-p)*(t-u+b)^(-q) across the excluded triangle."""
    ab=t+a+b
    if p<=0:
        return sum(comb(-p,j)*(-1)**j*ab**(-p-j)*harmonic(q-j,b,t) for j in range(1-p))
    if q<=0:return excluded_diagonal(q,p,b,a,t)
    return (sum(F(comb(p+q-h-1,q-1),ab**(p+q-h))*harmonic(h,a,t) for h in range(1,p+1))
            +sum(F(comb(p+q-h-1,p-1),ab**(p+q-h))*harmonic(h,b,t) for h in range(1,q+1)))


@cache
def excluded_sum(p,q,a,b,s):
    if s<=0:return F(0)
    return excluded_sum(p,q,a,b,s-1)+excluded_diagonal(p,q,a,b,s-1)


@cache
def pair_moment(p,q,a,b,s):
    if (q,b)<(p,a):return pair_moment(q,p,b,a,s)
    return harmonic(p,a,s)*harmonic(q,b,s)-excluded_sum(p,q,a,b,s)


@cache
def triple_moment(s,axes):
    """Sum three shifted reciprocal powers on i+j+k=s by a finite recurrence."""
    if axes!=tuple(sorted(axes)):return triple_moment(s,tuple(sorted(axes)))
    (p,a),(q,b),(r,c)=axes
    total=2*s+a+b+c
    if p==0:return pair_moment(q,r,b,c,s)
    if p<0:
        return sum(comb(-p,j)*comb(-p-j,k)*(-1)**(j+k)*total**(-p-j-k)
                   *pair_moment(q-j,r-k,b,c,s)
                   for j in range(1-p) for k in range(1-p-j))
    return (triple_moment(s,((p-1,a),(q,b),(r,c)))
            +triple_moment(s,((p,a),(q-1,b),(r,c)))
            +triple_moment(s,((p,a),(q,b),(r-1,c))))/total


def shell(p,q,s,partition=0):
    columns,scale=columns_and_scale(p,q,partition)
    lam=q-1
    total=F(0)
    for exterior in (False,True):
        radial=sum(b[2] for b in columns)-q if exterior else q
        factors=[]
        for b in columns:
            entries=defaultdict(F)
            for l,h,a,ex,c in rational_triangle(2*q-sum(b) if exterior else b[0],b[1],b[2],q):
                for j,w in enumerate(polynomial_binomial(s,a,lam)):
                    entries[l,h,a,ex-j]+=c*w
            factors.append([(*key,c) for key,c in entries.items() if c])
        for rows in product(*factors):
            l,h=sum(v[0] for v in rows),sum(v[1] for v in rows)
            integral=F((-1)**h*factorial(h),(2*(s+radial+l))**(h+1))
            total+=integral*prod(v[4] for v in rows)*triple_moment(s,tuple((v[3],v[2]) for v in rows))
    # The single-index binomials in the three directional factors were
    # included in polynomial_binomial. Only the remaining s weight is needed.
    angle=F(lam*factorial(lam-1),factorial(2*lam-1))*prod(s+d for d in range(lam+1,2*lam))
    return scale*F(2,factorial(q-1))*angle*total


def controls():
    import json
    from pathlib import Path
    from quartic_bipartite_angular import shell_term
    trials=((2,2,1,3),(2,1,1,3),(2,2,2,3),(0,2,2,3),(5,6,4,9),(-2,8,6,9),(4,10,9,12))
    for b0,b1,b2,q in trials:
        for n in (0,1,5,9):
            actual=defaultdict(F)
            for l,h,a,p,c in rational_triangle(b0,b1,b2,q):
                actual[n+2*l,h]+=c*F(n+a)**(-p)
            assert {k:v for k,v in actual.items() if v}=={(r,h):v for r,h,v in triangle_mode(n,b0,b1,b2,q)}
    for s in range(6):
        for axes in (((2,1),(1,3),(3,2)),((0,2),(1,1),(-2,2)),((-2,3),(-1,4),(1,2))):
            direct=sum(prod(F(s-u+a)**(-p) for u,(p,a) in zip((i,j,s-i-j),axes))
                       for i in range(s+1) for j in range(s-i+1))
            assert triple_moment(s,axes)==direct
    data=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())['K33']
    p=data[0]['powers']
    columns,scale=columns_and_scale(p,3,0)
    for s in (0,1,3):
        exact=scale*sum(shell_term(i,j,s-i-j,columns,3) for i in range(s+1) for j in range(s-i+1))
        assert shell(p,3,s)==exact
    print('Partial fractions at q=3,9,12 and single harmonic shell sums agree with direct radial/angular modes.',flush=True)


if __name__=='__main__':
    controls()
