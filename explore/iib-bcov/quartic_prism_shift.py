"""Symmetric determinant dimension lowering on the prism hyperplane L_14=0."""
from fractions import Fraction as F
from gzip import open as open_gzip
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
import json

from flint import fmpq, fmpq_mat
from quartic_planar_periods import core_data
from quartic_projective_ibp import PRISM, PERMS
from square_three import G


def symbols():
    cycle=G(6,PRISM).C
    ctx,u,_,_=core_data(tuple(PRISM))
    x=ctx.gens()
    slots=[(i,j) for i in range(4) for j in range(i,4) if (i,j)!=(0,3)]
    assert all(cycle[0][e]*cycle[3][e]==0 for e in range(9))
    matrix=fmpq_mat([[cycle[i][e]*cycle[j][e]*(1+(i!=j)) for i,j in slots] for e in range(9)])
    inverse=matrix.inv()
    ops=[[ctx.constant(0) for _ in range(4)] for _ in range(4)]
    for k,(i,j) in enumerate(slots):
        ops[i][j]=ops[j][i]=sum(inverse[k,e]*x[e] for e in range(9))
    def minor(ids):
        return sum((-1)**sum(a>b for i,a in enumerate(p) for b in p[i+1:])
                   *prod(ops[i][j] for i,j in zip(ids,p)) for p in permutations(ids))
    a=sum(x[e]*cycle[0][e]**2 for e in range(9))
    b=sum(x[e]*cycle[3][e]**2 for e in range(9))
    d=minor((1,2))
    h=minor((0,1,2))
    k=minor((1,2,3))
    assert a==x[1]+x[2]+x[5] and b==x[3]+x[4]+x[8]
    r,s,t=x[0],x[6],x[7]
    assert 4*d==-r*r-s*s-t*t+2*(r*s+r*t+s*t)
    return u,a,b,d,h,k


def taylor_derivative(u,operator,strength,times):
    """Return operator(U^strength)/U^strength by an independent local Taylor recurrence."""
    expansion={}
    for ex,c in u.terms():
        support=[i for i,k in enumerate(ex) if k]
        for length in range(len(support)+1):
            for subset in combinations(support,length):
                mask=sum(1<<i for i in subset)
                expansion[mask]=expansion.get(mask,0)+int(c)*prod(times[i] for i in support if i not in subset)
    @cache
    def coefficient(ex):
        if not any(ex):return F(1)
        i=next(i for i,a in enumerate(ex) if a)
        value=F(0)
        mask=sum(1<<j for j,a in enumerate(ex) if a)
        part=mask
        while part:
            if part in expansion:
                beta=tuple((part>>j)&1 for j in range(9))
                child=tuple(v-k for v,k in zip(ex,beta))
                value+=((strength+1)*beta[i]-ex[i])*expansion[part]*coefficient(child)
            part=(part-1)&mask
        return value/(ex[i]*expansion[0])
    return sum(F(str(c))*prod(factorial(int(k)) for k in ex)*coefficient(tuple(map(int,ex)))
               for ex,c in operator.terms())


def falling(k,r):
    return prod(range(k-r+1,k+1)) if k>=r else 0


def prism_twelve(u,d,h,k):
    ex=(5,4,4,4,4,4,5,5,6)
    cotrees=[tuple(map(int,t)) for t,c in u.terms()]
    # Every strict subset margin follows from a spanning forest of that
    # subset extended to a global tree: h(S)=min_T |T intersect S|.
    subsets=[(ids,min(sum(t[i] for i in ids) for t in cotrees)) for n in range(1,9)
             for ids in combinations(range(9),n)]
    required=set()
    sizes=[]
    for op,q in ((d,12),(h*k,11)):
        degree=tuple(map(int,op.degrees()))
        assert all(a>=b for a,b in zip(ex,degree))
        minimum=100
        for r,c in op.terms():
            powers=tuple(a-int(b)+1 for a,b in zip(ex,r))
            assert sum(powers)==4*q
            minimum=min(minimum,min(sum(powers[i] for i in ids)-q*cycles for ids,cycles in subsets))
            if q==12:
                required.add(min(tuple(powers[j] for j in p) for p in PERMS))
        assert minimum>0
        sizes.append((len(op),minimum))
    present=set()
    with open_gzip(Path(__file__).with_name('quartic_prism.weights.gz'),'rt') as bank:
        for line in bank:
            *p,kweight,jweight=line.split()
            powers=tuple(map(int,p))
            if powers in required and kweight!='0' and jweight!='0':
                present.add(powers)
    print('D2/D6 monomials and least margins',sizes,'physical high-dimensional orbits',len(present),
          'of',len(required),flush=True)
    return sizes,len(required),len(present)


def integral_rows(data,u,factor,operator,prefactor,degree):
    """Check integration by parts, with all original coordinate-boundary terms zero."""
    values=[tuple(map(F,row)) for row in data['coefficients']]
    lookup={tuple(row[j] for j in p):values[i] for i,row in enumerate(data['exponents']) for p in PERMS}
    ctx=u.context()
    terms=[(tuple(map(int,ex)),c) for ex,c in operator.terms()]
    minimum=tuple(map(int,operator.degrees()))
    ranges=[range(m,3) for m in minimum]
    assert all(ranges)
    tested=0
    for ex in product(*ranges):
        if sum(ex)!=degree:continue
        h=ctx.from_dict({ex:1})
        derivative=ctx.from_dict({tuple(v-k for v,k in zip(ex,r)):c*prod(falling(v,k) for v,k in zip(ex,r))
                                  for r,c in terms})
        left=operator.total_degree()
        first=prefactor*h*factor
        second=(-1)**int(left)*u*derivative
        # For each monomial derivative of h require every cotree child to
        # converge separately, even if some cancel algebraically in U*D(h).
        if any(tuple(v-k+int(t) for v,k,t in zip(ex,r,tree)) not in lookup
               for r,c in terms for tree,w in u.terms()):continue
        if any(tuple(map(int,r)) not in lookup for r,c in first.terms()):continue
        for i in range(3):
            def integrate(poly):
                return sum(F(str(c))*lookup[tuple(map(int,r))][i] for r,c in poly.terms())
            assert integrate(first)==integrate(second)
        tested+=1
    return tested


def certificate():
    u,a,b,d,h,k=symbols()
    for s in (-12,-6,2):
        for times in ([1]*9,[2,3,5,7,11,13,17,19,23]):
            tree=F(str(u(*times)))
            def apply(op,power):
                return taylor_derivative(u,op,power,times)
            lower=(s+1)*(F(s)+F(3,2))
            assert tree*apply(d,s+1)==lower*F(str((a*b)(*times)))
            assert tree*apply(h,s+1)==lower*(s+2)*F(str(b(*times)))
            assert tree*apply(k,s+1)==lower*(s+2)*F(str(a(*times)))
            assert tree*apply(h*k,s+1)==lower*(s+2)**2*apply(d,s)
    bank=Path(__file__).resolve().parents[1]/'nonSP/verify/prism_forcer/periods.json'
    data=json.loads(bank.read_text())
    q=6
    c2=fmpq(1-q)*fmpq(3-2*q,2)
    c3=c2*(2-q)
    counts=[integral_rows(data,u,a*b,d,c2,13),
            integral_rows(data,u,b,h,c3,14),integral_rows(data,u,a,k,c3,14)]
    assert counts==[84,0,0]
    print('Two-row shifts checked on the integrated q=6 bank:',counts[0],flush=True)
    print('Full sixth/second derivative identity checked at q=12, q=6, and positive exponent 2.',flush=True)
    assert prism_twelve(u,d,h,k)==([(6,3),(328,2)],4,4)


if __name__=='__main__':
    certificate()
