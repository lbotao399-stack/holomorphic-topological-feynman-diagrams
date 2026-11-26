"""The doubly marked two-loop fan on the one-loop triangle and a graviton."""
from math import factorial
from fractions import Fraction as F
from itertools import product
from functools import cache
from collections import defaultdict
import ternary_product as p

@cache
def moments(powers):
    a,b,c,d=(powers.count(i) for i in range(4))
    exact=F(factorial(b)*factorial(d),2*(a+1)*(c+1)*factorial(b+d+3))
    # Independent simplex integral after r4=(1-r2)*v. All factors have
    # degree <=3 on each interval for the single displaced input derivative.
    def simpson(f):
        return (f(F(0))+4*f(F(1,2))+f(F(1)))/6
    check=simpson(lambda r:r**a)*simpson(lambda r:r**c)
    check*=simpson(lambda r:r**b*(1-r)**(d+2))*simpson(lambda v:v**d*(1-v))/2
    assert exact==check
    return exact

def components(poly):
    out={}
    for v,y in product(range(3),repeat=2):
        for j2,partial in p.variational(poly,v,1).items():
            for j1,c in p.variational(partial,y).items():
                out[y,v,j1,j2]=c
    return out

def two_loop_three_ghosts(A,B):
    def one(color,indices):
        return {(p.letter('ghost',0,color,indices),):1}
    a,b=components(A),components(B)
    t=tuple(map(p.theta,range(6)))
    total=defaultdict(F)
    # (edge, differentiated site 1/2/3, factor, r-indices in monomial).
    shifts=(((0,1,(0,)),(0,-1,(0,1)),(1,-1,(1,2)),(2,-1,(1,))),
            ((0,1,()),(0,-1,(0,)),(0,1,(0,1)),(1,1,(1,2)),(2,1,(1,))),
            ((1,1,(2,)),(0,-1,(0,3)),(1,-1,(2,3)),(2,-1,(3,))),
            ((1,1,()),(1,-1,(2,)),(0,1,(0,3)),(1,1,(2,3)),(2,1,(3,))))
    for x,y,z,u1 in p.COLORS:
        for h,v,w,u2 in p.COLORS:
            for (yy,vv,j13,j23),left in a.items():
                if (y,v)!=(yy,vv) or not left: continue
                for (zz,ww,j14,j24),right in b.items():
                    if (z,w)!=(zz,ww) or not right: continue
                    # Only one external, non-interaction ghost may remain.
                    al={k:c for k,c in left.items() if all(l[0] in ('0theta','ghost') for l in k)}
                    br={k:c for k,c in right.items() if all(l[0] in ('0theta','ghost') for l in k)}
                    if not al or not br: continue
                    jet=tuple((edge,d) for edge,jj in enumerate((j13,j14,j23,j24)) for d in jj)
                    assert len(jet)<=1
                    for aa,bb,e in p.SPIN:
                        for dd,hh,f in p.SPIN:
                            for first,second in product(range(2),repeat=2):
                                initial=[[aa],[dd],[]]
                                initial[1 if first else 2].append(bb)
                                initial[0 if second else 2].append(hh)
                                pref=(2,)*first+(0,)*second
                                for assignment in product(*(shifts[edge] for edge,d in jet)):
                                    sites=[i.copy() for i in initial]
                                    signed=u1*u2*e*f
                                    weights=pref
                                    for (edge,d),(site,s,rr) in zip(jet,assignment):
                                        sites[site].append(d)
                                        signed*=s
                                        weights+=rr
                                    weight=-signed*moments(weights)
                                    raw=p.mul(one(x,sites[0]),one(h,sites[1]),p.with_derivatives(al,sites[2]),br)
                                    for word,coeff in raw.items():
                                        if word[:6]==t:
                                            total[word[6:]]+=weight*coeff
    return {w:v for w,v in total.items() if v}

def scalar():
    out=defaultdict(int)
    q=defaultdict(int)
    for (a,x,e),(b,y,f),(c,z,g) in product(p.SPIN,repeat=3):
        for i,j,k,s in p.COLORS:
            p.add(out,(p.letter('ghost',0,i,[a]),p.letter('ghost',0,j,[b,c]),
                       p.letter('ghost',0,k,[x,y,z])),s*e*f*g)
        for i in range(3):
            p.add(q,(p.letter('ghost',0,i,[a,b,c]),p.letter('ghost',0,i,[x,y,z])),e*f*g)
    assert not p.plus((p.classical(q),1),(out,-6))
    return out,q

def main():
    g,h,k=map(p.graviton,range(3))
    scalar0,q=scalar()
    for a,b,expected in [(p.triangle(g,h),k,F(1,24)),(g,p.triangle(h,k),F(1,72))]:
        o3=two_loop_three_ghosts(a,b)
        assert not p.plus((o3,1),(scalar0,-expected))
        assert not p.classical(o3)
        print('Coefficient in g^3/(32 pi^2)^3 units with trace factors:',expected/8)
    assert not p.first_fan_displaced(q)
if __name__=='__main__':
    main()
