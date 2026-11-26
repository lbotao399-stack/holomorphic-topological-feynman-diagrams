"""Three normals on one doubled edge and two on the other, with every angular integral."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from math import factorial, prod

from flint import fmpq, fmpq_mpoly_ctx

from quartic_seven_plane_double import (
    Arm, Q, SCALE, decompose, graph_data, local_moment, orders, polynomials, symmetries,
)


MARKS=((0,0,0,7,7),(0,0,0,7,8),(0,0,1,7,7),(0,0,1,7,8))
WEIGHTED=(F(-17678261,4400827200000),F(-3279091,1320248160000),
          F(55782439,264049632000000),F(-430110173,1188223344000000))


class SeparatedArm(Arm):
    def polynomial(self,indices):
        i,j,k=indices
        if (i,j) not in self.products:
            self.products[i,j]=self.polys[i]*self.polys[j]
        return self.products[i,j]*self.polys[k]

    def measure(self,indices):
        if indices not in self.moments:
            finite=0
            rejected={}
            for powers,coefficient in self.polynomial(indices).terms():
                exp=tuple(map(int,powers))
                i,j=self.split
                a=i+sum(exp[:i])
                b=j+sum(exp[i:i+j])
                if min(Q-a,Q-b,a+b-Q)>0:
                    finite+=int(coefficient)*local_moment(self.split,exp)
                else:
                    rejected[exp]=int(coefficient)
            self.moments[indices]=finite,rejected
        return self.moments[indices]


def integrated(data,order,original=False,**kinematics):
    ds,voltage=polynomials(data,order,original,**kinematics)
    return integrate_polynomials(data,ds,voltage)


def integrate_polynomials(data,ds,voltage):
    _,_,_,adj,_,groups,sizes=data
    arms=[SeparatedArm(g,n) for g,n in zip((groups[0],groups[2]),sizes)]
    def parts(poly):
        return [(arms[0].register(l),arms[1].register(r),m,num) for l,r,m,num in decompose(poly,groups)]
    dparts=[parts(d) for d in ds]
    bad=({}, {})
    total=F()
    for i,j in combinations_with_replacement(range(4),2):
        if not dparts[i] or not dparts[j] or not adj[i][j]:
            continue
        grouped=defaultdict(F)
        for u,v,w in product(dparts[i],dparts[j],parts(adj[i][j]*voltage)):
            middle=tuple(map(sum,zip(u[2],v[2],w[2])))
            key=tuple((*sorted((u[h],v[h])),w[h]) for h in range(2))
            grouped[key,middle]+=(1+(i!=j))*u[3]*v[3]*w[3]
        for (key,middle),coefficient in grouped.items():
            if not coefficient:
                continue
            degrees=[sum(arm.degrees[k] for k in ix)+len(arm.group)-Q for arm,ix in zip(arms,key)]
            degrees.append(sum(middle)+len(middle))
            assert sum(degrees)==2*Q and all(0<a<Q for a in degrees)
            factors=[arm.measure(ix) for arm,ix in zip(arms,key)]
            # A forbidden left monomial is retained as its full polynomial on the
            # right (and conversely). Cancel before taking any improper integrals.
            for h,(_,rejected) in enumerate(factors):
                if not rejected:
                    continue
                other=arms[1-h]
                ctx=fmpq_mpoly_ctx.get(tuple('x'+str(k) for k in range(len(other.group))))
                complement=ctx.from_dict(other.polynomial(key[1-h]).to_dict())
                for exp,intcoeff in rejected.items():
                    c=coefficient*intcoeff
                    term=complement*fmpq(c.numerator,c.denominator)
                    index=middle,exp
                    bad[h][index]=bad[h].get(index,0)+term
            outer=F(prod(factorial(Q-1-d) for d in degrees),factorial(Q-1))
            beta=F(prod(map(factorial,middle)),factorial(degrees[2]-1))
            total+=coefficient*outer*beta*prod(x for x,_ in factors)/SCALE**2
    assert not any(bool(p) for table in bad for p in table.values())
    return total


def summed(marks):
    images={tuple(sorted(a[i][0] for i in marks)) for a in symmetries()}
    weight=F(len(images),len(symmetries()))
    data=graph_data(marks)
    value=sum(mult*integrated(data,order,order==tuple(range(5))) for order,mult in orders(marks))
    return value,weight*value


if __name__=='__main__':
    from sys import argv
    pick=list(map(int,argv[1:])) or list(range(4))
    out=F()
    for i in pick:
        raw,weighted=summed(MARKS[i])
        assert weighted==WEIGHTED[i]
        print('three/two shape',i,MARKS[i],raw,weighted,flush=True)
        out+=weighted
    if pick==list(range(4)):
        assert out==F(-3161360627,475289337600000)
        print('three/two total',out)
