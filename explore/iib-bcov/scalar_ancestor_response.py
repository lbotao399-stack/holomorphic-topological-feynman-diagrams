"""Scalar-target intersection measurements of the fivefold elementary differential."""
from collections import Counter
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import sys

from elliptic_descendants import vertex as elliptic_vertex
from genus_two_graph_ward import vertex as two_loop_vertex

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'twisted-holography-N4'/'code'))
from one_target_dr import homogeneous_dr
from betti_power_graphs import scalar_power


def projected(genus,scalar_charges,levels):
    """Four groups of g unit transverse marks, polarized and labeled."""
    if sum(scalar_charges) or sum(levels)!=2*genus-3+len(levels) or max(levels)>4:
        raise ValueError('Check axial momentum conservation and the fivefold descendant levels.')
    coefficient=homogeneous_dr(genus,tuple(scalar_charges)).to_dict().get(tuple(levels),0)
    return factorial(genus)**4*F(str(coefficient))


def fused_totals(vertex,retained,levels,forgotten):
    """Forget labeled primary waves successively, lowering the chosen surviving mark."""
    totals=Counter()
    def visit(waves,ks,left,absorbed):
        if not left:
            totals[absorbed]+=vertex(waves,ks)
            return
        p,*rest=left
        visit((*waves,p),(*ks,0),rest,absorbed)
        for i,k in enumerate(ks):
            if not k: continue
            current=waves[:i]+(tuple(x+y for x,y in zip(waves[i],p)),)+waves[i+1:]
            visit(current,ks[:i]+(k-1,)+ks[i+1:],rest,absorbed+1)
    visit(tuple(retained),tuple(levels),forgotten,0)
    return [totals[j] for j in range(max(totals)+1)]


def check():
    # Compact intersections and independent fivefold Feynman vertices.
    first=[(1,0,0,0,0),(0,1,0,0,0),(0,0,1,0,0),(1,2,1,0,0),(-2,-3,-2,-1,-1)]
    extra=((0,0,0,1,0),(0,0,0,0,1))
    levels=(1,1,0,0,0)
    a=fused_totals(elliptic_vertex,first,levels,extra)
    assert sum((-1)**i*x for i,x in enumerate(a))==elliptic_vertex([p[:3] for p in first],levels)
    print('Two distinct elliptic descendant marks:',a,flush=True)

    extra=tuple(tuple(sign*int(i==j) for j in range(5)) for i in (2,3,4) for sign in (1,-1))
    cases=(((1,0),(0,1),(-1,-1)),((1,0),(0,1),(1,-1),(-2,0)))
    for charges,levels in zip(cases,((2,0,0),(1,2,0,0))):
        points=tuple((*p,0,0,0) for p in charges)
        a=fused_totals(two_loop_vertex,points,levels,extra)
        base=scalar_power(2,4,charges,levels)
        assert sum((-1)**i*x for i,x in enumerate(a))==8*base
        print('Genus-two fivefold graph totals by absorbed count:',a,'target:',8*base,flush=True)
    assert projected(2,(1,1,-2),(4,0,0))==F(91,360)
    expected=(F(5,24),F(91,360),F(985827,1120),F(19064476708,175))
    for g,want in enumerate(expected,1):
        m=(2*g+1)//3
        p=(*range(1,m+1),1,-(1+m*(m+1)//2))
        k=(*(4,)*(m-1),2*g+3-3*m,0,0)
        value=projected(g,p,k)
        assert value==want
        assert value==factorial(g)**4*scalar_power(g,g,p,k)
        print('Scalar ancestor, genus',g,'levels',k,'charges',p,'value',value,flush=True)


if __name__=='__main__':
    check()
