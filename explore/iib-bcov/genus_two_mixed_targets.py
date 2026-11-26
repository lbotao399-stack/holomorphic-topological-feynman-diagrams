"""Genus-two threefold vertices by Newton moments and independent fivefold Ward graphs."""
from fractions import Fraction as F
from scalar_ancestor_response import fused_totals, scalar_power
from genus_two_graph_ward import vertex


def check():
    charges=((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))
    assert scalar_power(2,6,charges,(1,0,0,0))==F(1,9720)
    charges=(*charges[:3],(1,-1,0),(-2,0,-1))
    hard=tuple((*p,0,0) for p in charges)
    forgotten=tuple(tuple(sign*int(j==i) for j in range(5)) for i in (3,4) for sign in (1,-1))
    for levels,wanted in (((2,0,0,0,0),F(569,77760)),((1,1,0,0,0),F(739,77760)),
                          ((0,0,0,2,0),F(37,8640))):
        weights=fused_totals(vertex,hard,levels,forgotten)
        value=scalar_power(2,6,charges,levels)
        assert 4*value==sum((-1)**k*v for k,v in enumerate(weights))
        assert value==wanted
        print(levels,weights,'gives',value,flush=True)


if __name__=='__main__':
    check()
