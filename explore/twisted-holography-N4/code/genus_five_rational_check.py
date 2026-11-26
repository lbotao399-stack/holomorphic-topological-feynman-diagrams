"""Original projected area charts for high-edge (2,2,1) stable-graph moments."""
from fractions import Fraction

from betti_series_parallel import coefficient
from genus_five_rational import LEVELS, capacity_moment, reduced_input
from planar_four_loop import stencil
from scalar_flow_recursion import coefficient as scalar_coefficient


# Each root is an edge index in the original graph, before removing any vertices.
CASES=(
    (((0,1),(0,2),(1,3),(1,4),(2,3),(2,4),(3,4)),
     ((1,0),(0,1),(-1,-1),(0,0),(0,0)),(0,1,0,0,0,0,0),
     6,'1/4608','30203431/73229764608000'),
    (((0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(2,3)),
     ((-1,-1),(0,0),(1,0),(0,1)),(0,0,0,0,0,1,2),
     5,'2/3','-11238390433/487066982036963328000'),
    (((0,1),(0,2),(1,3),(1,4),(2,3),(2,4),(3,4),(3,4)),
     ((0,0),(-1,-1),(0,0),(1,0),(0,1)),(1,1,0,0,0,0,0,0),
     0,'1/48','256200495893/13113341824072089600'),
    (((0,1),(0,1),(0,2),(1,2),(2,3),(2,4),(3,4),(3,4),(3,4)),
     ((-1,-1),(0,0),(0,0),(1,0),(0,1)),(0,0,0,0,0,1,0,0,0),
     2,'1/12','5219/134120448000'),
    (((0,1),(0,1),(0,2),(1,2),(2,4),(3,4),(3,5),(3,5),(4,5),(4,5)),
     ((-1,-1),(0,0),(0,0),(0,0),(0,1),(1,0)),(0,)*10,
     4,'1/8','-1/40642560'),
)


def check():
    for edges,charges,powers,root,weight,answer in CASES:
        table=stencil(5,len(edges),(LEVELS[-1],),shape_subset=(edges,))
        assert table[edges,charges,powers,0]==Fraction(weight)
        assert capacity_moment(edges,charges,powers)==Fraction(answer)
        assert coefficient(edges,charges,powers,root)==Fraction(answer)
        n=len(charges)
        relabeled=[(n-1-v,n-1-u) for u,v in edges]
        rotated=tuple((-y,x) for x,y in reversed(charges))
        assert capacity_moment(*reduced_input(relabeled,rotated,powers))==Fraction(answer)

    # Removing bridge 24 in the ten-edge graph leaves collinear loads on
    # each side. Their directions (1,1) and (1,0) give squared norms 2 and 1.
    # The four left edges therefore scale the scalar moment by 2^4.
    left=2**4*scalar_coefficient(((0,1),(0,1),(0,2),(1,2)),(-1,0,1),(0,)*4)
    right=scalar_coefficient(((0,1),(0,2),(0,2),(1,2),(1,2)),(0,-1,1),(0,)*5)
    assert (left,right)==(Fraction(1,168),Fraction(1,241920))
    assert -left*right==Fraction(CASES[-1][-1])
    print('Five original area charts and the independent scalar bridge product agree.')


if __name__=='__main__':
    check()
