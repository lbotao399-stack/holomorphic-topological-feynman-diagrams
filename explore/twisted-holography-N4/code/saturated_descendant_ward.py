"""Highest collinear multiplicity at fixed descendant degrees in the flat fivefold Ward equation."""
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import sys

IIB = Path(__file__).resolve().parents[2]/'iib-bcov'
sys.path.insert(0, str(IIB))
from elementary_scalar import elementary_three_loop
from genus_two_graph_ward import vertex as two_loop_vertex
from two_loop_single_descendant import primary as two_loop_primary

BASIS = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))


def multiplier(genus, strengths, levels):
    degree = sum(levels)
    assert len(strengths)==genus+degree
    return Fraction(factorial(degree),prod(factorial(k) for k in levels))*sum(
        prod(x*x for x in selected) for selected in combinations(strengths,genus))


@cache
def tree_recovery(genus,strengths,levels):
    """Evaluate all remaining genus-zero tree cuts at each raised hard marking."""
    assert len(strengths)==genus+sum(levels)
    if not any(levels):
        return prod(x*x for x in strengths)
    destination = next(i for i,k in enumerate(levels) if k)
    lower = tuple(k-int(i==destination) for i,k in enumerate(levels))
    # The hard charges modulo the line are independent of the chosen subgraph.
    # Retain separately the coefficient of p_i and the common axial potential.
    rows = []
    for i,k in enumerate(lower):
        rest = tuple(0 if j==i else h for j,h in enumerate(lower))
        hard = axial = 0
        for cut in combinations(range(len(strengths)),k+1):
            free = tuple(x for j,x in enumerate(strengths) if j not in cut)
            weight = tree_recovery(genus,free,rest)
            hard += weight
            axial += weight*sum(strengths[j] for j in cut)
        rows.append((hard,axial))
    assert len({r[1] for r in rows})==1
    return rows[destination][0]


def check_two_loops():
    for k in ((2,0,0,0,0),(1,1,0,0,0),(2,1,0,0,0),(1,1,1,0,0),(3,0,0,0,0)):
        for strengths in ((1,)*(2+sum(k)), (2,-1,3,-2,4)[:2+sum(k)]):
            hard = (*BASIS[1:],(-sum(strengths),-1,-1,-1,-1))
            line = tuple((a,0,0,0,0) for a in strengths)
            value = two_loop_vertex((*line,*hard), (0,)*len(line)+k)
            base = two_loop_primary((BASIS[0],BASIS[0],*BASIS[1:],(-2,-1,-1,-1,-1)))
            assert value == multiplier(2,strengths,k)*base == tree_recovery(2,strengths,k)*base
        print('Genus-two Ward control',k, value)


def check_three_loops():
    # Original five-fusion fixture: the third hard mark acquires one u in the primary cut.
    hard = [(1,0,-4,-4,-1),(-1,0,0,2,0),(0,1,1,0,0),(0,-1,0,1,0),BASIS[4]]
    scalars = [*hard,BASIS[2],BASIS[2]]
    base = tuple(-x for x in elementary_three_loop(scalars,BASIS[3],BASIS[2]))
    assert base==(Fraction(11180467393,94152554496000),-Fraction(220697,2270268000))
    for levels in ((0,0,2,0,0,0),(1,0,1,0,0,0),(1,1,1,0,0,0),(2,2,1,1,0,0)):
        for strengths in ((1,)*(3+sum(levels)), (2,-1,3,-2,4,1,2,-1,3)[:3+sum(levels)]):
            assert tree_recovery(3,strengths,levels)==multiplier(3,strengths,levels)
        scale = multiplier(3,(1,)*(3+sum(levels)),levels)
        print('Genus-three consecutive Ward trees',levels,*(scale*x for x in base))


if __name__=='__main__':
    check_two_loops()
    check_three_loops()
