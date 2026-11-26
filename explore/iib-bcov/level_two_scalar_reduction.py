"""Genus-two nine-wave worldsheet reductions on two saturated coordinate lines.

The four auxiliary primaries are forgotten. The two-psi vertex at a single
hard mark is its ancestor plus four lower one-psi vertices minus six primary
fusions. For two distinct hard first descendants, each forgotten point can
fuse into either mark. The double fusions use distinct auxiliary points and
distinct receiving marks. A rational tail containing both retained marks
survives stabilization, and gives no cotangent pullback correction.
"""
from fractions import Fraction as F
from itertools import combinations, permutations

from level_one_scalar_reduction import primary
from two_loop_moving_descendant import BASIS, descendant_vertices


def momenta(x,y,strengths):
    distinguished = (0,0,0,x,y)
    extras = tuple(tuple(a*t for t in u) for a,u in zip(strengths,(BASIS[3],BASIS[3],BASIS[4],BASIS[4]),strict=True))
    balance = tuple(-sum(col) for col in zip(*BASIS[:3],distinguished,*extras))
    return (*BASIS[:3],balance,distinguished),extras


def fuse(hard,extras,assignments):
    remaining = list(hard)
    for point,mark in assignments:
        remaining[mark] = tuple(a+b for a,b in zip(remaining[mark],extras[point]))
    removed = {j for j,_ in assignments}
    return (*remaining,*(p for j,p in enumerate(extras) if j not in removed))


def reduced_terms(x,y,split,strengths=(1,-1,1,-1)):
    """Return ancestor, sum of positive single fusions, and sum subtracted."""
    hard,extras = momenta(x,y,strengths)
    a,b,c,d = strengths
    fiber = 4*(a*b*c*d)**2
    if split:
        # The e^3 mark carries psi and the distinguished zero projected mark
        # is a dilaton on Mbar(2,4). The labeled dilaton factor is 2g-2+4.
        ancestor = fiber*F(6,9720)
        singles = sum(descendant_vertices(fuse(hard,extras,[(j,receiver)]))[other]
                      for receiver,other in ((2,4),(4,2)) for j in range(4))
        doubles = sum(primary(list(fuse(hard,extras,[(j,2),(k,4)])))
                      for j,k in permutations(range(4),2))
    else:
        # The distinguished mark with psi squared has zero projected wave.
        ancestor = fiber*F(59,77760)
        singles = sum(descendant_vertices(fuse(hard,extras,[(j,4)]))[4] for j in range(4))
        doubles = sum(primary(list(fuse(hard,extras,[(j,4),(k,4)])))
                      for j,k in combinations(range(4),2))
    return ancestor,singles,doubles


def vertex(x,y,split):
    # Each five-minor is affine in d=x e^4+y e^5: the two varying rows
    # contain d and -d, so the quadratic wedge vanishes. Both lower vertices
    # are quartic in these minors. Reflections and exchange of axes 4 and 5
    # give the four polynomial terms fixed by (0,0),(1,0),(2,0),(1,1).
    if split:
        coeffs = F(7,1296),F(49,19440),F(1,15552),F(1,3888)
    else:
        coeffs = F(5,864),F(23,19440),F(1,31104),F(1,7776)
    return sum(n*m for n,m in zip(coeffs,(1,x*x+y*y,x**4+y**4,x*x*y*y)))


def check():
    for split in (False,True):
        for x,y in ((0,0),(1,0),(2,0),(1,1),(3,-1),(F(-2,3),F(5,7))):
            a,b,c = reduced_terms(x,y,split)
            assert a+b-c==vertex(x,y,split)
    assert vertex(1,1,False)==F(649,77760) and vertex(1,1,True)==F(421,38880)
    assert reduced_terms(1,1,True)==(F(1,405),F(137,12960),F(43,19440))

    # One zero first-descendant input upstairs: forget it with the dilaton
    # equation, leaving the eight-wave moving first descendant at e^3.
    hard,extras = momenta(0,0,(1,-1,1,-1))
    assert vertex(0,0,True)==10*descendant_vertices((*hard[:4],*extras))[2]==F(7,1296)

    # One zero primary on the first auxiliary line: the string equation lowers
    # the single level-two mark or either of the two level-one marks.
    hard,extras = momenta(1,1,(0,-1,1,-1))
    lower = descendant_vertices((*hard,*extras[1:]))
    for split in (False,True):
        a,b,c = reduced_terms(1,1,split,(0,-1,1,-1))
        assert a==0 and a+b-c==lower[4]+(lower[2] if split else 0)

    print('Moving level-two scalar at e^4+e^5:',vertex(1,1,False))
    print('Moving level-one scalars at e^3 and e^4+e^5:',vertex(1,1,True))
    print('Zero-momentum level-two scalar:',vertex(0,0,False))
    print('Zero-momentum first descendant with hard level one:',vertex(0,0,True))


if __name__=='__main__':
    check()
