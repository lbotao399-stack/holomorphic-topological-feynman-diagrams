"""Fivefold scalar differential on a planar triangle and three transverse primary waves."""
from fractions import Fraction as F
from math import factorial, prod

from elliptic_descendants import vertex as elliptic
from genus_one_graph_ward import plus
from genus_two_graph_ward import vertex as genus_two
from scalar_ancestor_response import scalar_power
from two_loop_moving_descendant import BASIS
from betti_flow_polytope import weak
from mixed_sextic_vertex import sextic


TRIANGLE = ((1,0),(0,1),(-1,-1))


def graph_value(genus, levels, counts, directions=BASIS[2:]):
    """Unlabeled identical fields receive their factorials when evaluating the background."""
    waves = (*BASIS[:2], *(p for p,n in zip(directions,counts) for _ in range(n)))
    waves = (*waves,tuple(-x for x in plus(waves)))
    degrees = (*levels,*(0,)*(len(waves)-2))
    return (elliptic if genus==1 else genus_two)(waves,degrees)


def graph_polynomial(genus, levels):
    """Compute the homogeneous response from the small-genus fivefold graph banks."""
    from flint import fmpq, fmpq_mpoly_ctx
    axes = fmpq_mpoly_ctx.get(('a','b','c')).gens()
    total = axes[0]*0
    by_shape = {}
    for counts in weak(3*genus,3):
        shape = tuple(sorted(counts))
        if shape not in by_shape:
            by_shape[shape] = graph_value(genus,levels,shape)
        coefficient = by_shape[shape]/prod(map(factorial,counts))
        total += fmpq(str(coefficient))*prod(x**m for x,m in zip(axes,counts))
    return total


def check():
    from flint import fmpq, fmpq_mpoly_ctx
    a,b,c = fmpq_mpoly_ctx.get(('a','b','c')).gens()
    for genus, levels, expected in ((1,(1,0),F(1,24)), (2,(2,0),F(1,3456)), (2,(1,1),F(11,17280))):
        value = scalar_power(genus,2*genus,TRIANGLE,(*levels,0))
        assert value == expected
        poly = graph_polynomial(genus, levels)
        # Verify the entire plane a+b+c=0 symbolically, including vanishing amplitudes.
        assert poly.compose(a,b,-a-b) == fmpq(str(value))*(a*b*(-a-b))**genus
        print('Graph response in the action normalization:',genus,levels,poly,flush=True)
    difference = graph_polynomial(2,(1,1))-2*graph_polynomial(2,(2,0))
    assert difference == (a*b*c)**2/17280
    for levels, expected in (((3,0),F(7,1296000)), ((2,1),F(17,907200))):
        rational, zeta = sextic(3,TRIANGLE,(*levels,0))
        assert zeta == 0 and (-1)**3*(-2)**3*rational == expected
        print('Heat-normalized Q3 on amplitudes (1,1,-2):',levels,expected,flush=True)


def check_oblique():
    from flint import fmpz_mat
    directions = (*BASIS[2:], (0,0,1,1,1))
    amplitudes = ((1,-1,2,-2),(1,2,-4,1))
    for genus, levels, ground in ((1,(1,0),F(1,24)),(2,(1,1),F(11,17280))):
        responses = [F(),F()]
        for counts in weak(3*genus,4):
            value = graph_value(genus,levels,counts,directions)
            for i,weights in enumerate(amplitudes):
                responses[i] += value*prod(F(a**m,factorial(m)) for a,m in zip(weights,counts))
        for weights,measured in zip(amplitudes,responses):
            matrix = fmpz_mat([[sum(a*q[i]*q[j] for a,q in zip(weights,directions))
                                for j in (2,3,4)] for i in (2,3,4)])
            assert sum(weights)==0 and measured == int(matrix.det())**genus*ground
        print('Oblique transverse wave and singular Hessian:',genus,*responses,flush=True)


def check_hamiltonian():
    from flint import fmpq, fmpq_mpoly_ctx
    from mixed_sextic_vertex import sextic
    # These two planar test waves have determinant five; the unit-area bivector contracts
    # the scalar momentum into their Hamiltonian vector with the same determinant.
    pair=((2,1),(-1,2),(-1,-3))
    for levels,unit in (((3,0,0),F(7,10368000)),((2,1,0),F(17,7257600))):
        assert sextic(3,pair,levels)==(5**6*unit,0)
    a,b,c=fmpq_mpoly_ctx.get(('a','b','c')).gens()
    cyclic=((a,b-c),(b,c-a),(c,a-b))
    cubic=-fmpq(1,24)
    quintic=fmpq(1,3456)
    assert sum(x*y for x,y in cyclic)==0
    assert sum(cubic*(x**3*y+x*y**3) for x,y in cyclic)==0
    result=sum(quintic*(x**5*y+x*y**5)+cubic**2*x**3*y**3 for x,y in cyclic)
    assert result==a*b*c*(a-b)*(b-c)*(c-a)/432
    assert result(1,2,3)==fmpq(1,36)
    print('Projected two-loop Hamiltonian cyclic composition:',result(1,2,3),flush=True)


if __name__ == '__main__':
    check()
    check_oblique()
    check_hamiltonian()
