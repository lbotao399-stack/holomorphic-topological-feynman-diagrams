"""Fivefold genus-two vertex with one moving first descendant and seven primary waves."""
from fractions import Fraction
from itertools import combinations
from math import lcm

from bulk_eight import BULK, SHAPES, potentials
from genus_one_graph_ward import plus, separate_divergences
from longitudinal_theta import longitudinal
from marked_boundaries import quartics


def fused_symbols(momenta):
    """The seven-point polynomial P at the classical cubic cut on each unordered pair."""
    values = {}
    for i,j in combinations(range(8),2):
        fused = tuple(p for k,p in enumerate(momenta) if k not in (i,j))+(plus((momenta[i],momenta[j])),)
        values[i,j] = sum(c*y for c,y in zip(BULK,quartics(fused)))
    return values


def descendant_vertices(momenta):
    """All eight choices of the level-one wave, in the positive A_k convention."""
    if len(momenta)!=8 or any(len(p)!=5 for p in momenta) or any(plus(momenta)):
        raise ValueError('Eight conserved five-component waves are required.')
    denominator = lcm(*(Fraction(x).denominator for p in momenta for x in p))
    integer = tuple(tuple(int(denominator*Fraction(x)) for x in p) for p in momenta)
    fused = fused_symbols(integer)
    return tuple((longitudinal(integer[:i]+integer[i+1:])
                  - sum(value for group,value in fused.items() if i in group))/denominator**20
                 for i in range(8))


def graph_test(momenta, by_shape=False):
    """Integrate the original four theta shapes and use only the known primary trace in the tree."""
    if by_shape:
        shapes = [potentials(momenta,check=True,shapes=(shape,)) for shape in SHAPES]
        theta = [plus(group) for group in zip(*shapes)]
    else:
        theta = potentials(momenta,check=True)
    fused = fused_symbols(momenta)
    tree = [[Fraction(0)]*5 for _ in momenta]
    for (i,j),value in fused.items():
        charge = plus((momenta[i],momenta[j]))
        # The local action V_(2;0)=-P, with the other end of the tree grounded.
        for slot in (i,j):
            for k in range(5):
                tree[slot][k] -= value*charge[k]
    recovered = separate_divergences(momenta,[plus(pair) for pair in zip(theta,tree)])
    assert recovered==descendant_vertices(momenta)
    return recovered


BASIS = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))


def family(t):
    # Root e^1 and vector e^2 with polarization e_1. Wave 6 carries the descendant.
    p = (*BASIS,BASIS[1],tuple(t*x for x in BASIS[0]))
    return p+(tuple(-x for x in plus(p)),)


def check():
    generic = (*BASIS,(1,2,-1,3,2),(-2,1,3,1,-2))
    generic += (tuple(-x for x in plus(generic)),)
    graph_test(generic)
    p = family(1)
    result = graph_test(p,by_shape=True)
    assert all(any(row) for row in p) and result[6]==Fraction(19,12960)
    print('All four theta shapes plus the cut primary trace at finite momentum:',result[6])
    distinct = (*BASIS,plus(BASIS[:2]),plus(BASIS[2:4]))
    distinct += (tuple(-x for x in plus(distinct)),)
    assert all(any(p[i]*q[j]!=p[j]*q[i] for i,j in combinations(range(5),2))
               for p,q in combinations(distinct,2))
    assert graph_test(distinct)[6]==Fraction(191,17280)
    print('Eight pairwise nonparallel waves, descendant e^3+e^4:',Fraction(191,17280))
    for t in (-2,-1,0,1,2,3,4,Fraction(1,3)):
        value = descendant_vertices(family(t))[6]
        expected = Fraction(5*t**4+20*t**3+130*t**2+28*t+45,155520)
        assert value==expected
    # Five checks determine this polynomial: every five-minor is affine in t because
    # all varying momenta move in the same one-dimensional direction.
    soft = family(0)
    primary = tuple(q for i,q in enumerate(soft) if i!=6)
    primary_value = -sum(c*y for c,y in zip(BULK,quartics(primary)))
    eight = descendant_vertices(soft)
    assert eight[6]==9*primary_value==Fraction(1,3456)
    assert all(value==primary_value for i,value in enumerate(eight) if i!=6)
    assert descendant_vertices(tuple(tuple(0 if j==4 else x for j,x in enumerate(p)) for p in generic))==(0,)*8
    print('Moving amplitude, primary string limit, descendant dilaton limit and rank deficiency passed.')


if __name__=='__main__':
    check()
