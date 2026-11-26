"""Independent discrete Pixton tests of the scalar Newton moments on cycles of parallel bundles."""
from fractions import Fraction
from itertools import accumulate
from math import comb, factorial, prod

from betti_flow_polytope import charge_cuts, normal_coefficient


def convolution(a, b):
    """Unnormalized circular convolution of the integer weight lists."""
    r = len(a)
    return [sum(v*b[(j-i) % r] for i,v in enumerate(a)) for j in range(r)]


def weighting(r, charges, bundles, scale):
    """A single bundle is a bridge; otherwise bundles form an oriented cycle."""
    tables = []
    for powers in bundles:
        table = None
        for m in powers:
            edge = [(u*(r-u))**m for u in range(r)]
            table = edge if table is None else convolution(table,edge)
        tables.append(table)
    if len(bundles)==1:
        numerator = tables[0][scale*charges[0] % r]
        rank = len(bundles[0])-1
    else:
        offsets = tuple(accumulate(charges))
        numerator = sum(prod(t[(u+scale*p) % r] for t,p in zip(tables,offsets)) for u in range(r))
        rank = 1+sum(len(b)-1 for b in bundles)
    return Fraction(numerator,r**rank)


def discrete_top(charges, bundles):
    exponents = tuple(m for b in bundles for m in b)
    degree = 2*sum(exponents)
    scaled = []
    for t in range(degree+1):
        # Use r beyond the ordering walls. Interpolate r^rank times the residue
        # average; cutting at the edge residues gives unimodular integral cells.
        rank = sum(len(b)-1 for b in bundles)+int(len(bundles)>1)
        length = degree+rank+1
        start = 1+t*sum(p for p in charges if p>0)
        moduli = range(start,start+length)
        values = [weighting(r,charges,bundles,t)*r**rank for r in moduli]
        # Newton forward interpolation about start; coefficient r^rank.
        basis = [Fraction(1)]
        coefficient = Fraction(0)
        for k in range(length):
            if k>=rank: coefficient += values[0]*basis[rank]
            values = [b-a for a,b in zip(values,values[1:])]
            basis = [((basis[i-1] if i else 0) -
                      (start+k)*(basis[i] if i<len(basis) else 0))/(k+1)
                     for i in range(min(k+1,rank)+1)]
        scaled.append(coefficient)
    top = sum((-1)**(degree-t)*comb(degree,t)*v for t,v in enumerate(scaled))/factorial(degree)
    return top*(-1)**(sum(exponents)-len(exponents))/(2**sum(exponents)*prod(map(factorial,exponents)))


def check():
    examples = (
        ((2,-2),((1,1,1),)),
        ((2,-2),((1,2,1),)),
        ((1,2,-3),((1,),(1,),(1,))),
        ((1,2,-3),((2,),(1,),(1,))),
        ((1,2,-3),((1,1),(1,),(1,))),
        ((3,-1,2,-4),((1,),(1,),(1,),(1,))),
        ((1,1,-2),((1,1),(1,1),(1,))),
    )
    for charges,bundles in examples:
        if len(bundles)==1: edges = ((0,1),)*len(bundles[0])
        else:
            edges = tuple((i,(i+1)%len(bundles)) for i,b in enumerate(bundles) for _ in b)
        powers = tuple(m-1 for b in bundles for m in b)
        value = discrete_top(charges,bundles)
        assert value==normal_coefficient(charge_cuts(edges,charges),powers),(charges,bundles,value)
        print(charges,bundles,value,flush=True)


if __name__=='__main__':
    check()
