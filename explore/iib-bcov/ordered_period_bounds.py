"""Coefficientwise height bounds for ten labelled three-loop path insertions."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import ceil, comb, factorial, gcd, prod

import numpy as np
from flint import fmpz_mat

from core_cubature import supports
from fourfold_bulk import CORES, symmetries
from period_bounds import core_bounds

# Cycle columns in core-edge order, as in partial_cotrees_mod.cpp.
CYCLES = {
    'tetrahedron': ((1,1,0),(-1,0,1),(0,-1,-1),(1,0,0),(0,1,0),(0,0,1)),
    'double-square': ((-1,1,1),(1,0,0),(0,-1,-1),(0,1,1),(0,1,0),(0,0,1)),
}


def determinant_masses(name):
    cycles = CYCLES[name]
    group = symmetries(CORES[name])
    lengths = [sum(x*x for x in column) for column in cycles]
    minors = []
    for cols in ((0,1,2),(0,1),(0,2),(1,2)):
        terms = []
        for paths in combinations(range(6),len(cols)):
            value = abs(int(fmpz_mat([[cycles[e][i] for i in cols] for e in paths]).det()))
            if value:
                terms.append((paths,value))
        minors.append(terms)

    masses = Counter()
    for positions in combinations_with_replacement(range(6),10):
        counts = tuple(positions.count(e) for e in range(6))
        rep = min(tuple(counts[e] for e in permutation) for permutation in group)
        edge_terms = [[(paths,det*prod(counts[e]+1 for e in paths)) for paths,det in terms]
                      for terms in minors]
        candidates = []
        for scale in (1,2,3,4,5,8):
            # Scale each source row by the reciprocal of sqrt(scale) times
            # its half-flow bound. Fourteen column norms remain after jumps.
            norms = [scale*x+5 for x in lengths]
            squares = [sum(w for _,w in terms)*sum(w*prod(norms[e] for e in paths)
                                                   for paths,w in terms) for terms in edge_terms]
            value = prod(x*y for x,y in zip(lengths,norms))*squares[0]**2*min(squares[1:])
            candidates.append(F(value,scale**9))
        masses[rep] += factorial(10)*prod(x**n for x,n in zip(lengths,counts))*min(candidates)
    return masses


def path_heights(name,table_file):
    core = CORES[name]
    denominators,_ = core_bounds(core,table_file)
    orbits = supports(core,12)
    entries = {tuple(map(int,r[:6])):tuple(map(F,r[6:]))
               for r in (line.split() for line in table_file.read_text().splitlines())}
    moments = []
    for orbit in orbits:
        factorials = prod(factorial(e) for e in orbit[0])
        ns = [abs(d*x/factorials) for d,x in zip(denominators,entries[orbit[0]])]
        assert all(x.denominator==1 for x in ns)
        moments.extend([tuple(x.numerator for x in ns)]*len(orbit))
    powers = np.asarray([e for orbit in orbits for e in orbit],dtype=np.int8)
    local = [factorial(5)**(r//5)*factorial(r%5) for r in range(11)]
    total = [F(0),F(0)]
    for counts,mass in determinant_masses(name).items():
        difference = powers-np.array(counts,dtype=np.int8)
        valid = np.all((difference>=0)&(difference<=5*(np.array(counts)+1)),axis=1)
        largest = [0,0]
        for i,row in zip(np.flatnonzero(valid),difference[valid]):
            factor = prod(local[k] for k in row)
            for j,value in enumerate(moments[i]):
                largest[j] = max(largest[j],factor*value)
        for k in range(2):
            total[k] += mass*largest[k]
    return denominators,total


def ordered_bounds(name,table_file,source,tests,jumps=None):
    """Use the hard source for flows and voltages, and each jump once per minor."""
    if len(source)!=10 or any(sum(col) for col in zip(*source)):
        raise ValueError('Expected ten conserved five-dimensional waves.')
    if jumps is None:
        jumps = source
    halves = [F(sum(max(v,0) for v in col),2) for col in zip(*source)]
    if len(halves)!=5 or not all(halves):
        raise ValueError('All five source components must be present.')
    scale = prod(h**6 for h in halves)*prod(sum((x/h)**2 for x,h in zip(row,halves)) for row in jumps)
    content = prod(gcd(*row) for row in jumps)**2
    denominators,heights = path_heights(name,table_file)
    result = []
    for i,j,k in tests:
        v = [int(s==i)-int(s==j) for s in range(10)]
        weight = F(sum(abs(v[a]*source[b][k]+v[b]*source[a][k]) for a in range(10) for b in range(a)),2)
        factor = scale*weight*comb(16,4)/len(symmetries(CORES[name]))
        row = []
        for d,height in zip(denominators,heights):
            divisor = gcd(d,content)
            row.append((d//divisor,ceil(factor*height/divisor)))
        result.append(row)
    return result
