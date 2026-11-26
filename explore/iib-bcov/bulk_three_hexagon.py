"""Genus-two bulk trace in a one-loop hexagon on ten conserved primary waves."""
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import lcm
from flint import fmpz_mat,fmpz_mpoly_ctx

from bulk_eight import BULK
from marked_theta import green_moments
from marked_boundaries import quartics


def replaced(group,slot,row):
    return int(fmpz_mat([*group[:slot],row,*group[slot+1:]]).det())


@cache
def green(k):
    return tuple(green_moments(k,i,j)[1] for i in range(5) for j in range(i,5))


def marked_vertex(group,single):
    """Five scalar momenta at the trace, five singleton momenta around it."""
    volume = int(fmpz_mat(single).det())
    if not volume or fmpz_mat(group).rank()<4:
        return [(F(0),)*5]*5
    ctx = fmpz_mpoly_ctx.get(('s0','s1','s2','s3','s4'))
    s = ctx.gens()
    c = [sum(replaced(group,h,p)*t for p,t in zip(single,s)) for h in range(5)]
    c += [-int(fmpz_mat(group).det()),0]
    mean_numerator = sum(c)
    a = [7*x-mean_numerator for x in c]
    den = lcm(*(z.denominator for z in BULK))
    first,second = (int(den*z) for z in BULK)
    poly = first*sum(x*x for x in a)**2+second*sum(x**4 for x in a)
    moments = [F(0)]*15
    for exp,coefficient in poly.terms():
        for j,value in enumerate(green(tuple(map(int,exp)))):
            moments[j] += int(coefficient)*value
    # The five singleton assignments fill the entire cube; reflection has order two.
    # The polarized genus-two bulk vertex is -P times the trace.
    factor = -F(volume**2,2*den)
    rows = [[F(0)]*5 for _ in single]
    for (i,j),v in zip(combinations_with_replacement(range(5),2),moments):
        for h in range(5):
            rows[i][h] += factor*v*single[j][h]
            if i!=j:
                rows[j][h] += factor*v*single[i][h]
    return rows


def potentials(p):
    if len(p)!=10 or any(len(q)!=5 for q in p) or any(map(sum,zip(*p))):
        raise ValueError('Expected ten conserved five-component momenta.')
    total = [[F(0)]*5 for _ in p]
    for group in combinations(range(10),5):
        active = [i for i in range(10) if i not in group]
        result = marked_vertex([p[i] for i in group],[p[i] for i in active])
        for i,row in zip(active,result):
            for h in range(5):
                total[i][h] += row[h]
    return total


def cubature():
    import sys
    from pathlib import Path
    sys.path.insert(0,str(Path(__file__).resolve().parent/'operator'))
    from marked_hexagon import compositions
    from math import factorial
    # The rule exact through degree seven on the five-simplex, grouped over
    # orders, directly tests a Green diagonal with a five-slot polynomial.
    single = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    group = [tuple(-v for v in x) for x in single]
    ans = F(0)
    for j in range(4):
        scale = 12-2*j
        weight = F((-1)**j*scale**7,64*factorial(j)*factorial(12-j))
        for parts in compositions(3-j,6):
            lengths = [2*k+1 for k in parts]
            ordered = tuple(sum(lengths[:i]) for i in range(1,6))
            # Permuting axes preserves this group, hence P is symmetric in these tails.
            active = [tuple(scale*x for x in q) for q in group]
            pol = quartics(active+[ordered,tuple(scale-x for x in ordered)])
            pval = sum(x*y for x,y in zip(BULK,pol))/scale**20
            response = sum(x*(scale-x) for x in ordered)
            ans -= weight*12*pval*F(response,scale**2)
    assert ans == marked_vertex(group,single)[0][0]
    return ans


if __name__=='__main__':
    diagonal = cubature()
    assert diagonal == F(421,447897600)
    print('Independent hexagon diagonal:',diagonal)
    basis = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    p = basis+[(1,1,0,1,0),(0,1,-1,1,1),(1,0,1,2,-1),(1,-1,0,0,1)]
    p.append(tuple(-sum(col) for col in zip(*p)))
    rows = potentials(p)
    pairing = rows[0][2]-rows[1][2]
    assert pairing == F(44224151,1175731200)
    print('Ten-wave potential difference:',pairing)
