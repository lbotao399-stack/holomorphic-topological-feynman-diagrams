"""Reconstruct each bare geometric kernel as a quartic in the three external determinants."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import prod
from pathlib import Path
import json
import sympy as sp

from reconstruct import graphs, here, normalizer, reconstruct
from reconstruct_cores import Reduction, symmetries


def momenta(y,z):
    return [(1,0),(0,1),(-z,y),(z-1,-y-1)]


GRID = [(y,z) for y in range(5) for z in range(5-y)]
POWERS = [(4-y-z,y,z) for y,z in GRID]
INVERSE = sp.Matrix([[y**j*z**k for _,j,k in POWERS] for y,z in GRID]).inv()


def coefficients(index,binary,canonical):
    red = Reduction(graphs[index])
    split = [len(c)-2 for c in red.chains]
    def sample(sources):
        total = defaultdict(int)
        for exponent,c in reconstruct(red,sources,binary).items():
            total[canonical(tuple(i+j for i,j in zip(exponent,split)))]+=c
        return {e:c for e,c in total.items() if c}
    grid = [sample(momenta(y,z)) for y,z in GRID]
    keys = sorted(set().union(*grid))
    expressions = [{} for _ in POWERS]
    for powers in keys:
        row = INVERSE*sp.Matrix([values.get(powers,0) for values in grid])
        for polynomial,value in zip(expressions,row):
            if value:
                assert value.q==1
                polynomial[powers]=int(value)
    # Reserved integral momenta have [01]=2 and both other determinants nonzero.
    extra=[(2,0),(0,1),(-1,-1),(-1,0)]
    actual = sample(extra)
    weights = [2**i*(-2)**j for i,j,k in POWERS]  # [12]=1
    for key in set(keys)|set(actual):
        assert sum(w*p.get(key,0) for p,w in zip(expressions,weights))==actual.get(key,0)
    return expressions


def build(binary,output):
    output.mkdir(exist_ok=True,parents=True)
    for start,stop in ((0,3),(3,11),(11,39),(39,50)):
        red = Reduction(graphs[start])
        group = symmetries(red)
        @cache
        def canonical(e):
            return min(tuple(e[j] for j in p) for p in group)
        polynomials = []
        support = set()
        for index in range(start,stop):
            expressions = coefficients(index,binary,canonical)
            polynomials.append(expressions)
            support.update(set().union(*expressions))
        support = sorted(support)
        # One check per orbit, before performing any cut. The least cotree
        # intersection with a subset is its loop number.
        masks=[sum(1<<i for i in s) for s,v in red.minors.items() if v]
        subsets=[([i for i in range(9) if S>>i&1],min((m&S).bit_count() for m in masks))
                 for S in range(1,511)]
        assert all(2*sum(e[i]+1 for i in subset)>11*b for e in support for subset,b in subsets)
        ids={e:j for j,e in enumerate(support)}
        rows=[[[[ids[e],c] for e,c in sorted(poly.items())] for poly in graph] for graph in polynomials]
        record=dict(start=start,edges=red.core,denominator=normalizer,powers=POWERS,
                    support=support,graphs=rows)
        (output/f'core{start:02}.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
        print(start,stop-1,'monomials',len(support),'nonzero polynomials',
              sum(bool(p) for g in rows for p in g),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    build(args.binary.resolve(),args.output)
