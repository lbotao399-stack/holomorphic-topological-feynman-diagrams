"""Exact coefficients of 1, zeta(3), zeta(5) for six conserved rational sources."""
from fractions import Fraction
from functools import cache
from itertools import permutations
import json
from math import prod
from pathlib import Path


@cache
def table():
    data = json.loads(Path(__file__).with_name('kernel.json').read_text())
    edges = list(map(tuple,data['edges']))
    automorphisms = []
    for g in permutations(range(6)):
        mapped = [tuple(sorted((g[u],g[v]))) for u,v in edges]
        if set(mapped)==set(edges):
            p = [edges.index(e) for e in mapped]
            sign = (-1)**sum(x>y for i,x in enumerate(p) for y in p[i+1:])
            automorphisms.append((g,sign))
    assert len(automorphisms)==12
    return automorphisms,[(seed,list(map(Fraction,row))) for seed,row in data['rows']]


def evaluate(source):
    source = [tuple(map(Fraction,p)) for p in source]
    if len(source)!=6 or any(len(p)!=2 for p in source) or any(map(sum,zip(*source))):
        raise ValueError('Expected six conserved two-component momenta.')
    automorphisms,rows = table()
    brackets = [[a*d-b*c for c,d in source] for a,b in source]
    result = [Fraction(0) for _ in range(3)]
    for seed,coefficients in rows:
        weight = sum(sign*prod(brackets[g[i]][g[j]] for i,j in seed) for g,sign in automorphisms)
        for k,c in enumerate(coefficients):
            result[k] += weight*c
    return tuple(result)


if __name__=='__main__':
    import sys
    source = json.loads(sys.argv[1],parse_float=Fraction)
    print(json.dumps(list(map(str,evaluate(source)))))
