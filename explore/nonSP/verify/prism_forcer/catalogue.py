"""Exact zero-displacement periods of all thirteen C2 four-loop graphs."""
from fractions import Fraction as F
from functools import cache
import json
from math import prod
from pathlib import Path
import sympy as sp

from four_loop import K33, PRISM, WHEEL, checked_cut, current_terms, k4orbit
from wheel import orbit as wheel_orbit
from htk.catalog import tight_graphs
from evaluate import evaluate as prism
from all_sources import SEEDS, exact_coefficients
from check import GRAPH as bipartite, EDGE_INDEX, permutation_sign
import networkx as nx

Z3,Z5=sp.symbols('z3 z5')


@cache
def graph_catalogue():
    graphs=tight_graphs(2,4)
    assert len(graphs)==13 and graphs[11].edges==PRISM and graphs[12].edges==K33
    return [g.edges for g in graphs]


@cache
def wheel_table():
    record=json.loads(Path(__file__).with_name('wheel.json').read_text())
    assert record['edges']==list(map(list,WHEEL))
    return {tuple(p):list(map(F,c)) for p,c in record['cut_integrals']}


def evaluate(index,source):
    if type(index)!=int or not 0<=index<13:
        raise ValueError('Use the catalog graph number 0 through 12.')
    source=[tuple(map(F,row)) for row in source]
    if len(source)!=6 or any(len(row)!=2 for row in source) or any(map(sum,zip(*source))):
        raise ValueError('Expected six conserved two-component momenta.')
    if index==11:
        return prism(source)
    if index==12:
        brackets=[[a*d-b*c for c,d in source] for a,b in source]
        values=[0]*len(SEEDS)
        for g in nx.isomorphism.GraphMatcher(bipartite,bipartite).isomorphisms_iter():
            sign=permutation_sign([EDGE_INDEX[frozenset((g[i],g[j]))] for i,j in K33])
            for k,seed in enumerate(SEEDS):
                values[k]+=sign*prod(brackets[g[i]][g[j]] for i,j in seed)
        v=sp.expand(sum(n*c for n,c in zip(values,exact_coefficients())))
        return F(v.subs({Z3:0,Z5:0})),F(v.coeff(Z3)),F(v.coeff(Z5))
    name,q,weights=current_terms(graph_catalogue()[index],source,())
    assert q==6
    if name=='edge':
        return weights.get((0,),F(0)),F(0),F(0)
    result=[F(0),F(0),F(0)]
    for p,c in weights.items():
        if name=='wheel':
            row=wheel_table()[wheel_orbit(p)]
        else:
            assert name=='kite'
            row=checked_cut(k4orbit(p))
        for j,x in enumerate(row):
            result[j]+=c*x
    return tuple(result)


if __name__=='__main__':
    import sys
    if len(sys.argv)==3:
        print(json.dumps(list(map(str,evaluate(int(sys.argv[1]),json.loads(sys.argv[2],parse_float=F))))))
    else:
        record=json.loads(Path(__file__).with_name('catalogue.json').read_text())
        for i,row in enumerate(record['graphs']):
            assert graph_catalogue()[i]==list(map(tuple,row['edges']))
            assert evaluate(i,record['source'])==tuple(map(F,row['coefficients']))
        witnesses=[
            (8,[(1,0),(0,1),(0,0),(1,1),(0,0),(-2,-2)],(-F(239,75600),F(3,1400),0)),
            (10,[(1,0),(0,1),(0,0),(1,1),(-2,-2),(0,0)],(-F(3209,604800),F(13,4200),0)),
            (11,[(1,0),(0,1),(-1,0),(0,-1),(0,0),(0,0)],(F(169,151200),-F(3,2800),0)),
        ]
        for i,source,period in witnesses:
            assert evaluate(i,source)==period
        assert evaluate(12,[(1,0),(0,1),(1,1),(-2,-2),(0,0),(0,0)])[2]==F(9,10)
        print('All thirteen graph periods and the coefficient-field witnesses reproduced.')
