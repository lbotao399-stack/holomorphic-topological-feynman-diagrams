"""First displacement derivative on the reflection-symmetric six-source prism."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod
import subprocess
import os

from flint import fmpz_mpoly_ctx
from periods import CUTS, HERE, PERMS, edge_order, form, read_vectors, tagged
from shift import period_evaluator
from sources import GRAPH, PRISM, density_exact, grid_input, numerator

# The odd edge permutation induced by (23)(45) fixes all six sources.
SOURCE = [(-2,0),(0,-2),(1,0),(1,0),(0,1),(0,1)]


def density_and_currents(source,binary):
    trees,records = grid_input(source,PRISM)
    raw = numerator(records,len(trees)**6*density_exact(GRAPH,[1]*9,source[:-1]),
                    symmetrize=False,executable=binary)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(9)))
    density = ctx.from_dict({tuple(map(int,row)):int(raw[row]) for row in zip(*raw.nonzero())})
    currents = [[ctx.from_dict({ex:record[4+2*e+i] for ex,record in zip(trees,records)
                               if record[4+2*e+i]}) for i in range(2)] for e in range(9)]
    graph = ctx.from_dict(dict.fromkeys(trees,1))
    for v in range(6):
        for i in (0,1):
            assert sum((int(v==a)-int(v==b))*row[i] for (a,b),row in zip(PRISM,currents))==source[v][i]*graph
    for cycle in GRAPH.C:
        for i in (0,1):
            assert sum(s*t*c[i] for s,t,c in zip(cycle,ctx.gens(),currents))==0
    return density,currents


def aggregate(polynomial,order):
    cache_rep = {}
    terms = defaultdict(int)
    for ex,c in polynomial.terms():
        ex = tuple(map(int,ex))
        if ex not in cache_rep:
            orbit = {tuple(ex[e] for e in p) for p in PERMS}
            cache_rep.update(dict.fromkeys(orbit,min(orbit)))
        terms[cache_rep[ex]] += int(c)
    assert all(sum(ex)==15+4*order for ex in terms)
    return {ex:v for ex,v in terms.items() if v}


def expressions(terms,order,cut):
    expr = []
    for ex,c in terms.items():
        powers = max(tuple(ex[p[e]] for e in edge_order(cut)) for p in PERMS)
        norm = Fraction(c*prod(map(factorial,ex)),factorial(5+order))
        den = '*'.join(f'1/(p{j}.p{j})^{k+1}' for j,k in enumerate(powers) if j)
        expr.append(f'+({norm})*'+den)
    return '('+''.join(expr)+')'


def moment_jobs(folder,executable):
    binary=folder/'grid'
    subprocess.run([os.getenv('CXX','c++'),'-std=c++17','-O3',str(HERE.parent/'k33_positive/grid.cpp'),
                    '-o',str(binary)],check=True)
    density,currents = density_and_currents(SOURCE,binary)
    assert density and not aggregate(density,0)
    # Edge 23 is reversed by the odd reflection, so its two moments can survive.
    adjacent = currents[PRISM.index((2,3))]
    rung = currents[PRISM.index((2,5))]
    rows = [aggregate(density*v,1) for v in (*adjacent,rung[0])]
    # Verify the full current vector modulo time permutations before integration.
    # Coordinates below multiply the two components on 23 and the first on 25.
    a = ((1,0,0),(0,1,0))
    b = a[::-1]
    def add(x,y):
        return tuple(tuple(i+j for i,j in zip(row,col)) for row,col in zip(x,y))
    def neg(x):
        return tuple(tuple(-i for i in row) for row in x)
    r = ((0,0,1),(0,0,1))
    zero = ((0,0,0),(0,0,0))
    predicted = [zero,add(a,r),neg(add(a,r)),add(b,r),neg(add(b,r)),a,r,neg(r),b]
    for edge in range(9):
        for j,coef in enumerate(predicted[edge]):
            poly = density*(currents[edge][j] - coef[0]*adjacent[0]
                            - coef[1]*adjacent[1] - coef[2]*rung[0])
            assert not aggregate(poly,1),(edge,j)
    print('moment orbit counts',list(map(len,rows)),flush=True)
    pair = []
    for cut in CUTS:
        text = form(folder,executable,f'moment{cut}',
                    tagged(enumerate(expressions(row,1,cut) for row in rows)),cut)
        pair.append(read_vectors(text))
    evaluate = period_evaluator(folder, executable, 4, shifts=5)
    result = [[evaluate(row[i]) for i in range(3)] for row in pair]
    assert result[0]==result[1]
    print('edge 23 and the first component on 25',result[0],flush=True)
    import sympy as sp
    z3 = sp.Symbol('z3')
    a,b,c = result[0]
    assert sp.expand(28*a+20*b+35*c) == (sp.Rational(588751,907200)-sp.Rational(5123,9450)*z3)
    return result[0]


if __name__=='__main__':
    from replay import run
    run(moment_jobs)
