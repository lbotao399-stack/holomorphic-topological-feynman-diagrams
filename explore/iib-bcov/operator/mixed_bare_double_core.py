"""All two-double four-loop contributions to the mixed fourfold Ward projection."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, prod
import json
from pathlib import Path

import numpy as np

from mixed_bare_double_family import CORE, CYCLES, G, chamber_values, fixture, numerator, tree_flow
from mixed_bare_ranks import census, chamber_orbits, hard_waves, integral_orbits
from quartic_seven_plane_double import Graph, fmpz_mpoly_ctx, forest_data, plan_generic, sp_monomial


def verify_zeros(geometric, nonzero):
    """Independently check all six integer cofactors on every remaining chamber."""
    graphs = {}
    count = 0
    for marks,color,_,_ in chamber_orbits(CORE, geometric):
        if (marks,color) in nonzero:
            continue
        if marks not in graphs:
            graph = G(12, fixture(CORE, marks))
            graph.C = [[row[e] for e in range(9) for _ in range(1+marks.count(e))] for row in CYCLES]
            graphs[marks] = graph
        graph = graphs[marks]
        _, waves = hard_waves(color)
        source = [(0,)*4]*6 + [(int(a),int(b+c),int(d),int(c)) for a,b,c,d in waves]
        flow = tree_flow(graph, source)
        for i,j in combinations(range(4), 2):
            assert not numerator(graph,source,flow,[(3,i),(3,j)])
        count += 1
    return count


def expanded_check():
    """Fold the entire numerator onto nine paths before integrating the core."""
    marks, color = (0,0,0,1,7,7), (0,0,2,1,1,2)
    edges = fixture(CORE,marks)
    graph = G(12,edges)
    graph.C = [[r[e] for e in range(9) for _ in range(1+marks.count(e))] for r in CYCLES]
    (u,t,z),waves = hard_waves(color)
    src = [(0,)*4]*6 + [(int(a),int(b+c),int(d),int(c)) for a,b,c,d in waves]
    flow = tree_flow(graph,src)
    minors = {(i,j):numerator(graph,src,flow,[(3,i),(3,j)]) for i,j in combinations(range(4),2)}
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    paths = []
    for e in range(9):
        start = sum(map(len,paths))
        paths.append(range(start,start+marks.count(e)+1))
    times = ctx.gens()
    path_time = [sum(times[j] for j in group) for group in paths]
    square = ctx.constant(0)
    for p,q in combinations(range(9),2):
        minor = sum((CYCLES[i][p]*CYCLES[j][q]-CYCLES[j][p]*CYCLES[i][q])*x
                    for (i,j),x in minors.items())
        square += path_time[p]*path_time[q]*minor**2
    _,_,res = forest_data(edges,12)
    def cross(a,b):
        return sum(res[6+i][6+j] for i,j in product(a,b))
    voltage = 4*res[6+u[0]][6+u[1]] + cross(t,z)-cross(u,z)-cross(t,u)
    # The additive residue splits the product into disjoint monomial sets.
    def slices(poly):
        by_residue = defaultdict(dict)
        for ex,c in poly.terms():
            by_residue[sum((i+1)**2*k for i,k in enumerate(ex))%37][ex] = int(c)
        return {r:ctx.from_dict(v) for r,v in by_residue.items()}
    left,right = slices(square),slices(voltage)
    scale = factorial(22)**9
    folded = defaultdict(int)
    term_count = 0
    for r in range(37):
        poly = sum(f*right[(r-i)%37] for i,f in left.items() if (r-i)%37 in right)
        term_count += len(poly)
        for ex,c in poly.terms():
            ex = tuple(map(int,ex))
            powers = tuple(len(path)-1+sum(ex[j] for j in path) for path in paths)
            factor,rem = divmod(scale*prod(map(factorial,ex)),prod(map(factorial,powers)))
            assert not rem
            folded[powers] += int(c)*factor
    assert term_count == 1218668 and len(folded) == 8773
    plan = plan_generic(Graph(6,CORE))
    value = F(128,2*16*scale)*sum(k*sp_monomial(ex,13,plan) for ex,k in folded.items() if k)
    assert value == F(9409,5589402610176000000)
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--expanded-check',action='store_true',help='Direct 15-segment integral on one chamber')
    parser.add_argument('--zeros',action='store_true',help='Check the six minors on the entire zero complement')
    parser.add_argument('--integrate',action='store_true',help='Reintegrate the selected shapes')
    parser.add_argument('shapes',nargs='*',type=int)
    args = parser.parse_args()
    data = json.loads(Path(__file__).with_suffix('.json').read_text())
    assignments = np.indices((9,)*6,dtype=np.int8).reshape(6,-1)
    geometric,_,physical = census(CORE,assignments)
    groups = defaultdict(list)
    for marks,color,weight,_ in integral_orbits(CORE,physical):
        groups[marks].append((color,weight))
    assert len(groups) == len(data['shapes']) == 32
    nonzero = set()
    grand = F()
    signs = Counter()
    for index,((marks,rows),saved) in enumerate(zip(groups.items(),data['shapes'],strict=True)):
        assert marks == tuple(saved['marks'])
        expected = [(tuple(map(int,c)),h,F(v)) for c,h,v in saved['chambers']]
        assert [(c,h) for c,h,_ in expected] == rows
        subtotal = sum(v for _,_,v in expected)
        assert subtotal == F(saved['subtotal'])
        grand += subtotal
        for c,_,v in expected:
            signs[(v>0)-(v<0)] += 1
            if v:
                nonzero.add((marks,c))
        if args.integrate and (not args.shapes or index in args.shapes):
            for result,reference in zip(chamber_values(marks,rows,all_pairs=True),expected,strict=True):
                assert result == reference
            print('Integrated',index,marks,subtotal,flush=True)
    assert len(nonzero) == 797 and signs[0] == 24
    assert grand == F(data['total'])
    if args.expanded_check:
        print('Independent 15-segment integration:',expanded_check())
    if args.zeros:
        assert verify_zeros(geometric,nonzero) == 1744
    print('Signed two-double core:',grand,'positive, negative, zero retained orbits:',
          signs[1],signs[-1],signs[0])


if __name__ == '__main__':
    main()
