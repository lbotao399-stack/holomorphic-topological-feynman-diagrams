"""Calibrate the two nonbubble masters at d=11 and integrate individual K3,3 numerators."""
from collections import Counter
from fractions import Fraction
from functools import cache
import json
import networkx as nx
import sympy as sp

from bipartite_cut import cuts
from bipartite_masters import period
from odd_cut import HERE, beta_product_limit, replay


@cache
def orbit(exponent):
    return tuple(sorted({tuple(exponent[i] for i in permutation)
                         for permutation in replay.AUTOMORPHISMS}))


def canonical_record(row):
    terms = Counter()
    for powers, coefficient in row['terms']:
        exponent = tuple(i+j for i,j in zip(powers,row['h'],strict=True))
        terms[orbit(exponent)[0]] += coefficient
    return dict(edges=row['edges'], h=[0]*9, denominator=row['denominator'],
                orbit_representatives=True, terms=sorted((k,v) for k,v in terms.items() if v))


@cache
def proper_subgraphs():
    result = []
    for mask in range(1,511):
        edges = tuple(i for i in range(9) if mask>>i&1)
        graph = nx.Graph()
        graph.add_nodes_from(range(6))
        graph.add_edges_from(replay.K33[i] for i in edges)
        loops = len(edges)-6+nx.number_connected_components(graph)
        if loops:
            result.append((edges,loops))
    return result


@cache
def converges(exponent):
    return (len(exponent)==9 and sum(exponent)==13 and min(exponent)>=0
            and all(2*sum(exponent[i]+1 for i in edges)>11*loops
                    for edges,loops in proper_subgraphs()))


def reduce_row(directory, executable, tag, expression):
    result = replay.run_form(directory, executable, tag, expression)
    if result=='0':
        return [replay.FIELD.zero]*5
    return replay.read_masters(result,sp.Rational(7,2))


def calibration_rows(directory, executable, choices):
    record = json.loads((HERE/'reconstruct_cyclic/core39-49.json').read_text())[0]
    assert record['N_power']==0
    assert all(converges(tuple(e)) for e,_ in record['terms'])
    return [reduce_row(directory,executable,f'calibrate_{i}',e)
            for i,e in enumerate(cuts(record)) if i in choices]


def eliminator(rows):
    first,second = ([v-u for u,v in zip(rows[0],row)] for row in rows[1:])
    determinant = first[0]*second[1]-first[1]*second[0]
    assert determinant

    def eliminate(row,minimum_order=0):
        # Each calibration difference is holomorphic and zero at d=11.
        a = (row[1]*second[0]-row[0]*second[1])/determinant
        b = (row[0]*first[1]-row[1]*first[0])/determinant
        assert all(min(replay.laurent(x,minimum_order),default=minimum_order)>=minimum_order
                   for x in (a,b))
        vector = [v+a*f+b*s for v,f,s in zip(row,first,second)]
        assert not any(vector[:2])
        return vector[2:]
    return eliminate


def calibration(directory, executable):
    eliminate = eliminator(calibration_rows(directory,executable,range(3)))

    def integrate(row):
        value = sp.expand(beta_product_limit(eliminate(row)))
        assert period(row)==value
        coeff = [value.coeff(sp.pi,k) for k in (2,4)]
        assert all(q.is_Rational for q in coeff) and value==sum(q*sp.pi**k for q,k in zip(coeff,(2,4)))
        return tuple(map(Fraction,coeff))
    return integrate


def two_cuts(directory, executable, integrate, tag, record):
    assert all(converges(tuple(e)) for e,_ in record['terms'])
    expressions = list(cuts(record))
    assert expressions[0]!='()'
    first = integrate(reduce_row(directory,executable,tag+'_min',expressions[0]))
    second = next((e for e in expressions[1:] if e!=expressions[0]),None)
    if second is not None:
        assert first==integrate(reduce_row(directory,executable,tag+'_alt',second)),tag
    return first


def cyclic(directory, executable, integrate):
    table = {}
    for graph in range(39,50):
        rows = json.loads((HERE/f'reconstruct_cyclic/graph{graph:02}.json').read_text())
        coefficients = {str(k):['0','0'] for k in (0,2,4)}
        for row in rows:
            rank = row['N_power']
            record = canonical_record(row)
            if record['terms']:
                value = two_cuts(directory,executable,integrate,f'graph{graph}_n{rank}',record)
                coefficients[str(rank)] = list(map(str,value))
        table[str(graph)] = coefficients
        print(graph,coefficients,flush=True)
    expected = {0:(Fraction(93253402,12247200),Fraction(-9461505,12247200)),
                2:(Fraction(734213669,183708000),Fraction(-74453325,183708000)),4:(0,0)}
    for k,pair in expected.items():
        assert tuple(sum(Fraction(row[str(k)][j]) for row in table.values()) for j in (0,1))==pair
    return dict(pi_powers=[2,4],graphs=table)


def solve(directory,executable):
    result = cyclic(directory,executable,calibration(directory,executable))
    assert result==json.loads((HERE/'bipartite_cyclic.json').read_text())


if __name__=='__main__':
    replay.run(solve)
