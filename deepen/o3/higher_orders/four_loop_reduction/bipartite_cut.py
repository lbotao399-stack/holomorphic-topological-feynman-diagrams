"""Exact cyclic color coefficients for the four-loop K3,3 core in five-dimensional CS."""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod
import json

from odd_cut import HERE, beta_product_limit, replay
import sympy as sp


def cuts(record):
    assert record['edges']==[list(e) for e in replay.K33]
    assert record['orbit_representatives'] and not any(record['h'])
    expressions = [defaultdict(Fraction) for _ in range(4)]
    for exponents, numerator in record['terms']:
        orbit = sorted({tuple(exponents[i] for i in p) for p in replay.AUTOMORPHISMS})
        choices = (orbit[0],orbit[-1],orbit[len(orbit)//2],orbit[len(orbit)//3])
        weight = Fraction(numerator*prod(map(factorial,exponents)),record['denominator'])
        for coefficients,exponent in zip(expressions,choices):
            coefficients[exponent[1:]] += weight
    for coefficients in expressions:
        terms = [f'+({w})*'+replay.propagators([k+1 for k in e])
                 for e,w in coefficients.items() if w]
        yield '('+'\n'.join(terms)+')'


def eliminate(rows):
    def minor(i,j):
        return rows[i][0]*rows[j][1]-rows[i][1]*rows[j][0]
    cofactors = [minor(1,2),minor(2,0),minor(0,1)]
    denominator = sum(cofactors,replay.FIELD.zero)
    weights = [c/denominator for c in cofactors]
    constants = [replay.laurent(weight,0) for weight in weights]
    assert all(min(c,default=0)>=0 for c in constants)
    assert sum(c.get(0,0) for c in constants)==1
    vector = [sum((w*row[j] for w,row in zip(weights,rows)),replay.FIELD.zero) for j in range(5)]
    assert not any(vector[:2])
    return beta_product_limit(vector[2:])


def solve(directory, executable):
    records = json.loads((HERE/'reconstruct_cyclic/core39-49.json').read_text())
    expected = {0:sp.pi**2*(93253402-9461505*sp.pi**2)/12247200,
                2:sp.pi**2*(734213669-74453325*sp.pi**2)/183708000}
    assert [row['N_power'] for row in records]==[0,2]
    for record in records:
        n = record['N_power']
        rows = []
        for j,expression in enumerate(cuts(record)):
            out = replay.run_form(directory,executable,f'bipartite{n}_{j}',expression)
            rows.append(replay.read_masters(out,sp.Rational(7,2)))
        for extra in (2,3):
            assert sp.expand(eliminate([rows[0],rows[1],rows[extra]])-expected[n])==0
        print(n,expected[n],flush=True)
    return expected


if __name__=='__main__':
    replay.run(solve)
