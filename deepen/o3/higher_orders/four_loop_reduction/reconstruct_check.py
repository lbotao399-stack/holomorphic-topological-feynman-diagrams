"""Independent exact integrations of reconstructed numerators on the series-parallel core."""
from fractions import Fraction as F
from functools import cache
from pathlib import Path
import json
import sympy as sp

here=Path(__file__).resolve().parent
measures=here/'reconstruct_cyclic'

@cache
def simplex(params):
    assert min(params)>0
    return sp.prod(sp.gamma(s) for s in params)/sp.gamma(sum(params))


def monomial_period(d):
    a_outer=sum(d[i] for i in (0,3,4,5))
    b_outer=sum(d[i] for i in (2,6,7,8))
    half=sp.Rational(1,2)
    value=simplex((7-a_outer,7-b_outer,9*half-d[1]))
    for left,right,p,q in [(0,5,3,4),(2,6,7,8)]:
        value*=simplex((d[left]+1,d[right]+1,d[p]+d[q]-7*half))
        value*=simplex((9*half-d[p],9*half-d[q]))
    return sp.factor(value)


def series_parallel_coefficients():
    totals={}
    for idx in range(3):
        for entry in json.loads((measures/f'graph{idx:02}.json').read_text()):
            summed=sum(c*monomial_period(tuple(x+y for x,y in zip(power,entry['h'])))
                       for power,c in entry['terms'])
            totals[idx,entry['N_power']]=sp.factor(summed/entry['denominator'])
    assert totals=={(0,4):sp.pi**2/10500},totals
    return {4:totals[0,4]}

if __name__=='__main__':
    print(series_parallel_coefficients())
