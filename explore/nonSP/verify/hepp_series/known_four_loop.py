"""Check 9-edge Hepp moments against the independently integrated K3,3 period."""
import sys
from fractions import Fraction as F
from pathlib import Path

from period import HeppSeries

sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'k33_positive'))
from check import K33, AUTOMORPHISMS, grid_input, numerator


def zeta_interval(s,n):
    partial=sum(F(1,k**s) for k in range(1,n+1))
    return partial+F(1,(s-1)*(n+1)**(s-1)),partial+F(1,(s-1)*n**(s-1))


def replay():
    _,records=grid_input()
    weights=numerator(records)
    seen=set()
    polynomial={}
    for row in zip(*weights.nonzero()):
        exponent=tuple(map(int,row))
        if exponent in seen:
            continue
        orbit={tuple(exponent[e] for e in p) for p in AUTOMORPHISMS}
        seen.update(orbit)
        polynomial[exponent]=F(int(weights[row])*len(orbit),6*len(AUTOMORPHISMS))
    assert len(polynomial)==175
    graph=HeppSeries(K33,6,polynomial)
    z3,z3_upper=zeta_interval(3,50)
    z5,z5_upper=zeta_interval(5,50)
    lower=(-539-1414*z3_upper+2160*z5)/14400
    upper=(-539-1414*z3+2160*z5_upper)/14400
    for k,term,partial,bound in graph.terms(2,permutations=AUTOMORPHISMS):
        partial,bound=F(str(partial)),F(str(bound))
        assert partial-bound<lower<upper<partial+bound
        print('K3,3',k,'partial',partial,'tail <=',bound,flush=True)


if __name__=='__main__':
    replay()
