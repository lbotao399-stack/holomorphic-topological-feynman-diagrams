"""Exact series and error bounds for the non-series-parallel cyclic I4 probe.

The core numerators already contain the orbit sums and graph symmetry factors.
"""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

from period import HeppSeries

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT/'deepen/o3/higher_orders/four_loop_reduction/reconstruct_cyclic'


def engine(record):
    assert record['orbit_representatives'] and len(record['edges'])==9
    terms={tuple(a+b for a,b in zip(e,record['h'])): F(c,record['denominator'])
           for e,c in record['terms']}
    g=HeppSeries(record['edges'],F(11,2),terms)
    assert len(g.cotree_permutations())==record['symmetries']
    return g


def series(core, order):
    for record in json.loads((SOURCE/(core+'.json')).read_text()):
        g=engine(record)
        values=list(g.terms(order,permutations='matroid'))
        row=dict(core=core,N_power=record['N_power'],monomials=len(g.numerator),
                 cotrees=len(g.cotrees),center=str(g.expansion('gegenbauer')[0]),
                 radius=str(g.expansion('gegenbauer')[1]),terms=[str(t) for k,t,p,b in values],
                 partial=str(values[-1][2]),tail=str(values[-1][3]))
        print(json.dumps(row,separators=(',',':')),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--order',type=int,default=2)
    p.add_argument('--cores',nargs='+',default=['core03-10','core11-38','core39-49'])
    args=p.parse_args()
    for core in args.cores:
        series(core,args.order)
