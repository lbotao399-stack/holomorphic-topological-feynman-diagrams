"""All external momenta for bare geometric graphs 0–10, prior to color and automorphisms."""
import argparse
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import json
import sympy as sp

from doubled_kite import finite_value, kite, shifted_bubble
from reconstruct_check import monomial_period


def doubled(exponents):
    a,b,c,d,e = [exponents[i]+1 for i in (0,2,3,6,4)]
    pair = kite(a,b,c,d,e)
    assert pair==kite(c,d,a,b,e)
    factor = prod(map(factorial,exponents))*shifted_bubble(exponents[7]+1,exponents[8])
    return finite_value(factor*pair[0],factor*pair[1])


def values(source,output):
    record=json.loads(source.read_text())
    table=[]
    for exponent in record['support']:
        value=monomial_period(tuple(exponent)) if record['start']==0 else doubled(exponent)
        value=sp.expand(value)
        h,k=value.coeff(sp.pi,2),value.coeff(sp.pi,4)
        assert h.is_Rational and k.is_Rational and value==h*sp.pi**2+k*sp.pi**4
        table.append((F(h),F(k)))
    graphs=[]
    for columns in record['graphs']:
        row=[]
        for terms in columns:
            row.append([str(sum((n*table[i][j] for i,n in terms),F())/record['denominator'])
                        for j in range(2)])
        graphs.append(row)
    result=dict(start=record['start'],powers=record['powers'],graphs=graphs)
    output.write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(output.name,len(table),'periods',len(graphs),'graphs',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    values(args.source,args.output)
