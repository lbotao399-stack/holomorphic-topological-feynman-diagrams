"""Exact denominator and numerator bounds for reconstructing sparse bulk projections."""
from fractions import Fraction as F
from math import ceil, comb, factorial, lcm, prod
from pathlib import Path
import argparse

from core_cubature import supports
from fourfold_bulk import CORES,symmetries


def core_bounds(core,file):
    entries = {tuple(map(int,k)):tuple(map(F,(r,s))) for *k,r,s in
               (row.split() for row in file.read_text().splitlines())}
    assert set(entries)=={orb[0] for orb in supports(core,12)}
    normalized = [tuple(v/prod(factorial(i) for i in ex) for v in values)
                  for ex,values in entries.items()]
    den = [lcm(*(r[k].denominator for r in normalized))*2*len(symmetries(core))
           for k in (0,1)]
    maximum = [max(abs(v[k]) for v in entries.values()) for k in (0,1)]
    return den,maximum


def bounds(p,tests,core,table):
    den,maximum = table
    heights = [sum(max(0,x) for x in axis) for axis in zip(*p)]
    # Parseval on the unit torus, then Hadamard on the integral 24-column
    # current determinant. Scaling each source component rescales three rows.
    square = comb(23,8)*24**16*16**16*8**8*prod(v**6 for v in heights)
    orders = factorial(10)*comb(15,5)
    out = []
    for left,right,axis in tests:
        v = [int(k==left)-int(k==right) for k in range(10)]
        weight = F(sum(abs(v[i]*p[j][axis]+v[j]*p[i][axis])
                       for i in range(10) for j in range(i)),2)
        ceiling = F(square*comb(16,4)*orders,len(symmetries(core)))*weight
        out.append([(d,ceil(d*ceiling*m)) for d,m in zip(den,maximum)])
    return out


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('tables',type=Path)
    args=cli.parse_args()
    from trace_reconstruction import waves,TESTS
    for name in ('tetrahedron','double-square'):
        core=CORES[name]
        table=core_bounds(core,args.tables/f'core-periods-{name}.tsv')
        sizes = [item for case,p in waves().items() for item in bounds(p,TESTS[case],core,table)]
        print(name,'denominators',*table[0],flush=True)
        print('maximum numerator bits:',*(max(n[k][1].bit_length() for n in sizes) for k in (0,1)))
