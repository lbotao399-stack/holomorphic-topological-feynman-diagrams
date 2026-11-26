"""Reconstruct all thirteen four-loop bare source polynomials in one bracket basis."""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import json
import os
import subprocess
import sys

from flint import fmpz_mat
from sampling import PAIRS,SEEDS,evaluation,sampling
from cubature import HERE,write_input
from catalogue import evaluate,graph_catalogue


def compile_grid(folder):
    output=folder/'cubature'
    prefix=os.environ.get('GMP_PREFIX','/opt/homebrew/opt/gmp')
    subprocess.run([os.environ.get('CXX','c++'),'-std=c++17','-O3','-I'+prefix+'/include',
                    '-L'+prefix+'/lib',str(HERE/'cubature.cpp'),'-lgmpxx','-lgmp','-o',str(output)],check=True)
    return output


def compute(index,folder,binary):
    probes,matrix=sampling()
    if index in (11,12):
        vals=[evaluate(index,p) for p in probes]
    else:
        inp=folder/f'graph{index}.in'
        write_input(index,inp)
        pts=folder/'probes'
        pts.write_text(str(len(probes))+'\n'+'\n'.join(' '.join(str(x) for v in p[:5] for x in v)
                                                     for p in probes)+'\n')
        path=folder/f'graph{index}.out'
        subprocess.run([binary,inp,pts,path],check=True)
        vals=[tuple(map(F,line.split())) for line in path.read_text().splitlines()]
        assert len(vals)==len(probes)
        for j in (-2,-1):assert vals[j]==evaluate(index,probes[j])
    denominators=[lcm(*(row[j].denominator for row in vals)) for j in range(3)]
    ints=[[int(v*d) for v,d in zip(row,denominators)] for row in vals]
    base=fmpz_mat(matrix[:490].tolist())
    coeff=base.solve(fmpz_mat(ints[:490]))
    assert fmpz_mat(matrix.tolist())*coeff==fmpz_mat(ints)
    rows={str(i):[str(coeff[i,j]/denominators[j]) for j in range(3)] for i in range(490)
          if any(coeff[i,j] for j in range(3))}
    print(index,len(rows),'nonzero standard monomials',flush=True)
    return rows


if __name__=='__main__':
    import argparse
    from tempfile import TemporaryDirectory
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save',type=Path)
    args=parser.parse_args()
    with TemporaryDirectory() as temp:
        folder=Path(temp);binary=compile_grid(folder)
        rows=[compute(i,folder,binary) for i in range(13)]
    result={'pairs':PAIRS,'seeds':SEEDS,'edges':graph_catalogue(),'coefficients':rows}
    data=(json.dumps(result,separators=(',',':'))+'\n').encode()
    if args.save:args.save.write_bytes(data)
    else:assert data==(HERE/'polynomials.json').read_bytes()
