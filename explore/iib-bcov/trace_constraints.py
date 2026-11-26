"""Exact evaluated three-loop trace boundary rows in the 18 soft-symbol basis."""
from fractions import Fraction as F
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import argparse
import json
import sys
from flint import fmpz_mat

from fourfold_bulk import fixture
from scaled_bulk import trace_rows
from planar_bulk import response
from bulk_three_hexagon import potentials


def measurements(exe):
    directory = Path(__file__).resolve().parent
    sys.path.insert(0,str(directory/'operator'))
    from three_loop_soft_jet import exact_boundary
    rows = []
    periods = []
    k = (-F(4059837269,70614415872000),F(1285913,27243216000))
    j = (-F(2078212919,635529742848000),F(131371,49037788800))
    for u in (((1,0),(1,0),(1,0),(0,1),(0,1),(0,1)),
              ((1,0),(0,1),(1,1),(1,-1),(1,0),(0,1))):
        p=[(*s,0,0,0) for s in u]+[(0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)]
        p.append(tuple(-sum(col) for col in zip(*p)))
        rows.append([a[0][1]-a[6][1] for a in exact_boundary(exe,p)])
        hx=potentials(p)
        assert hx[0][1]-hx[6][1]==0
        periods.append([response(u,x,y)[0][1] for x,y in zip(k,j)])
    p=fixture(3)
    arr=exact_boundary(exe,p)
    rows.extend([[a[i][2]-a[j][2] for a in arr] for i,j in ((3,7),(0,3))])
    periods.extend(([-F(23310136147,141228831744000),F(459733,3405402000)],
                    [-F(37217287,15406781644800),F(23449,14859936000)]))
    return rows,periods


def build(root):
    from fourfold_bulk import CORES
    from scaled_bulk import periods as core_sum
    with TemporaryDirectory() as td:
        exe=Path(td)/'nine'
        cpp=Path(__file__).resolve().parent/'operator'/'nine_sextics.cpp'
        run(['c++','-std=c++17','-O3',str(cpp),'-o',str(exe)],check=True)
        rows,rhs=measurements(exe)
        for strength in (1,2,3,4):
            new,hexagon=trace_rows(strength,exe)
            rows.extend(new)
            graphs=[]
            for name in ('tetrahedron','double-square'):
                file=root/f'scaled-{name}-{strength}.tsv' if root else None
                graphs.append(core_sum(CORES[name],strength,file)[1])
            for i,contact in enumerate(hexagon):
                rhs.append([contact+sum(g[2*i] for g in graphs),sum(g[2*i+1] for g in graphs)])
    matrix=fmpz_mat(rows)
    assert matrix.rank()==len(rhs)==11
    reduced,den,rank=matrix.rref()
    pivots=[next(j for j in range(18) if reduced[i,j]) for i in range(rank)]
    minor=fmpz_mat([[row[j] for j in pivots] for row in rows]).det()
    assert minor
    return dict(basis='three_loop_soft_jet.DOUBLE_SOFT',rows=rows,
                reduced_periods=[[str(c) for c in row] for row in rhs],
                pivot_columns=pivots,minor=str(minor))


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--checkpoints',type=Path)
    cli.add_argument('--save',action='store_true')
    args=cli.parse_args()
    result=build(args.checkpoints)
    file=Path(__file__).with_name('operator')/'three_loop_constraints.json'
    if args.save:
        file.write_text(json.dumps(result,indent=2)+'\n')
    else:
        assert result==json.loads(file.read_text())
    print('Eleven independent constraints; minor modulo 1009:',int(result['minor'])%1009)
    for row in result['reduced_periods']:print(*row)
