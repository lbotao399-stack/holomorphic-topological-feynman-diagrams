"""Eighteen single-component transverse measurements on two sparse ten-wave fixtures."""
from itertools import combinations
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from flint import nmod_mat
from three_loop_soft_jet import DOUBLE_SOFT, native


def fixtures():
    basis=[tuple(int(i==j) for j in range(5)) for i in range(5)]
    def scale(n,p):
        return tuple(n*x for x in p)
    first=[basis[0],scale(2,basis[0]),basis[1],scale(3,basis[1]),
           basis[2],scale(-1,basis[2]),basis[3],basis[4],(0,-2,-1,-1,1)]
    second=[basis[0],scale(2,basis[0]),basis[1],scale(3,basis[1]),*basis[2:],
            (1,0,1,-2,-2),(-2,0,1,1,-2)]
    return [p+[tuple(-sum(row) for row in zip(*p))] for p in (first,second)]


PROBES=(((0,1,4),(0,2,4),(0,3,4),(0,4,4),(0,5,4),(0,6,4),
         (0,1,2),(0,2,2),(0,3,2),(0,6,2),(0,7,2),
         (0,1,1),(0,4,1),(0,5,1),(0,6,1),(2,3,0),(2,4,0)),((0,1,1),))


def check(exe):
    rows=[]
    for p,probes in zip(fixtures(),PROBES):
        pairs=list(combinations(range(10),2))
        keys=[tuple(sorted([x for k,x in enumerate(p) if k not in (i,j)]
                           +[tuple(a+b for a,b in zip(p[i],p[j]))])) for i,j in pairs]
        unique=list(dict.fromkeys(keys))
        vals=dict(zip(unique,native(exe,unique,DOUBLE_SOFT,False)))
        for i,j,k in probes:
            assert p[i][k]==p[j][k]==0
            row=[0]*18
            for (a,b),key in zip(pairs,keys):
                charge=(int(a==i)-int(a==j))*p[b][k]+(int(b==i)-int(b==j))*p[a][k]
                row=[x+charge*y for x,y in zip(row,vals[key])]
            rows.append(row)
    assert int(nmod_mat(rows,1009).det())==963
    print('17 three-pair measurements and one two-pair measurement: determinant 963 modulo 1009.')


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'nine'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(exe)],check=True)
        check(exe)
