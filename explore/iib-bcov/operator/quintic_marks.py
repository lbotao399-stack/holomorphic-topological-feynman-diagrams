"""Scalar counterterms with nine differentiated labels for the integrated quintic."""
from itertools import combinations
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from fractions import Fraction as F

from higher_genus_two import add, det


SUBSETS = list(combinations(range(9),5))


QUINTIC_ORBITS = ((17,26,85,104),(16,32,103,119),(50,102,106,119))


def fused_scalar_orbits(p, executable, words):
    pairs = list(combinations(range(len(p)),2))
    rows = []
    for i,j in pairs:
        q = [p[h] for h in range(len(p)) if h not in (i,j)]+[add(p[i],p[j])]
        rows.append([det([q[h] for h in t]) for t in SUBSETS])
    text = '\n'.join(' '.join(map(str,row)) for row in [(len(words),len(rows)),*words,*rows])
    response = run([executable],text=True,input=text,capture_output=True,check=True)
    vals = [list(map(int,line.split())) for line in response.stdout.splitlines()]
    out = [[[0]*5 for h in p] for _ in words]
    for (i,j),v in zip(pairs,vals):
        for h,c in enumerate(v):
            for k in range(5):
                out[h][i][k] += c*p[j][k]
                out[h][j][k] += c*p[i][k]
    return out


QUINTIC_MARK = (F(17,3110400),F(97,933120),-F(1,622080))


if __name__=='__main__':
    from higher_genus_two import terms,transverse,QUARTIC_MARK,dot,general_fused_marked
    from quintic_unisolvence import FIXTURES
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'scalar-mark-orbits'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('scalar_mark_orbits.cpp')),
             '-o',str(exe)],check=True)
        for p in FIXTURES:
            basis=fused_scalar_orbits(p,exe,QUINTIC_ORBITS)
            raw=terms(p)
            lower=general_fused_marked(p)[2:]
            rows=transverse(p,[*basis,*lower,*raw])
            assert all(dot(QUINTIC_MARK,r[:3])+dot(QUARTIC_MARK,r[3:7])+sum(r[7:])==0 for r in rows)
            print('Forty exact transverse quintic equations passed.',flush=True)
