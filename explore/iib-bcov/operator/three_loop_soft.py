"""Constraint from a constant scalar on the candidate three-loop bulk trace."""
from itertools import combinations
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from genus_two_probes import det
from quartic_unisolvence import FIXTURES,pivot_rows,multiplicities


WORDS=((29,29,29,29,29,29),(4,4,4,4,34,34),(40,44,44,44,52,52),(3,3,10,18,31,53),
       (2,12,18,29,31,50),(25,25,53,53,53,53),(20,20,20,20,23,26),(7,18,18,35,49,49),
       (11,21,27,31,33,33),(3,17,17,37,45,48),(27,27,27,41,53,53),(3,8,8,15,15,34),
       (10,24,44,44,47,49),(14,21,25,27,27,45),(5,5,40,40,44,44),(8,8,28,28,43,49))


def columns(exe,p,orbits):
    subsets=list(combinations(range(8),5))
    pairs=list(combinations(range(9),2))
    data=[(len(orbits),len(pairs)),*orbits]
    for i,j in pairs:
        q=[v for h,v in enumerate(p) if h not in (i,j)]
        q.append(tuple(x+y for x,y in zip(p[i],p[j])))
        data.append([det([q[h] for h in s]) for s in subsets])
    t='\n'.join(' '.join(map(str,row)) for row in data)
    result=run([exe],input=t,text=True,capture_output=True,check=True)
    all_fused=[list(map(int,line.split())) for line in result.stdout.splitlines()]
    out=[]
    for values in zip(*all_fused):
        A=[[0]*5 for _ in p]
        for (i,j),v in zip(pairs,values):
            for k in range(5):
                A[i][k]+=v*p[j][k]
                A[j][k]+=v*p[i][k]
        # Anchor the four-vector at the last leg. The other two marked choices carry vectors.
        vals=[]
        for i,j in combinations(range(8),2):
            for a in range(5):
                if p[i][a]:continue
                for b in range(5):
                    if p[j][b]:continue
                    vals.append(((A[i][a]-A[8][a])*p[8][b]-(A[j][b]-A[8][b])*p[8][a])%1009)
        out.append(vals)
    return list(zip(*out))


if __name__=='__main__':
    assert multiplicities(8,(6,)*5)[0]==16 and multiplicities(9,(6,)*5)[0]==50
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'sextic-traces'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('sextic_trace_orbits.cpp')),
             '-o',str(exe)],check=True)
        words_=WORDS
        rows=[]
        for p in FIXTURES:
            rows+=columns(exe,p,words_)
            sel=pivot_rows(rows)
            if len(sel)==16:break
    mat=[rows[i] for i in sel]
    labels=pivot_rows(zip(*mat))
    assert len(sel)==16
    assert det([[r[i] for i in labels] for r in mat])%1009==243
    print('The restricted trace boundary has rank 16; its determinant is 243 modulo 1009.')
