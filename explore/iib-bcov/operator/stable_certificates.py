"""Soft transverse certificates for the stable genus-two primary equation."""
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from quartic_unisolvence import pivot_rows
from stable_lifts import SIXTH_FIXTURES, SEVENTH_FIXTURES, EIGHTH_FIXTURE, final_fixture, cost
from stable_marked import exterior_shapes, soft_test


def select(exe,fixtures,target,chunk=24):
    n=len(fixtures[0])-1
    shapes=sorted(exterior_shapes(n,True),key=cost)
    chosen=[]
    cols=[]
    for start in range(0,len(shapes),chunk):
        subset=shapes[start:start+chunk]
        rows=[]
        for p in fixtures:
            rows.extend(soft_test(exe,subset,p))
        addcols=list(zip(*rows))
        ids=pivot_rows(cols+addcols)
        chosen += [start+i-len(cols) for i in ids if i>=len(cols)]
        cols=[(cols+addcols)[i] for i in ids]
        if len(chosen)==target:
            break
    rr=pivot_rows(zip(*cols))
    from genus_two_probes import det
    residue=det([[col[i] for col in cols] for i in rr])%1009
    print('New-target size',n,len(chosen),'shapes',chosen,'rows',rr,'residue',residue,flush=True)
    return chosen,rr,residue


CERTIFICATES = (
    (SIXTH_FIXTURES, (1,4,5,8,13,18,19,22,37,40,43,50,65,87,89,106,107,118,120,128,130,
                      182,194,209,212,243,255,274,300,307,309,326), tuple(range(31))+(44,),41),
    (SEVENTH_FIXTURES[:1],(0,2,4,10,18,21,25,29,45,61,86,90,108,118,166,197),tuple(range(16)),203),
    ((EIGHTH_FIXTURE,), (0,4,11,22), (0,1,2,3),122),
    ((final_fixture(),), (0,), (0,),865))


def check(exe):
    from genus_two_probes import det
    for fixtures,ids,row_ids,residue in CERTIFICATES:
        shapes=sorted(exterior_shapes(len(fixtures[0])-1,True),key=cost)
        chosen=[shapes[i] for i in ids]
        rows=[]
        for p in fixtures:
            rows.extend(soft_test(exe,chosen,p))
        assert det([rows[i] for i in row_ids])%1009==residue
        print('Exact minor of size',len(ids),'is',residue,'modulo 1009.',flush=True)


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'stable-orbits'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('stable_orbits.cpp')),
             '-o',str(exe)],check=True)
        check(exe)
