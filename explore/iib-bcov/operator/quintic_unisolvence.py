"""Unisolvent transverse tests that vanish when any scalar momentum is zero."""
from itertools import combinations, permutations
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from genus_two_probes import det
from quartic_unisolvence import FIXTURES as QUARTIC, pivot_rows


def expand(p, extra):
    p = [*p[:-1],extra]
    return tuple([*p,tuple(-sum(z) for z in zip(*p))])


FIXTURES = (expand(QUARTIC[0],(2,-1,1,-1,1)),
            expand(QUARTIC[1],(-1,2,0,1,-2)),
            expand(QUARTIC[2],(1,0,-1,2,2)))


def soft_orbits():
    subsets={s:i for i,s in enumerate(combinations(range(9),5))}
    permutations4 = list(permutations(range(4)))
    masks={order:[sum(1<<order[i] for i in range(4) if h>>i&1) for h in range(16)]
           for order in permutations4}
    def fillings(h, n, rem):
        if h==16:
            if n==0 and not any(rem):
                yield ()
            return
        places=[i for i in range(4) if h>>i&1]
        upper=min(n,*(rem[i] for i in places))
        if h.bit_count()%2:
            upper=min(upper,1)
        for m in range(upper+1):
            rest=tuple(x-m*(i in places) for i,x in enumerate(rem))
            if min(rest)<0 or max(rest)>n-m:
                continue
            for tail in fillings(h+1,n-m,rest):
                yield (h,)*m+tail
    words=[]
    for h0 in (0,1,3,7,15):
        stabilizer=[images for images in masks.values() if images[h0]==h0]
        for tail in fillings(1,8,tuple(5-(h0>>i&1) for i in range(4))):
            if any(tuple(sorted(images[h] for h in tail))<tail for images in stabilizer):
                continue
            flags=(h0,*tail)
            words.append(tuple(subsets[tuple(i for i,h in enumerate(flags) if h>>j&1)]
                               for j in range(4)))
    return words


def eval_soft(exe,p,orbits):
    vals=[]
    for i,k in enumerate(p):
        other=[q for j,q in enumerate(p) if j!=i]
        inp='\n'.join(' '.join(map(str,r)) for r in [*other,(len(orbits),),*orbits])
        res=run([exe,'sum'],text=True,input=inp,capture_output=True,check=True)
        cols=[list(map(int,row.split())) for row in res.stdout.splitlines()]
        pivot=next(j for j,x in enumerate(k) if x)
        for j in range(5):
            if j!=pivot:
                vals.append([(k[pivot]*a[j]-k[j]*a[pivot])%1009 for a in cols])
    return vals


def check():
    chosen_shapes=(0,1,2,10,14,16,17,19,20,21,22,23,27,31,38,39,43,44,49,51,52,53,54,55,66,67,82,84,93,97,98,100,102,103,120,121,128,153,158,159,164,183,184)
    all_shapes=soft_orbits()
    orbits=[all_shapes[i] for i in chosen_shapes]
    rows=[]
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'quintic-certificate'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('quartic_probe_orbits.cpp')),
             '-o',str(exe)],check=True)
        for p in FIXTURES:
            rows.extend(eval_soft(exe,p,orbits))
            ids=pivot_rows(rows)
            print('Soft rank',len(ids),flush=True)
            if len(ids)==43:
                break
    assert len(ids)==43
    chosen=[rows[i] for i in ids]
    cols=pivot_rows(zip(*chosen))
    assert det([[r[j] for j in cols] for r in chosen])%1009==175
    print('The exact 43-row minor is 175 modulo 1009.')


if __name__=='__main__':
    check()
