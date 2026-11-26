"""Stable primary genus-two scalar-mark basis and remaining soft equations."""
from collections import Counter
from fractions import Fraction as F
from math import factorial,prod
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

from higher_genus_two import add,dot,terms,transverse,general_fused_marked,QUARTIC_MARK
from quintic_marks import QUINTIC_ORBITS,QUINTIC_MARK
from quintic_unisolvence import FIXTURES as QUINTIC_FIXTURES,expand
from quartic_unisolvence import pivot_rows
from stable_marked import exterior_shapes,flags_of,orbit_weight,native,fused


def cost(shape):
    return factorial(len(shape))//prod(factorial(m) for m in Counter(shape).values())


def next_basis(n,exe):
    rng=Random(3104)
    pts=[]
    for _ in range(7):
        rows=[tuple(rng.randrange(-2,3) for j in range(5)) for i in range(n-1)]
        pts.append(rows+[tuple(-x for x in add(*rows))])
    shape=sorted(exterior_shapes(n,False),key=cost)
    all_cols=native(exe,shape,pts,modular=True)
    ids=pivot_rows(all_cols)
    print('Next scalar dimension',n,len(ids),[shape[i] for i in ids],flush=True)
    return [shape[i] for i in ids]


def base_mark(p,exe):
    lower=general_fused_marked(p)[2:6]
    shapes=[flags_of(w) for w in QUINTIC_ORBITS]
    extra=fused(exe,shapes,p)
    coeffs=[c/orbit_weight(s) for s,c in zip(shapes,QUINTIC_MARK)]
    return [*lower,*extra],[*QUARTIC_MARK,*coeffs]


SIXTH = (((3,3,3,3,3,12,12,12,12,12), F(5,1536)),
         ((1,2,3,3,3,12,12,12,12,15), -F(5,15552)),
         ((1,2,3,3,4,8,12,12,15,15), -F(149,37324800)))
SIXTH_FIXTURES = (expand(QUINTIC_FIXTURES[0],(-1,1,2,0,-2)),
                  expand(QUINTIC_FIXTURES[1],(2,0,-1,-1,1)))


SEVENTH = (((1,2,3,3,3,3,12,12,12,12,12), F(1,1280)),
           ((1,2,3,3,3,4,8,12,12,12,15), -F(133,2985984)))
SEVENTH_FIXTURES = (expand(SIXTH_FIXTURES[0],SIXTH_FIXTURES[0][0]),
                    expand(SIXTH_FIXTURES[1],SIXTH_FIXTURES[1][2]))
EIGHTH_FIXTURE = expand(SEVENTH_FIXTURES[0],SEVENTH_FIXTURES[0][1])

EIGHTH = (((1,2,3,3,3,3,4,8,12,12,12,12), F(1,23040)),)

def marked_equation(p,exe,known,new):
    from flint import fmpq_mat
    lower,weights=base_mark(p,exe)
    other=fused(exe,[s for s,w in known],p) if known else []
    arr=fused(exe,new,p) if new else []
    raw=terms(p)
    rows=transverse(p,[*arr,*lower,*other,*raw])
    r=len(new); q=len(weights)+len(known)
    coeffs=[*weights,*(c for s,c in known)]
    data=[[*z[:r],dot(coeffs,z[r:r+q])+sum(z[r+q:])] for z in rows]
    if not new:
        assert all(row==[0] for row in data)
        print('All',len(rows),'transverse equations vanish at',len(p),'external legs.',flush=True)
        return
    m=fmpq_mat([[str(x) for x in row] for row in data])
    reduced,rank=m.rref()
    print('Rank and fitted residual',rank,reduced.tolist()[:r+1],flush=True)



def final_fixture():
    basis=[tuple(int(i==j) for i in range(5)) for j in range(5)]
    rows=basis*2+[basis[0],basis[1],(2,-1,-2,1,0)]
    return tuple([*rows,tuple(-x for x in add(*rows))])


def final_row(exe, known, final_shape):
    from higher_genus_two import overlap_symbols
    p=final_fixture()
    # Leg 0 has momentum e^1; the first contraction uses the transverse vector e_2.
    raw=sum(a[0][1] for a in terms(p))
    q=Counter()
    for j in range(1,len(p)):
        charge=p[j][1]
        if charge:
            rows=[p[h] for h in range(len(p)) if h not in (0,j)]+[add(p[0],p[j])]
            q[tuple(sorted(rows))]+=charge
    q={k:c for k,c in q.items() if c}
    lower=sum(w*dot(QUARTIC_MARK,overlap_symbols(k)[2:6]) for k,w in q.items())
    first=[flags_of(w) for w in QUINTIC_ORBITS]
    weights=[c/orbit_weight(s) for c,s in zip(QUINTIC_MARK,first)]
    shapes=[*first,*(s for s,c in known),final_shape]
    vals=[dot(x,q.values()) for x in native(exe,shapes,list(q))]
    base=raw+lower+dot(vals[:-1],[*weights,*(c for s,c in known)])
    print('Final row before the twelve-slot term',base,'its multiplier',vals[-1],flush=True)
    return base,vals[-1]

if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'stable-orbits'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('stable_orbits.cpp')),
             '-o',str(exe)],check=True)
        for p in SIXTH_FIXTURES:
            marked_equation(p,exe,SIXTH,[])
        for p in SEVENTH_FIXTURES:
            marked_equation(p,exe,SIXTH+SEVENTH,[])
        marked_equation(EIGHTH_FIXTURE,exe,SIXTH+SEVENTH+EIGHTH,[])
        a,b=final_row(exe,SIXTH+SEVENTH,EIGHTH[0][0])
        assert (a,b)==(F(419,40),-241344)
        assert a+EIGHTH[0][1]*b==0
