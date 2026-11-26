"""Row and column SL(2) integration by parts on the K3,3 Schwinger periods."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations,permutations
from math import factorial,prod

from flint import nmod_mat
from quartic_projective_ibp import AUTOMORPHISMS,K33
from quartic_planar_periods import forest_data


@cache
def convergence(q):
    cotrees=[tuple(map(int,e)) for e,c in forest_data(tuple(K33),6)[1].terms()]
    tests=[]
    for k in range(1,9):
        for ids in combinations(range(9),k):
            h=min(sum(e[i] for i in ids) for e in cotrees)
            if h==0:continue
            # With nonnegative weights, a subset of the same cycle rank
            # implies the inequality on any superset up to its next cycle.
            if any(set(old)<=set(ids) and limit==q*h for old,limit in tests):continue
            tests.append((ids,q*h))
    return tuple(tests)


def finite(p,q):
    return max(p)<q and all(sum(p[i] for i in ids)>h for ids,h in convergence(q))


@cache
def orientations():
    edges={e:i for i,e in enumerate(K33)}
    rows=(0,2,4);cols=(1,3,5)
    return tuple(tuple(tuple(edges[tuple(sorted((i,j)))] for j in b) for i in a)
                 for a,b in ((rows,cols),(cols,rows)) for a in permutations(a) if a[0]<a[1])


def canonical(p):
    return min(tuple(p[i] for i in g) for g in AUTOMORPHISMS)


def orbits(q):
    rows=set()
    def build(prefix,total):
        if len(prefix)==8:
            if 0<total<q:
                v=prefix+(total,)
                if finite(v,q):rows.add(canonical(v))
            return
        for k in range(max(1,total-(8-len(prefix))*(q-1)),min(q-1,total-(8-len(prefix)))+1):
            build(prefix+(k,),total-k)
    build((),4*q)
    return sorted(rows)


def equations(p):
    for first,second,middle in orientations():
        row=Counter({p:sum(p[i]-p[j] for i,j in zip(first,second))})
        for i,j,k in zip(first,second,middle):
            for a,sgn in ((i,1),(j,-1)):
                new=list(p);new[a]-=1;new[k]+=1
                row[tuple(new)]+=sgn*p[k]
        yield {v:c for v,c in row.items() if c}


def system(q,all_rows=None):
    full=orbits(q) if all_rows is None else all_rows
    lookup={p:i for i,p in enumerate(full)}
    contacts=set()
    equations_out=[]
    for p in full:
        for row in equations(p):
            normal=Counter()
            for v,c in row.items():normal[canonical(v)]+=c
            normal={v:c for v,c in normal.items() if c}
            if not normal or not all(finite(v,q) for v in normal):continue
            contacts.update(v for v in normal if not min(v))
            equations_out.append(normal)
    return full,sorted(contacts),equations_out


def verify(q=3):
    from pathlib import Path
    import json
    from quartic_projective_ibp import faces,reduction
    full,contacts,identities=system(q)
    index={v:i for i,v in enumerate(full)}
    rank=nmod_mat([[r.get(v,0)%1009 for v in full] for r in identities],1009).rank()
    print('dimension, positive orbits, wheel contacts, valid equations, interior rank',
          q,len(full),len(contacts),len(identities),rank,flush=True)
    if q!=3:return full,contacts,identities
    reference=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())
    values={}
    for v in reference['K33']:
        p=tuple(v['powers'])
        normalized=F(factorial(q-1),prod(factorial(a-1) for a in p))
        values[canonical(p)]=tuple(normalized*F(j) for j in v['1,zeta3,zeta5'])
    _,_,keep=reduction(faces(K33,5)[1])
    for v in reference['wheel']:
        short=[0]*8
        for i,j in enumerate(keep):short[j]=v['powers'][i]
        p=short[:5]+[0]+short[5:]
        normalized=F(factorial(q-1),prod(factorial(a-1) for a in p if a))
        values[canonical(p)]=tuple(normalized*F(j) for j in v['1,zeta3,zeta5'])
    for row in identities:
        assert all(sum(c*values[p][i] for p,c in row.items())==0 for i in range(3))
    return full,contacts,identities


if __name__=='__main__':
    for q in (3,4,5):verify(q)
