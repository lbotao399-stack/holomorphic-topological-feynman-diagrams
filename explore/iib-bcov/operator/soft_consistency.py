"""Eleven-primary Chevalley--Eilenberg consistency tensor for a soft bulk potential."""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

import numpy as np
from flint import nmod_mat
from genus_two_probes import det
from soft_target_basis import SUBSETS, basis_data, fixtures


def eleven(n=0):
    p=fixtures(n+1)[n][:-1]
    p.append((2,-1,3,1,-2) if n==0 else (-1,2,2,-3,1))
    p.append(tuple(-sum(t) for t in zip(*p)))
    return tuple(p)


def fused_arrays(exe,p):
    pairs=list(combinations(range(11),2));ts,ws,_=basis_data()
    data=[(len(pairs),len(ts),len(ws)),*ts,*ws]
    for a,b in pairs:
        q=[v for i,v in enumerate(p) if i not in (a,b)]
        q.append(tuple(x+y for x,y in zip(p[a],p[b])))
        data.extend(q)
        data.append([det([q[i] for i in s])%1009 for s in SUBSETS])
    inp='\n'.join(' '.join(map(str,x)) for x in data)
    r=run([exe],input=inp,text=True,capture_output=True,check=True)
    return np.fromstring(r.stdout,dtype=np.int16,sep=' ').reshape(len(pairs),len(ws),10,5)


def wedge_sign(a,b):
    sign=1
    while a:
        v=a&-a;a-=v
        if (b&(v-1)).bit_count()%2:sign=-sign
    return sign


@lru_cache(maxsize=12000)
def wedge(a,b):
    z={}
    for i,v in a:
        for j,w in b:
            if i&j==0:z[i|j]=z.get(i|j,0)+wedge_sign(i,j)*v*w
    return tuple((h,v) for h,v in sorted(z.items()) if v)


@lru_cache(maxsize=8000)
def interior(q,a):
    z={}
    for mask,v in a:
        for j,x in enumerate(q):
            b=1<<j
            if x and (b&mask):
                sign=(-1)**((mask&(b-1)).bit_count())
                z[mask-b]=z.get(mask-b,0)+sign*x*v
    return tuple((m,x) for m,x in sorted(z.items()) if x)


@lru_cache(maxsize=1000)
def polarized_weights(fields):
    one=((0,1),)
    before=[one];after=[one]*(len(fields)+1)
    for f in fields:before.append(wedge(before[-1],f))
    for i in range(len(fields)-1,-1,-1):after[i]=wedge(fields[i],after[i+1])
    out=[];deg=0
    basis=[tuple(int(s==j) for s in range(5)) for j in range(5)]
    for pos,f in enumerate(fields):
        if f[0][0]:
            for k,c in enumerate(basis):
                z=wedge(before[pos],wedge(interior(c,f),after[pos+1]))
                if z:out.append((pos,k,(-1)**deg*dict(z).get(31,0)))
        deg+=f[0][0].bit_count()
    return tuple((a,k,v) for a,k,v in out if v)


def vector(v):
    return tuple((1<<a,x) for a,x in enumerate(v) if x)


def transverse_vector(p,rng):
    a=next(i for i,x in enumerate(p) if x)
    while True:
        z=[rng.randrange(-2,3) if j!=a else 0 for j in range(5)]
        if any(z):break
    return tuple(-sum(v*x for v,x in zip(p,z)) if j==a else p[a]*z[j] for j in range(5))


def cochain_row(p,anchor,legs,pol):
    f=[((0,1),) for _ in p]
    f[anchor]=interior(p[anchor],((31,1),))
    for i,u in zip(legs,pol):f[i]=vector(u)
    ranks=[v[0][0].bit_count() for v in f]; out=[]
    for pair,(a,b) in enumerate(combinations(range(11),2)):
        if not(ranks[a] or ranks[b]):continue
        q=tuple(x+y for x,y in zip(p[a],p[b]))
        merged=interior(q,wedge(f[a],f[b]))
        if not merged:continue
        remaining=[t for t in range(11) if t not in (a,b)]
        # Move the pair to the end; the odd divergence crosses the preceding slots.
        exponent=sum(ranks[t] for t in remaining)
        exponent+=sum(ranks[t]*ranks[s] for s in (a,b) for t in remaining if s<t)
        sign=(-1)**exponent
        kept=[t for t in remaining if ranks[t]]
        fields=tuple([f[t] for t in kept]+[merged])
        offsets=[remaining.index(t) for t in kept]+[9]
        for j,k,c in polarized_weights(fields):
            out.append((50*pair+5*offsets[j]+k,c*sign))
    return out


def consistency(arrays,p,seed,count):
    rng=Random(seed);cols=arrays.transpose(0,2,3,1).reshape(55*50,-1).astype(np.int64)
    out=[]
    for _ in range(count):
        anchor,*legs=rng.sample(range(11),4)
        pol=[transverse_vector(p[i],rng) for i in legs]
        row=cochain_row(p,anchor,legs,pol)
        if row:
            ids,vals=zip(*row)
            out.append((np.asarray([int(v)%1009 for v in vals],dtype=np.int64)@cols[list(ids),:])%1009)
    return np.array(out)


def check_signs(p):
    """Test a symmetric scalar trace boundary and the gauge on each fused fiber."""
    arrays=[]
    for a,b in combinations(range(11),2):
        q=[v for i,v in enumerate(p) if i not in (a,b)]
        q.append(tuple(x+y for x,y in zip(p[a],p[b])))
        trace=np.zeros((10,5),dtype=np.int64)
        for i,j in combinations(range(10),2):
            merged=tuple(x+y for x,y in zip(q[i],q[j]))
            # A symmetric scalar polynomial on the nine fused momenta.
            values=[v for h,v in enumerate(q) if h not in (i,j)]+[merged]
            scalar=sum((v[0]+2*v[1]-v[2])**6 for v in values)%1009
            trace[i]+=scalar*np.asarray(q[j])
            trace[j]+=scalar*np.asarray(q[i])
        common=np.tile(q[0],(10,1))
        longitudinal=np.asarray(q)*np.arange(1,11)[:,None]
        arrays.append((trace%1009,common%1009,longitudinal%1009))
    assert any(np.any(a[0]) for a in arrays)
    zero=consistency(np.asarray(arrays),p,9163,64)
    assert len(zero)>=60 and not np.any(zero)


def check(exe, second_fixture=False):
    p=eleven()
    check_signs(p)
    arrays=fused_arrays(exe,p)
    rows=np.concatenate([consistency(arrays,p,5432+i,200) for i in range(3)])
    assert rows.shape==(598,518)
    minor=nmod_mat(rows[:500,:500].tolist(),1009)
    assert int(minor.det())==27
    print('First 500 rows and columns: determinant 27 modulo 1009.')
    assert nmod_mat(rows.tolist(),1009).rank()==500
    print('598 transverse rows, 518 soft potentials. Consistency rank 500.')
    if second_fixture:
        p=eleven(1)
        rows=consistency(fused_arrays(exe,p),p,19875,720)
        assert nmod_mat(rows.tolist(),1009).rank()==500
        print('Independent eleven-wave fixture: consistency rank 500.')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--second-fixture',action='store_true',
                        help='repeat the rank calculation at independent momenta')
    args=parser.parse_args()
    with TemporaryDirectory() as tmp:
        exe=Path(tmp)/'soft-target-orbits'
        run(['c++','-O3','-std=c++17','-DACCELERATE_NEW_LAPACK','-framework','Accelerate',
             str(Path(__file__).with_name('soft_target_orbits.cpp')),'-o',str(exe)],check=True)
        check(exe,args.second_fixture)
