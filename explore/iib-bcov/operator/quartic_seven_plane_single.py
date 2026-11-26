"""Exact seven-planar quadratic-lift periods of the single-double four-loop core.

The two parallel paths reduce to one effective line of a tetrahedron. Ordinary
segment subdivisions integrate on simplices. The two propagator cuts of each
convergent tetrahedron moment return its rational and zeta(3) coefficients.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
import gzip
import hashlib
import json
import sys

from flint import fmpq_mpoly_ctx, fmpz_mpoly_ctx
from quartic_seven_plane_double import (
    G, Q, SCALE, fixture, forest_data, local_moment, polynomial_det, polynomials,
)
from square_three import checked_cut, orbit

CORE=((0,1),(0,1),(0,2),(1,3),(2,4),(2,5),(3,4),(3,5),(4,5))
CYCLES=((1,-1,0,0,0,0,0,0,0),(0,0,0,0,1,-1,0,0,1),(0,0,0,0,0,0,1,-1,1),(1,0,-1,1,-1,0,1,0,0))


@cache
def graph_data(marks):
    edges=fixture(CORE,marks)
    paths=[]
    for i in range(9):
        start=sum(map(len,paths))
        paths.append(tuple(range(start,start+marks.count(i)+1)))
    graph=G(11,edges)
    graph.C=[[row[i] for i,path in enumerate(paths) for _ in path] for row in CYCLES]
    ctx=fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
    times=ctx.gens()
    gram=[[sum(x*y*t for x,y,t in zip(row,col,times)) for col in graph.C] for row in graph.C]
    adj=[[(-1)**(i+j)*polynomial_det([row[:i]+row[i+1:] for k,row in enumerate(gram) if k!=j],ctx)
          for j in range(4)] for i in range(4)]
    _,u,dist=forest_data(edges,11)
    u=ctx.from_dict({k:int(v) for k,v in u.terms()})
    assert u==polynomial_det(gram,ctx)
    return graph,ctx,u,adj,dist,paths,(len(paths[0]),len(paths[1]))


def reduce_polynomial(data,integ):
    paths,split=data[5:]
    weights=defaultdict(F)
    bubble=sum(paths[:4],())
    cache_reduced={}
    for ex,coef in integ.terms():
        m=tuple(map(int,ex))
        eL=tuple(m[i] for i in bubble)
        powerL=sum(eL)+len(eL)-Q
        local=F(local_moment(split,eL),SCALE)
        rest=tuple(sum(m[i]+1 for i in path) for path in paths[4:])
        if rest not in cache_reduced:
            signature=(powerL,*rest)
            assert sum(signature)==3*Q
            cache_reduced[rest]=orbit(signature)
        k=cache_reduced[rest]
        # The five ordinary series-chain factorial denominators cancel those of the reduced K4 moment.
        local*=F(factorial(powerL-1)*prod(factorial(m[i]) for path in paths[4:] for i in path),factorial(Q-1))
        weights[k]+=int(coef)*local
    return {k:v for k,v in weights.items() if v}


def numerator_poly(data,ds,v):
    adj=data[3]
    return v*sum(ds[i]*adj[i][j]*ds[j]*(1+(i!=j)) for i,j in combinations_with_replacement(range(4),2))


def reduce_full(marks,order=tuple(range(5))):
    data=graph_data(marks)
    return reduce_polynomial(data,numerator_poly(data,*polynomials(data,order)))


@cache
def tetrahedron_bank():
    raw=gzip.decompress(Path(__file__).with_name('quartic_seven_plane_tetra.weights.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest()=='c96fa16bda07cf42925e381fab8d1e2dd93cde7c5d2d8b1574d30c473c875776'
    data={tuple(p):(F(a),F(b)) for p,a,b in map(json.loads,raw.splitlines())}
    assert len(data)==3872
    return data


def verify_bank():
    from kite_native import coefficients
    seen=set()
    for first in product(range(1,Q),repeat=5):
        last=3*Q-sum(first)
        if not 0<last<Q:
            continue
        p=(*first,last)
        triangles=((0,1,3),(0,2,4),(1,2,5),(3,4,5))
        if all(sum(p[i] for i in triangle)>Q for triangle in triangles):
            seen.add(orbit(p))
    assert seen==tetrahedron_bank().keys()
    for i,p in enumerate(sorted(seen),1):
        assert checked_cut(p)==tetrahedron_bank()[p]
        if i%128==0:
            coefficients.cache_clear()


def cuts(weights):
    total=[F(),F()]
    for powers,coef in weights.items():
        for i,value in enumerate(tetrahedron_bank()[powers]):
            total[i]+=coef*value
    return tuple(total)


def shape_sum(marks,**kinematics):
    data=graph_data(marks)
    adj=data[3]
    names=tuple(f't{i}' for i in range(14))
    qq=fmpq_mpoly_ctx.get(names)
    bases=[[] for _ in range(4)]
    voltages={}
    for order in permutations(range(5)):
        ds,v=polynomials(data,order,original=order==tuple(range(5)),**kinematics)
        coords=[]
        for r,d in enumerate(ds):
            rem=qq.from_dict(d.to_dict())
            for j,(pivot,base) in enumerate(bases[r]):
                mul=rem[pivot]
                if mul:
                    rem-=mul*base
                    coords.append((r,j,mul))
            if rem:
                pivot,lead=next(iter(rem.terms()))
                coords.append((r,len(bases[r]),lead))
                bases[r].append((pivot,rem/lead))
        pot=qq.from_dict(v.to_dict())
        for x,y in combinations_with_replacement(coords,2):
            key=(*x[:2],*y[:2])
            voltages[key]=voltages.get(key,0)+pot*(x[2]*y[2]*(1+(x[:2]!=y[:2])))
    # Combine the rational coefficient polynomials before the local convergence tests.
    total=qq.constant(0)
    for (i,k,j,h),v in voltages.items():
        if v:
            total+=bases[i][k][1]*bases[j][h][1]*qq.from_dict(adj[i][j].to_dict())*v
    assert all(c.denominator==1 for _,c in total.terms())
    integer=data[1].from_dict({m:int(c) for m,c in total.terms()})
    return cuts(reduce_polynomial(data,integer))


if __name__=='__main__':
    if '--verify-bank' in sys.argv:
        verify_bank()
        print('verified both cuts of all 3872 tetrahedron moments',flush=True)
    marks=(0,0,8,8,8)
    data=graph_data(marks)
    identity=tuple(range(5))
    assert not reduce_full(marks,identity)
    ds,v=polynomials(data,(0,2,1,3,4),original=True)
    source=numerator_poly(data,ds,v)
    weights=reduce_polynomial(data,source)
    assert (len(source),len(weights))==(350728,185)
    assert cuts(weights)==(F(27051556137099721,1903058507750400000),-F(2480717909,209772763200))
    print('single retained order',cuts(weights))
    if '--orbit' in sys.argv:
        targets={(0,0,8,8,8):(-F(7652526749891,33713479800000),F(1403408,7432425)),
                 (0,0,0,8,8):(F(14566533501799,58998589650000),-F(4579976,22297275))}
        for shape,expected in targets.items():
            result=shape_sum(shape)
            assert result==expected
            print('120 ordered normal assignments',shape,*result,flush=True)
