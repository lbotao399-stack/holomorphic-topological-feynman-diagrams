"""Quadratic tangential variation of the seven-planar bare graph density."""
from fractions import Fraction
from itertools import combinations, permutations
from math import prod
from pathlib import Path
import sys

from flint import fmpq_mat, fmpz_mat

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from point_graph import solve_graph
from square_three import G
from quartic_planar_cores import cubic_cores

NORMAL=((1,2,0),(3,-1,0),(-2,1,1),(0,2,2),(-2,-4,-3))


def fixture(core,paths):
    edges=[]
    for e,(a,b) in enumerate(core):
        chain=[a,*[6+i for i,h in enumerate(paths) if e==h],b]
        edges.extend(zip(chain,chain[1:]))
    return edges


def density(edges,times,normal=NORMAL,slopes=(1,-1,0,0,0),orders=None):
    cycle=G(11,edges).C
    cycle_mat=fmpz_mat(cycle)
    et=fmpz_mat([[c*t for c,t in zip(row,times)] for row in cycle])
    gram=et*cycle_mat.transpose()
    u=int(gram.det())
    adj=gram.inv()*u
    source=[(0,)*5]*6+[tuple(int(i==j) for i in range(5)) for j in range(5)]
    source[0]=(-1,)*5
    flow,potential=solve_graph(edges,source,times,0)
    flow=[[int(u*x) for x in row] for row in flow]
    voltage=[[int(u*x) for x in potential[j]] for j in range(6,11)]
    total=Fraction()
    runs=[]
    for order in (permutations(range(5)) if orders is None else orders):
        q=[normal[j] for j in order]
        f=[[[sum(x*v[i] for x,v in zip(e,q))*c for e,c in zip(flow,row)] for row in cycle] for i in range(3)]
        y=[sum(x*slopes[j] for x,j in zip(e,order)) for e in flow]
        coef=[]
        for row in cycle:
            coef.append(int(fmpz_mat([[1]*14,[a*b for a,b in zip(row,y)],*f[0],*f[1],*f[2]]).det()))
        b=fmpz_mat([coef])
        signal=int((b*adj*b.transpose())[0,0])
        i,j=order.index(0),order.index(1)
        delta=sum((x-y)*v[2] for x,y,v in zip(voltage[i],voltage[j],q))
        runs.append(Fraction(signal*delta,u**27*sum(times)**2))
        total+=runs[-1]
    return runs,total


def original_edges_check(edges,times):
    cycles=G(11,edges).C
    m=fmpz_mat([[sum(a*b*t for a,b,t in zip(x,y,times)) for y in cycles] for x in cycles])
    u=int(m.det())
    adj=m.inv()*u
    # s is the derivative of the one tangential component, supported at normal slots.
    slopes=(1,-1,0,0,0)
    source=[(0,)*4]*6+[(z,*p) for z,p in zip(slopes,NORMAL)]
    flow,_=solve_graph(edges,source,times,0)
    y=[[int(u*v) for v in q] for q in flow]
    normals=[[c*q[j] for c,q in zip(cs,y)] for j in (1,2,3) for cs in cycles]
    b=[int(fmpz_mat([[1]*14,[c*q[0] for c,q in zip(cs,y)],*normals]).det()) for cs in cycles]
    bc=fmpz_mat([b]).transpose()

    def columns(ids):
        return [tuple(row[i] for row in cycles) for i in ids]
    four=three=0
    candidates=[]
    for ids in combinations(range(14),4):
        minor=int(fmpz_mat(columns(ids)).det())
        four+=prod(times[i] for i in ids)*minor**2
        if minor:
            candidates.append((ids,minor))
    assert four==u

    flags=[]
    for ids in combinations(range(14),3):
        cof=[(-1)**j*int(fmpz_mat([q[:j]+q[j+1:] for q in columns(ids)]).det()) for j in range(4)]
        contracted=sum(x*z for x,z in zip(cof,b))
        three+=prod(times[i] for i in ids)*contracted**2
        if contracted:
            flags.append((ids,contracted))
    assert three==int((bc.transpose()*adj*bc)[0,0]) and three
    # Restore seven actual five-dimensional waves and take the linear coefficient of the signed graph cofactor.
    triple,exterior=flags[0]
    quartet,volume=next((ids,z) for ids,z in candidates if set(ids).isdisjoint(triple))
    places=triple+quartet
    eps=(1,1,-2,1,1,-1,-1)
    planar=[(q,0,0,0,0) if i<3 else (0,q,0,0,0) for i,q in enumerate(eps)]
    full_edges,full_t,parents=[],[],[]
    for e,(a,c) in enumerate(edges):
        if e in places:
            i=places.index(e)
            v=11+i
            full_edges.extend(((a,v),(v,c)))
            full_t.extend((Fraction(times[e],3),Fraction(2*times[e],3)))
            parents.extend((e,e))
        else:
            full_edges.append((a,c))
            full_t.append(times[e])
            parents.append(e)
    values=[]
    for lift in (-1,0,1,2):
        src=[(0,)*5]*6+[(lift*z,0,*p) for z,p in zip(slopes,NORMAL)]+planar
        flow,_=solve_graph(full_edges,src,full_t,0)
        rows=[[str(cs[e]*q[j]) for e,q in zip(parents,flow)] for j in range(5) for cs in cycles]
        values.append(Fraction(str(fmpq_mat([[1]*21,*rows]).det())))
    assert values[0]+values[2]==2*values[1] and values[1]+values[3]==2*values[2]
    assert abs(values[2]-values[1])==abs(Fraction(volume*exterior*prod(eps),u**13))


def kinematic_kernel():
    # Three 5-minors carry pairs of light slots, the fourth one light slot.
    missing=((0,1),(0,2),(3,4))
    heavy_groups=[tuple(i for i in range(5) if i not in pair) for pair in missing]
    heavy_groups.append((0,1,2,4))
    counts=[sum(i in g for g in heavy_groups) for i in range(5)]
    assert min(counts)>=2 and 7+min(counts)>2*4
    z=(1,-1,0,0,0)
    d=int(fmpz_mat([(a,*p) for a,p in zip(z,NORMAL)][:4]).det())
    outside=0
    for order in permutations(range(5)):
        a,b=order[:2]
        charge=(int(a==0)-int(a==1))*NORMAL[b][2]+(int(b==0)-int(b==1))*NORMAL[a][2]
        outside+=charge*prod(int(fmpz_mat([NORMAL[order[i]] for i in group]).det())**2 for group in heavy_groups[:3])
    assert outside==29701056 and d*d==16
    plane=((1,0),(1,0),(-2,0),(0,1),(0,1),(0,-1),(0,-1))
    def br(x,y):
        return x[0]*y[1]-x[1]*y[0]
    inside=sum(prod(br(x,y)**2 for x,y in zip(p[::2],(*p[1::2],(1,0)))) for p in permutations(plane))
    assert inside==4608 and inside*d*d*outside==2189799456768


def independent_models():
    q=((1,0,0),(0,1,0),(0,0,-3),(-4,Fraction(-11,3),8),(3,Fraction(8,3),-5))
    assert tuple(map(sum,zip(*q)))==(0,0,0)
    pairs=((0,1),(0,2),(0,3))
    minors={(i,j,k):Fraction(str(fmpq_mat([list(map(str,q[h])) for h in (i,j,k)]).det()))
            for i in range(5) for j in range(5) for k in range(5)}
    assert all(x.denominator==1 for x in minors.values())
    rows=[]
    for omissions,singleton,marked in (
        (((0,1),(0,1),(2,3)),2,(0,1)),
        (((0,1),(0,1),(2,3)),2,(0,2)),
        (((0,1),(0,2),(1,2)),3,(0,1)),
    ):
        assert all(sum(k in m for m in omissions)+int(k==singleton)<=2 for k in range(5))
        values=[0]*3
        for order in permutations(range(5)):
            weight=prod(minors[tuple(order[k] for k in range(5) if k not in miss)]**2 for miss in omissions)
            a,b=(order[h] for h in marked)
            for col,(i,j) in enumerate(pairs):
                charge=(int(a==i)-int(a==j))*minors[i,j,b]+(int(b==i)-int(b==j))*minors[i,j,a]
                values[col]+=int(weight*charge)
        rows.append(values)
    assert rows==[[6083520,-3024384,-133056],[-1750800,8443776,-4291056],[-1704960,-3091968,6386688]]
    assert int(fmpz_mat(rows).det())==188775118544044032000


if __name__=='__main__':
    core=cubic_cores()[-2]
    times=[2+(i*3)%11 for i in range(14)]
    edges=fixture(core,(0,2,3,6,8))
    original_edges_check(edges,times)
    kinematic_kernel()
    independent_models()
    values,total=density(edges,times)
    assert all(values)
    assert total==Fraction(21356867491118736770664566300207257078494806450,5479804936105143233773111885248362676926355573)
    print('fixed-time normal-label sum',total)
    # The one positive and one negative charge force a positive voltage drop between the marked interior slots.
    normal=((1,1,0),(2,-1,0),(-1,2,1),(1,0,-1),(-3,-2,0))
    edges=fixture(cubic_cores()[-1],(0,0,0,1,3))
    values,total=density(edges,times,normal,orders=[(2,0,1,3,4)])
    assert total==Fraction(29382104764424983171591305241051193671680000,12126555610398186200391867570813863461438965986513)
    print('single labeled positive density',total)
