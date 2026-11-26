"""Sum polynomial wheel cuts using a common angular momentum on one rim edge."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import factorial,prod
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import json

from mixed_bipartite_numerator_bank import (
    HERE,certified,compile_program,completed,direct_period,height,prism_charts,reflected,
    geometry,paths,subtotal,unpack,zeta_interval,
)


def divergent_cuts(key,q,degree):
    """List all failed collision inequalities for the two extreme scalar terms."""
    spokes,ring=unpack(key)
    result=[]
    for end in (0,2):
        graph=[[0]*6 for _ in range(6)]
        for v in range(4):
            graph[v][4]=graph[4][v]=spokes[v]-degree*int(v==end)
            w=(v-1)%4
            graph[v][w]=graph[w][v]=ring[v]
        for v in range(5):
            graph[v][5]=graph[5][v]=2*q-sum(graph[v])
        assert sum(graph[5])==2*q
        result.append(tuple((vertices,excess) for size in (2,3) for vertices in combinations(range(6),size)
                            if (excess:=sum(graph[i][j] for i,j in combinations(vertices,2))-(size-1)*q)>=0))
    return tuple(result)



def verify_chart(key,q,d):
    """Verify the unsplit signed diagram, and every pole used in its radial sums."""
    u,b=unpack(key)
    graph=[[0]*6 for _ in range(6)]
    for i in range(4):
        graph[4][i]=graph[i][4]=u[i]
        graph[i][(i-1)%4]=graph[(i-1)%4][i]=b[i]
    graph[0][2]=graph[2][0]=-d
    for i in range(5):
        graph[i][5]=graph[5][i]=2*q-sum(graph[i])
    assert sum(graph[5])==2*q
    assert max(map(max,graph))<q
    assert all(graph[i][j]+graph[j][k]+graph[k][i]<2*q for i,j,k in combinations(range(6),3))
    for r,j in paths(q,d):
        for mask in ((6,9) if j else (0,)):
            for vs,cross,internal in geometry()[1].values():
                extra=r*(0 in vs)+(2*d-r)*(2 in vs)+sum(j*((mask>>i)&1) for i in cross)
                assert extra%2==0
                low=sum(q-u[i] for i in vs)-sum(b[i] for i in internal)+extra//2
                assert 1<=low and low+sum(q-1-b[i] for i in cross)<=2*q+d-3

def evaluate(q,d,keys,exe,tmp,prime_workers=1):
    return certified(q,keys,exe,tmp,degree=d,proof=height(q,d,displacement=d),fixed_spin=True,prime_workers=prime_workers)


def controls(complete,wheel,tmp):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    inputs=[]
    targets={}
    for ex,cols in zip(ref['exponents'],ref['coefficients'],strict=True):
        seed=tuple(x+1 for x in ex)
        inputs.append((*seed,1,0))
        values=tuple(F(str(v))*factorial(5)/prod(map(factorial,ex)) for v in cols)
        assert targets.setdefault(min(prism_charts(completed(seed,6),6)),values)==values
    for d in range(1,6):
        selected=reflected(complete,6,inputs,d,formal=True)
        wanted={}
        actual_divergent=0
        for seed,(_,key,scale) in selected.items():
            values=tuple(v/scale for v in targets[seed])
            assert wanted.setdefault(key,values)==values
            actual_divergent+=any(divergent_cuts(key,6,d))
        for key in wanted:verify_chart(key,6,d)
        if wanted:
            assert evaluate(6,d,wanted,wheel,tmp)==wanted
        print('Independent q=6 seeds, including divergent monomials:',d,len(selected),actual_divergent,flush=True)


def physical(complete,wheel,tmp,check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as stream:
        data={tuple(map(int,r[:9])):r[:11] for line in stream if (r:=line.split())[11]=='0'}
    for kind in ('wheel',1,2,3,4,5):
        old=reflected(complete,13,data.values(),kind)
        for seed in old:del data[seed]
    assert len(data)==4769
    total=list(subtotal())
    for d,expected in enumerate((87,232,271,218,109,12),1):
        selected=reflected(complete,13,data.values(),d,formal=True)
        assert len(selected)==expected
        for seed in selected:del data[seed]
        cols=[v for v in selected.values()]
        keys={k for _,k,_ in cols}
        excess=Counter()
        for _,key,_ in cols:
            verify_chart(key,13,d)
            a,b=divergent_cuts(key,13,d)
            assert a==b and len(a)==1 and a[0][0]==(4,5) and 0<=a[0][1]<d
            excess[a[0][1]]+=1
        values=evaluate(13,d,keys,wheel,tmp)
        if check:
            key=min(keys,key=lambda k:(prod(13-v for v in unpack(k)[1]),k))
            assert direct_period(*unpack(key),13,d,fixed_spin=True)==values[key]
        add=tuple(sum(w*fac*values[key][j] for w,key,fac in cols) for j in range(3))
        total=[a+b for a,b in zip(total,add,strict=True)]
        print('degree, excess counts, additional three columns:',d,sorted(excess.items()),*add,flush=True)
    assert len(data)==3840
    lo=hi=total[0]
    for v,n in zip(total[1:],(3,5),strict=True):
        a,b=sorted(v*z for z in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    assert total==[F(1875840319482653063081961833881,8933882473602043013726208000000000),
                   F(192571615355536381505743,342912864219655137408000000),
                   -F(119708316672311,140256167203152000)]
    assert -198448<10**20*lo<=10**20*hi<-198447
    print('33551 components:',*total,flush=True)
    print('strict bounds, units of 10^-20:',(10**20*lo)//1,-((-10**20*hi)//1),flush=True)
    return total


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--controls-only',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        complete,wheel=tmp/'complete',tmp/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_wheel_bank.cpp',wheel)
        if args.check or args.controls_only:controls(complete,wheel,tmp)
        if not args.controls_only:physical(complete,wheel,tmp,args.check)
