"""Polarize the soft-pair cofactor over the four permutable hard charges."""
from argparse import ArgumentParser
from fractions import Fraction
from collections import Counter
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import json
from time import perf_counter

from flint import fmpq, fmpq_mat
from mixed_bare_primitive_pole import PrimitiveChamber
from mixed_bare_ranks import WAVES
from mixed_unordered_hard import CTX, X, SOFT, grounded_response, integrate_unit

PARTITIONS=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
HARD=(*WAVES[:3],WAVES[5])


def common_minors(matrix):
    """The 45 maximal minors of the ten remaining rows and eight mean columns."""
    previous={0:CTX.constant(1)}
    # Columns can be exchanged together; their common orientation disappears
    # from the squared cofactor.
    columns=sorted(zip(*matrix),key=lambda col:sum(bool(x) for x in col))
    for column in columns:
        current={}
        for mask,value in previous.items():
            for i,x in enumerate(column):
                if not x or mask&(1<<i):continue
                final=mask|(1<<i)
                add=value*x
                if (mask>>(i+1)).bit_count()%2:add=-add
                if final in current:current[final]+=add
                else:current[final]=add
        previous=current
    full=(1<<len(matrix))-1
    return [(i,j,(-1)**(i+j+1)*previous[full^(1<<i)^(1<<j)])
            for i,j in combinations(range(len(matrix)),2)
            if previous.get(full^(1<<i)^(1<<j))]


class HardFrame:
    def __init__(self,probe,marks,times,z_sites):
        self.jacobian=prod(times[e] for e in marks)
        self.soft_scale=fmpq(2,sum(times)**2)
        cycles=probe.cycles
        currents,green=grounded_response(tuple(probe.core),tuple(times))
        self.four=tuple(i for i in range(6) if i not in z_sites)
        waves=[WAVES[3]]*6
        for i,p in zip(self.four,HARD,strict=True):waves[i]=p
        def endpoint(row,i):
            a,b=probe.core[marks[i]]
            return fmpq(str(row[a]))+X[i]*fmpq(str(row[b]-row[a]))
        flow=[[sum((waves[i][mu]*endpoint(row,i) for i in range(6) if waves[i][mu]),CTX.constant(0))
               for mu in range(4)] for row in currents]
        means=[[CTX.constant(1)]*9]+[[row[e]*flow[e][mu] for e in range(9)] for mu in range(4) for row in cycles]
        jumps=[[0]*6]+[[row[e]*p[mu] for e,p in zip(marks,waves)] for mu in range(4) for row in cycles]
        echelon,rank=fmpq_mat(jumps).transpose().rref()
        self.minors=[]
        if rank<6:return
        selected=[next(j for j in range(17) if echelon[i,j]) for i in range(6)]
        pivot=fmpq_mat([jumps[i] for i in selected])
        self.determinant=pivot.det()
        remaining=[i for i in range(17) if i not in selected]
        transform=fmpq_mat([jumps[i] for i in remaining])*pivot.inv()
        assert remaining[0]==0
        reduced=[[means[i][j]-means[i][0]-sum(transform[k,l]*(means[s][j]-means[s][0])
                          for l,s in enumerate(selected) if transform[k,l]) for j in range(1,9)]
                 for k,i in enumerate(remaining) if k]
        self.minors=common_minors(reduced)
        if not self.minors:return
        def soften(direction):
            rows=[[0]*4]+[[int(k==i)*v for k in range(4)] for v in direction for i in range(4)]
            return (fmpq_mat([rows[i] for i in remaining])-transform*fmpq_mat([rows[i] for i in selected])).tolist()[1:]
        self.projections=[soften(p) for p in HARD]
        gram=[[sum(t*x*y for t,x,y in zip(times,a,b)) for b in cycles] for a in cycles]
        self.metric=[[gram[a][c]*gram[b][d]-gram[a][d]*gram[b][c] for c,d in SOFT] for a,b in SOFT]
        def site_green(i,j):
            a,b=probe.core[marks[i]]
            value=(1-X[i])*endpoint(green[a],j)+X[i]*endpoint(green[b],j)
            if i==j:value+=times[marks[i]]*X[i]*(1-X[i])
            return value
        self.resistance={(i,j):site_green(i,i)+site_green(j,j)-2*site_green(i,j) for i,j in combinations(range(6),2)}
        self.z=z_sites
        self.bridges={(i,j):times[marks[i]] for i,j in combinations(range(6),2) if marks[i]==marks[j]}

    def cofactor(self,pair):
        i,j=pair
        block=[[x-y for x,y in zip(a,b)] for a,b in zip(self.projections[i],self.projections[j])]
        return [self.determinant*sum((poly*(block[r][a]*block[s][b]-block[r][b]*block[s][a])
                 for r,s,poly in self.minors if block[r][a]*block[s][b]!=block[r][b]*block[s][a]),CTX.constant(0))
                for a,b in SOFT]

    def inner(self,left,right):
        return self.soft_scale*sum((x*sum((c*y for c,y in zip(row,right) if c and y),CTX.constant(0))
                                    for x,row in zip(left,self.metric) if x),CTX.constant(0))

    def voltage(self,pair,complement):
        i,j=(self.four[k] for k in pair)
        k,l=(self.four[r] for r in complement)
        coeff=Counter()
        def add(m,n,v):coeff[tuple(sorted((m,n)))]+=v
        add(k,l,4)
        for t in (i,j):
            for z in self.z:add(t,z,1)
        for t in (k,l):
            for z in self.z:add(t,z,-1)
        for t in (i,j):
            for u in (k,l):add(t,u,-1)
        poly=sum(c*self.resistance[p] for p,c in coeff.items())/4
        bridge={p:-fmpq(t*coeff[p],2) for p,t in self.bridges.items() if coeff[p]}
        return poly,bridge

    def channels(self):
        even=CTX.constant(0)
        odd=CTX.constant(0)
        bridges=({}, {})
        if not self.minors:return (even,odd),bridges
        voltage_sum=CTX.constant(0)
        bridge_sum=Counter()
        for pair,complement in PARTITIONS:
            v,vr=self.voltage(pair,complement)
            w,wr=self.voltage(complement,pair)
            voltage_sum+=v+w
            bridge_sum.update(vr)
            bridge_sum.update(wr)
            a=self.cofactor(pair)
            b=self.cofactor(complement)
            symmetric=[(x+y)/2 for x,y in zip(a,b)]
            alternating=[(x-y)/2 for x,y in zip(a,b)]
            gram_even=self.inner(symmetric,symmetric)+self.inner(alternating,alternating)
            gram_odd=2*self.inner(symmetric,alternating)
            even+=(v+w)*gram_even
            odd+=(v-w)*gram_odd
            for p in vr.keys()|wr.keys():
                left,right=vr.get(p,0),wr.get(p,0)
                for data,c,poly in zip(bridges,(left+right,left-right),(gram_even,gram_odd)):
                    if c:data[p]=data.get(p,CTX.constant(0))+c*poly
        if voltage_sum or any(bridge_sum.values()):raise ArithmeticError('Trivial hard permutation channel survived')
        return (even,odd),bridges


def integrated(poly,bridges):
    result=fmpq(integrate_unit(poly))
    for (i,j),kernel in bridges.items():
        for ex,c in kernel.terms():
            a,b=int(ex[i]),int(ex[j])
            denominator=prod(int(d)+1 for k,d in enumerate(ex) if k!=i and k!=j)
            result+=c/fmpq(denominator*(a+2)*(b+2)*(a+b+3))
    return result


def placement(word,times,progress=False):
    marks=tuple(map(int,word))
    if len(marks)!=6 or any(i not in range(9) for i in marks):raise ValueError('Expected six core path labels')
    probe=PrimitiveChamber(4,marks)
    stabilizer=sum(tuple(sorted(g[e][0] for e in marks))==marks for g in probe.actions)
    result=[fmpq(0),fmpq(0)]
    repeated=prod(map(factorial,Counter(marks).values()))
    for number,z in enumerate(combinations(range(6),2),1):
        start=perf_counter()
        frame=HardFrame(probe,marks,times,z)
        channels,bridges=frame.channels()
        for i,(poly,correction) in enumerate(zip(channels,bridges,strict=True)):
            result[i]+=8*frame.jacobian*integrated(poly,correction)/(stabilizer*repeated)
        if progress:print('Z pair',number,z,'common complementary minors',len(frame.minors),'seconds',round(perf_counter()-start,3),flush=True)
    return tuple(map(str,result))


def all_unit():
    records=json.loads(Path(__file__).with_name('mixed_bipartite_unit_metric.json').read_text())
    totals=[Fraction(),Fraction()]
    channels=Counter()
    for word,record in sorted(records.items()):
        values=tuple(map(Fraction,placement(word,(1,)*9)))
        if sum(values)!=Fraction(record['density']):
            raise ArithmeticError(f'Original label-summed density differs at {word}')
        for i,v in enumerate(values):totals[i]+=v
        channels.update({len(set(word)):1})
    print('Complete fixed-metric occupancies, complementary sectors and sum:',dict(channels),*totals,sum(totals),flush=True)
    return totals


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('word',nargs='?',default='012345')
    parser.add_argument('--times',nargs=9,type=int,default=(1,)*9)
    parser.add_argument('--progress',action='store_true')
    parser.add_argument('--all-unit',action='store_true',help='compare all 54 bipartite occupancies')
    args=parser.parse_args()
    if args.all_unit:all_unit()
    else:
        values=placement(args.word,args.times,args.progress)
        print('Even and odd complementary channels:',*values,'total',sum(map(Fraction,values)),flush=True)
