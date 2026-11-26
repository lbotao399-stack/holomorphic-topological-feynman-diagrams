"""Reconstruct the primitive bare numerator from the exact quadratic grid R_8=1."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, gcd, prod
from pathlib import Path
import json
import subprocess
import sys

from chain_moments import Reduction, graphs, parity, sp

here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parents[1]/'code'))
import ribbon

normalizer=4*factorial(8)
generic=[(2,1),(1,4),(-2,3),(-1,-8)]
probe=[(1,0),(0,1),(1,1),(-2,-2)]


def pack(red,sources):
    assert len(sources)==4 and tuple(map(sum,zip(*sources)))==(0,0)
    out=[]
    def append(row):
        out.extend(str(int(x)) for x in row)
    for j in range(9):
        append(red.C[:,j])
    charges=[[0]*2 for _ in red.core]
    offsets=[[[0]*5 for _ in range(2)] for _ in red.core]
    vertical=[]
    for i,(e,j,v) in enumerate(red.interior):
        charges[e]=[a+b for a,b in zip(charges[e],sources[v-6])]
        vertical.append((e,tuple(charges[e])))
        for k in range(2):
            offsets[e][k][i+1]=charges[e][k]
    for k in range(2):
        injections=[sum(charges[e][k] for e,(_,v) in enumerate(red.core) if v==j) for j in range(5)]
        for e,n in zip(red.tree,red.tree_inverse*sp.Matrix(injections)):
            offsets[e][k][0]=int(n)
    for row in offsets:
        for u in row:
            append(u)
    append([(-1)**sum((8-e)*len(g) for e,g in enumerate(red.groups))])
    pairs=list(combinations_with_replacement(range(5),2))
    for left in pairs:
        for right in pairs:
            moment=factorial(8)*red.moment(tuple((left+right).count(i) for i in range(1,5)))
            assert moment.denominator==1 and 0<=moment<=factorial(8)
            append([moment])
    terms=[]
    for J in combinations(range(8),4):
        if not red.minors[J]:
            continue
        H=[e for e in range(8) if e not in J]
        labels=H+[e for e,_ in vertical]
        coef=[0]*16
        for S in combinations(range(8),4):
            T=tuple(i for i in range(8) if i not in S)
            d=red.minor([labels[i] for i in S])*red.minor([labels[i] for i in T])
            if d:
                coef[sum(1<<i for i in T if i<4)]+=parity(S+T)*d*prod(p[int(4+i in T)] for i,(_,p) in enumerate(vertical))
        if any(coef):
            terms.append(list(J)+H+[red.minors[J]*parity(list(J)+H)]+coef)
    append([len(terms)])
    for term in terms:
        append(term)
    return ' '.join(out)


def interpolate(values):
    assert len(values)==3**8
    for axis in range(8):
        stride=3**axis
        for start in range(0,len(values),3*stride):
            for off in range(stride):
                i=start+off
                a,b,c=(values[i+j*stride] for j in range(3))
                assert (a+c)%2==0
                quadratic=(a-2*b+c)//2
                linear=b-a-3*quadratic
                values[i],values[i+stride],values[i+2*stride]=a-linear-quadratic,linear,quadratic
    coeff={}
    for powers,v in zip(product(range(3),repeat=8),values):
        if v:
            last=9-sum(powers)
            assert last in (0,1,2)
            coeff[powers+(last,)]=v
    return coeff


def reconstruct(red,sources,binary):
    data=pack(red,sources)
    p=subprocess.run([str(binary)],input=data,text=True,capture_output=True,check=True)
    return interpolate([int(x) for x in p.stdout.split()])


def check(red,sources,coeff,R):
    reduced,U=red.density(R,sources)
    expected=reduced*U**5/prod(R[e] for e,_,_ in red.interior)*normalizer
    value=sum(c*prod(r**j for r,j in zip(R,powers)) for powers,c in coeff.items())
    assert value==expected,(value,expected)


def serialize(index,red,coeff,denominator,**metadata):
    coeff={t:c for t,c in coeff.items() if c}
    common=gcd(denominator,*coeff.values())
    return {'graph':index,'edges':red.core,'h':[len(g) for g in red.groups],
            'denominator':denominator//common,'terms':[[list(t),c//common] for t,c in sorted(coeff.items())],
            **metadata}


def colour_data(red):
    edges=[edge for ch in red.chains for edge in zip(ch,ch[1:])]
    tensor=ribbon.colour_weight(10,edges,list(range(6,10)))
    modes={}
    raw=[]
    for perm in permutations(range(4)):
        ranks=defaultdict(int)
        for cycles,colour in tensor.items():
            ordered=ribbon.canon_tau([tuple(perm[i] for i in c) for c in cycles])
            if ordered==((0,1,2,3),):
                for k,v in colour.items():
                    assert v.denominator==1
                    ranks[k]+=parity(perm)*int(v)
        ranks={k:v for k,v in ranks.items() if v}
        if not ranks:
            continue
        sources=[probe[i] for i in perm]
        raw.append((sources,ranks))
        minors=[sources[i][0]*sources[j][1]-sources[i][1]*sources[j][0] for i,j in ((0,1),(0,2),(1,2))]
        scale=gcd(*minors)
        if next(i for i in minors if i)<0:
            scale=-scale
        # Any two such rank-two source lists with proportional minors are related by GL(2).
        # The two sets of four holomorphic forms give the fourth determinant power.
        key=tuple(i//scale for i in minors)
        if key not in modes:
            modes[key]=(sources,scale,defaultdict(int))
        for rank,v in ranks.items():
            modes[key][2][rank]+=v*scale**4
    return modes,raw


def cyclic(index,red,binary):
    modes,raw=colour_data(red)
    by_rank=defaultdict(lambda:defaultdict(F))
    for sources,scale,ranks in modes.values():
        if not any(ranks.values()):
            continue
        coeff=reconstruct(red,sources,binary)
        for rank,v in ranks.items():
            for powers,c in coeff.items():
                by_rank[rank][powers]+=F(v*c,scale**4)
    for rank,poly in by_rank.items():
        assert all(c.denominator==1 for c in poly.values())
    by_rank={k:{t:int(c) for t,c in poly.items() if c}
             for k,poly in by_rank.items() if any(poly.values())}
    # Off-grid replay without combining any projectively equivalent momentum assignments.
    R=[2,3,5,7,11,13,17,19,23]
    reference=defaultdict(F)
    for sources,ranks in raw:
        density,U=red.density(R,sources)
        for k,v in ranks.items():
            reference[k]+=v*density*U**5*normalizer/prod(R[e] for e,_,_ in red.interior)
    for k in set(reference)|set(by_rank):
        val=sum(c*prod(r**n for r,n in zip(R,t)) for t,c in by_rank.get(k,{}).items())
        assert val==reference[k],(index,k,val,reference[k])
    return [serialize(index,red,poly,normalizer*graphs[index]['aut'],N_power=k,automorphism_order=graphs[index]['aut'])
            for k,poly in sorted(by_rank.items())]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--index',type=int,default=3)
    parser.add_argument('--cyclic',action='store_true')
    parser.add_argument('--through',type=int)
    args=parser.parse_args()
    if args.cyclic:
        for idx in range(args.index,(args.through if args.through is not None else args.index)+1):
            red=Reduction(graphs[idx])
            rows=cyclic(idx,red,args.binary)
            args.output.mkdir(exist_ok=True)
            (args.output/f'graph{idx:02}.json').write_text(json.dumps(rows,separators=(',',':'))+'\n')
            print(idx,[(r['N_power'],len(r['terms']),r['denominator']) for r in rows],flush=True)
        return
    red=Reduction(graphs[args.index])
    coeff=reconstruct(red,generic,args.binary)
    for R in ([2,3,5,7,11,13,17,19,23],[11,7,5,4,8,2,9,3,6]):
        check(red,generic,coeff,R)
    record=serialize(args.index,red,coeff,normalizer,sources=generic)
    args.output.write_text(json.dumps(record,separators=(',',':'))+'\n')
    print(args.index,len(coeff),record['denominator'],flush=True)


if __name__=='__main__':
    main()
