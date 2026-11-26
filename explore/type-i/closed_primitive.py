"""Factor the two-loop six-plus-two trace component through the gauge variation of a closed field."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys
import sympy as sp

from open_ce import coordinates,differential
from open_theta import HERE,parity,reduce_minor,straight


def add_brackets(output,words,c):
    if all(len(set(w))==5 for w in words):
        key=tuple(sorted(tuple(sorted(w)) for w in words))
        output[key]+=c*parity(words[0])*parity(words[1])


def substitute(polynomial,mapping):
    result=Counter()
    for (u,v),c in polynomial.items():
        options=[mapping.get(i,((i,1),)) for i in u+v]
        for assignment in product(*options):
            labels,weights=zip(*assignment)
            coefficient=c
            for x in weights:coefficient*=x
            add_brackets(result,(labels[:5],labels[5:]),coefficient)
    return {k:v for k,v in result.items() if v}


def to_standard(polynomial):
    result=Counter()
    for (u,v),c in polynomial.items():
        for i,w in reduce_minor(u):
            for j,z in reduce_minor(v):
                for key,s in straight(i,j):
                    result[key]+=c*w*z*s
    return {k:v for k,v in result.items() if v}


def lift_pair(row):
    original={}
    for a,b,c in row['terms']:
        brackets=[tuple(k for k in range(7) if k not in pair) for pair in (a,b)]
        original[tuple(brackets)]=F(c)*parity((*brackets[0],*a))*parity((*brackets[1],*b))
    # Labels 7,8,9 denote q=p0+p1, and the two arguments a,b of a closed two-form a wedge b.
    relative=substitute(original,{0:((8,1),),1:((7,1),(8,-1))})
    derivative=Counter()
    for (u,v),c in relative.items():
        labels=u+v
        if labels.count(8)!=1:
            continue
        degree=labels.count(7)
        for i,k in enumerate(labels):
            if k==7:
                term=list(labels)
                term[i]=9
                add_brackets(derivative,(term[:5],term[5:]),c/(degree+1))
    swapped=substitute(derivative,{8:((9,1),),9:((8,1),)})
    derivative.subtract(swapped)
    primitive={k:v for k,v in derivative.items() if v}
    projected=substitute(primitive,{7:((0,1),(1,1)),8:((0,1),),9:((1,1),)})
    assert to_standard(projected)==to_standard(original)
    return [[list(u),list(v),str(c)] for (u,v),c in sorted(primitive.items())]


def preimages(data):
    sources=[s for s in coordinates() if s[0] in ((3,4),(2,2,3))]
    columns=[differential(*s) for s in sources]
    targets=[{} for _ in range(3)]
    for row in data:
        lengths=tuple(row['lengths'])
        if lengths in ((4,4),(2,2,4)):
            targets[row['rank']].update({(lengths,(tuple(a),tuple(b))):F(c) for a,b,c in row['terms']})
    keys=sorted(set().union(*map(set,[*columns,*targets])))
    matrix=sp.Matrix([[c.get(k,0) for c in columns] for k in keys])
    target=sp.Matrix([[c.get(k,0) for c in targets] for k in keys])
    _,ids=matrix.rref()
    basis=matrix[:,list(ids)]
    _,pivots=basis.T.rref()
    x=basis[list(pivots),:].inv()*target[list(pivots),:]
    assert basis*x==target
    return [[rank,[list(item) for item in sources[s]],str(x[i,rank])]
            for i,s in enumerate(ids) for rank in range(3) if x[i,rank]]


def main():
    data=json.loads((HERE/'open_theta.json').read_text())['sum']
    result={'pure_gauge_preimages':preimages(data),
            'closed_two_form_preimages':{str(row['rank']):lift_pair(row) for row in data if row['lengths']==[2,6]}}
    print('seven-field preimages',len(result['pure_gauge_preimages']),
          'closed two-form polynomial sizes',[len(p) for p in result['closed_two_form_preimages'].values()])
    path=HERE/'multitrace_preimages.json'
    if '--save' in sys.argv:path.write_text(json.dumps(result,separators=(',',':'))+'\n')
    assert result==json.loads(path.read_text())


if __name__=='__main__':
    main()
