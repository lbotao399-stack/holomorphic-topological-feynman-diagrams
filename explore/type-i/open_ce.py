"""Exact CE image for all translation- and SL(5)-invariant pure-open seven-field counterterms."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from pathlib import Path
import json
import sympy as sp
from open_theta import HERE,reduce_minor,relabelings,straight


def trace_tags(lengths):
    start=0
    words=[]
    for n in lengths:
        words.append(tuple(range(start,start+n)))
        start+=n
    return words


def split_minor(indices,v):
    entries=[(u,) if u<v else (u,u+1) if u==v else (u+1,) for u in indices]
    return list(product(*entries))


def differential(lengths,first,second):
    """d Alt_7 = -sum_v (-1)^v Alt_8 inserting a product in trace slot v."""
    source=trace_tags(lengths)
    out=Counter()
    for v in range(7):
        words=[]
        for word in source:
            w=[]
            for i in word:
                w.extend([i] if i<v else [i,i+1] if i==v else [i+1])
            words.append(tuple(w))
        words.sort(key=lambda w:(len(w),w))
        target=tuple(map(len,words))
        a,b=split_minor(first,v),split_minor(second,v)
        for perm,sgn in relabelings(tuple(words)):
            c=(-1)**(v+1)*sgn
            for x in a:
                m1=reduce_minor(tuple(perm[i] for i in x))
                for y in b:
                    m2=reduce_minor(tuple(perm[i] for i in y))
                    for i,wi in m1:
                        for j,wj in m2:
                            for term,ws in straight(i,j):
                                out[target,term]+=c*wi*wj*ws
    return {k:v for k,v in out.items() if v}


def coordinates():
    minors=list(combinations(range(6),5))
    for lengths in [(7,),(2,5),(3,4),(2,2,3)]:
        for a,b in combinations_with_replacement(minors,2):
            yield lengths,a,b


def solve():
    targets={}
    for rec in json.loads((HERE/'open_theta.json').read_text())['sum']:
        rank=rec['rank']
        targets.setdefault(rank,{}).update({
            (tuple(rec['lengths']),(tuple(a),tuple(b))):F(c) for a,b,c in rec['terms']})
    sources=list(coordinates())
    table=[differential(*s) for s in sources]
    rows=sorted(set().union(*(set(i) for i in [*table,*targets.values()])))
    M=sp.Matrix([[col.get(k,0) for col in table] for k in rows])
    B=sp.Matrix([[targets[n].get(k,0) for n in range(3)] for k in rows])
    _,independent=M.rref()
    basis=M[:,list(independent)]
    _,eqs=basis.T.rref()
    square=basis[list(eqs),:]
    x=square.inv()*B[list(eqs),:]
    residual=B-basis*x
    witness=((2,6),((0,1),(2,3)))
    assert all(col.get(witness,0)==0 for col in table)
    assert [targets[i].get(witness,0) for i in range(3)]==[F(13,43200),F(1,43200),0]
    assert all(not any(lengths==(2,6) for lengths,_ in col) for col in table)
    assert len(sources)==84 and len(independent)==5 and residual.rank()==3
    print('84 candidates, CE rank 5, residual rank 3 across formal rank powers')
    print('all (2,6)-trace boundary coefficients vanish; raw [01][23] at SO(32) = 1/960')


if __name__=='__main__':
    solve()
