"""Literal adjoint-index contraction for all 40320 external assignments on the orthogonal probe."""
from fractions import Fraction as F
from itertools import combinations,permutations
from math import factorial,lcm,prod

import numpy as np
from open_theta import HERE,SHAPES,det,parity,theta,twoloop
from hexagon import edge_forms
from htk.graph import Graph
import json

EXT=((0,1),(0,1),(2,3),(2,3),(4,5),(5,6),(6,7),(7,4))
P=[[int(i==j) for j in range(5)] for i in range(5)]+[[1,-1,2,0,3],[0,2,-2,1,1]]
P.append([-sum(x) for x in zip(*P)])
ASSIGNMENTS=np.array(list(permutations(range(8))),dtype=np.int16)
SIGNS=np.array([parity(p) for p in ASSIGNMENTS.tolist()],dtype=np.int16)


def numerators(shape):
    g=theta(shape)
    denominator=lcm(*shape)
    weights=[int(denominator*twoloop.mu_weights(g)[i]) for i in range(1,9)]
    minors=np.zeros((8,)*5,dtype=np.int64)
    for s in combinations(range(8),5):
        value=int(det([P[i] for i in s]))
        for p in permutations(s):
            minors[p]=parity(p)*value
    labels=ASSIGNMENTS.T
    def with_mean(ids):
        return sum(w*minors[(labels[i],*(labels[j] for j in ids))] for i,w in enumerate(weights))
    A,B,C=([i-1 for i in g.interior[k]] for k in 'ABC')
    overall=np.zeros(len(ASSIGNMENTS),dtype=np.int64)
    for chosen in combinations(B,4-len(A)):
        top=A+list(chosen)
        bottom=[i for i in B if i not in chosen]+C
        shuffle=[0,*[2+i for i in top],1,*[2+i for i in bottom]]
        overall+=parity(shuffle)*with_mean(top)*with_mean(bottom)
    coefficient=twoloop.c0(shape,5)*(-1)**(shape[1]+1)/denominator**2
    # Integer bound on every product and intermediate sum above.
    assert 10*(sum(weights)*int(abs(minors).max()))**2<2**62
    return overall,coefficient


def literal_colors(n,shape,external_edges=EXT):
    basis=list(combinations(range(n),2))
    idx={b:i+1 for i,b in enumerate(basis)}
    def slot(i,j):
        return idx[(i,j)] if i<j else -idx[(j,i)] if j<i else 0
    def comm(a,b):
        i,j=a
        k,l=b
        return (int(j==k)*slot(i,l)-int(j==l)*slot(i,k)
                -int(i==k)*slot(j,l)+int(i==l)*slot(j,k))
    table=np.zeros((len(basis)+1,len(basis)+1),dtype=np.int16)
    for a,b in permutations(basis,2):
        table[idx[a],idx[b]]=comm(a,b)
    external=np.array([slot(*pair) for pair in external_edges],dtype=np.int16)
    def act(sequence):
        states=np.broadcast_to(np.arange(len(basis)+1,dtype=np.int16),(len(ASSIGNMENTS),len(basis)+1))
        for v in sequence:
            x=external[ASSIGNMENTS[:,v-1,None]]
            states=np.sign(x)*np.sign(states)*table[abs(x),abs(states)]
        return states
    g=theta(shape)
    a,b,c=[act(g.interior[k]) for k in 'ABC']
    weight=np.zeros(len(ASSIGNMENTS),dtype=np.int32)
    for u,v in zip(*np.nonzero(table)):
        uv=table[u,v]
        i,j,k=a[:,u],b[:,v],c[:,abs(uv)]
        ij=table[abs(i),abs(j)]
        weight+=(abs(ij)==abs(k))*np.sign(ij)*np.sign(i)*np.sign(j)*np.sign(k)*int(np.sign(uv))
    # f = -2 sign, the three initial inverse metrics contribute (-1/2)^3,
    # and the odd-slot permutation has sign -1. Therefore weight/2 is the color.
    return weight


def evaluate_terms(entry,momenta):
    pair={}
    for i,j in combinations(range(7),2):
        rest=[k for k in range(7) if k not in (i,j)]
        pair[i,j]=parity([*rest,i,j])*det([momenta[k] for k in rest])
    return sum(F(coefficient)*pair[tuple(a)]*pair[tuple(b)] for a,b,coefficient in entry['terms'])


def stored_witness(record):
    return 4*evaluate_terms(next(d for d in record if d['lengths']==[2,2,4]),P)


def ce_witness(data):
    p=P[:5]+P[:3]
    p.append([-sum(x) for x in zip(*p)])
    cycle=list(range(2,9))
    by_rank=[F(0)]*3
    for i in range(7):
        s=cycle[i:]+cycle[:i]
        fused=[x+y for x,y in zip(p[s[0]],p[s[1]])]
        waves=p[:2]+[fused]+[p[k] for k in s[2:]]
        for row in data:
            if row['lengths']==[2,6]:
                by_rank[row['rank']]-=2*evaluate_terms(row,waves)
    assert sum(32**i*c for i,c in enumerate(by_rank))==F(1,240)
    print('nine-input CE boundary in the (2,7) trace sector',by_rank)


def check_heat(shape):
    g=theta(shape)
    graph=Graph(10,g.edges)
    t=[2+i for i in range(11)]
    y=twoloop.currents(10,g.edges,9,t,[[0]*5,*P,[0]*5])
    cycles=graph.cycle_matrix()
    rows=[[cycle[e]*y[e][i] for e in range(11)] for i in range(5) for cycle in cycles]
    minors=[(-1)**j*det([row[:j]+row[j+1:] for row in rows]) for j in range(11)]
    assert edge_forms(graph,y)==minors


def check():
    expected=[F(-7,270),F(-7,270),F(-7,810),F(-7,1620)]
    data=json.loads((HERE/'open_theta.json').read_text())
    for shape,target in zip(SHAPES,expected):
        check_heat(shape)
        num,factor=numerators(shape)
        sym=2*prod(factorial(shape.count(i)) for i in set(shape))
        for n in (8,9,10):
            color=literal_colors(n,shape)
            value=factor*sum(map(int,(num*SIGNS*color).tolist()))/(2*sym)
            assert max(map(abs,num.tolist()))*max(map(abs,color.tolist()))<2**62
            assert value==target,(n,shape,value,target)
        assert stored_witness(data[''.join(map(str,shape))])==target
        print(shape,target,flush=True)
    assert sum(expected)==stored_witness(data['sum'])==F(-7,108)
    ce_witness(data['sum'])


if __name__=='__main__':
    check()
