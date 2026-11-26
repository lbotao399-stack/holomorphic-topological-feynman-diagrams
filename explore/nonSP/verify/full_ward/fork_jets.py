"""First nonzero displacement order of k three-loop forks on one marked edge."""
from fractions import Fraction as F
from math import factorial

from forest_jets import edge, evaluate, sign


def check(k):
    a,b=3*k,3*k+1
    block_edges=[]
    chosen=[]
    spokes=[]
    source=[(0,0)]*(3*k)+[(1,0),(-1,-k)]
    for j in range(k):
        u,p,q=3*j,3*j+1,3*j+2
        source[u]=(0,1)
        block_edges.extend(((p,u),(p,a),(p,b),(q,u),(q,a),(q,b)))
        chosen.append((p,u))
        spokes.append((p,a))
    block_edges.append((a,b))
    edges=tuple(sorted(edge(*e) for e in block_edges))
    order=sign([edges.index(edge(*e)) for e in block_edges])
    insertions=[(edges.index(edge(*e)),0) for e in chosen]
    n=3*k
    paper_sign=(-1)**(n*(n+1)//2)*order
    assert evaluate(edges,a,b,source,insertions[:-1])==0
    value=paper_sign*(-1)**sum(i>j for i,j in chosen)*evaluate(edges,a,b,source,insertions)
    assert value==F((-1)**k,36**k*factorial(4*k+1))
    physical=paper_sign*(-1)**sum(i>j for i,j in spokes)*evaluate(
        edges,a,b,source,[(edges.index(edge(*e)),0) for e in spokes])
    assert physical==F(1,72**k*factorial(4*k+1))
    print('loops',n,'first current degree',k,'internal and A-spoke moments',value,physical,flush=True)


if __name__=='__main__':
    for k in range(1,5):
        check(k)
