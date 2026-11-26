"""Independent physical single-trace probe with eight SO(32) generators on an eight-cycle."""
from itertools import combinations,permutations
from fractions import Fraction as F
from math import factorial,prod
import json
from check_open_theta import P,SIGNS,evaluate_terms,literal_colors,numerators
from open_theta import HERE,SHAPES,det,parity


def trace(edges):
    """T_ij has nonzero matrix entries +1 at (i,j), -1 at (j,i)."""
    first=edges[0]
    result=0
    for start in first:
        current=start
        value=1
        for i,j in edges:
            if current==i:
                current=j
            elif current==j:
                current=i
                value=-value
            else:
                break
        else:
            if current==start:
                result+=value
    return result


def bracket(a,b):
    """Nonzero commutator of oriented edge matrices is another edge matrix; orientation includes the sign."""
    if a is None:return None
    i,j=a;k,l=b
    if {i,j}=={k,l}:return None
    if j==k:return (i,l)
    if i==l:return (j,k)
    if j==l:return (k,i)
    if i==k:return (l,j)
    return None


def alternating_mixed(edges):
    cache_minors={s:det([P[i] for i in s]) for s in combinations(range(8),5)}
    def minor(s):
        return parity(s)*cache_minors[tuple(sorted(s))]
    value=0
    for s in permutations(range(8)):
        x=edges[s[0]]
        for j in s[1:5]:
            x=bracket(edges[j],x) if x is not None else None
        if x is None:continue
        color=trace([x,*[edges[i] for i in s[5:]]])
        if color:
            second=sum(w*minor((s[0],*s[5:],s[k])) for k,w in ((1,3),(2,2),(3,1)))
            value+=parity(s)*color*minor(s[:5])*second
    return value


def seven_cochain(edges,p,terms):
    def evaluation(s,first,second):
        return det([p[s[i]] for i in first])*det([p[s[i]] for i in second])
    total=0
    for s in permutations(range(7)):
        color=trace([edges[i] for i in s])
        if color:
            total+=color*parity(s)*sum(F(c)*evaluation(s,a,b) for (_,a,b),c in terms)
    return total


def main():
    edges=[(i,(i+1)%8) for i in range(8)]
    source=json.loads((HERE/'single_trace_preimages.json').read_text())
    exact=0
    for i,j in combinations(range(8),2):
        joined=bracket(edges[i],edges[j])
        if joined is None:continue
        rest=[k for k in range(8) if k not in (i,j)]
        p=[[x+y for x,y in zip(P[i],P[j])]]+[P[k] for k in rest]
        a=[joined]+[edges[k] for k in rest]
        exact+=(-1)**(i+j)*seven_cochain(a,p,source['pure_gauge_preimages'])
    mixed=alternating_mixed(edges)
    graphs=json.loads((HERE/'open_theta.json').read_text())['sum']
    pure=sum(32**d['rank']*evaluate_terms(d,P) for d in graphs if d['lengths']==[8])
    direct=0
    for shape in SHAPES:
        num,factor=numerators(shape)
        sym=2*prod(factorial(shape.count(j)) for j in set(shape))
        for n,lagrange in ((8,253),(9,-528),(10,276)):
            colors=literal_colors(n,shape,edges)
            direct+=factor*lagrange*sum(map(int,(num*SIGNS*colors).tolist()))/(2*sym)
    assert direct==pure
    assert (pure,mixed,exact)==(F(-299,17280),-364,F(143,5760))
    assert pure-F(1,8640)*mixed==exact
    print('SO(32) eight-cycle currents: pure theta',pure,'raw mixed',mixed,'gauge differential',exact)


if __name__=='__main__':
    main()
