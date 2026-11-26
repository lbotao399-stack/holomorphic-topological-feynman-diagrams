"""Exact signed loop-raising identity from adding a parallel H-edge path."""
from fractions import Fraction as F
from math import comb,factorial
from pathlib import Path
from random import Random
import sys

VERIFY=Path(__file__).resolve().parents[2]/'nonSP/verify'
sys.path.insert(0,str(VERIFY))
from vk import G,currents_exact,density_exact,fdet


def graft(graph,e,h):
    u,v=graph.edges[e]
    def old(i):
        return i+h-1 if i==graph.root else i
    path=[old(u),*range(graph.root,graph.root+h-1),old(v)]
    return G(graph.nv+h-1,[(old(a),old(b)) for a,b in graph.edges]+list(zip(path,path[1:])))


def child_sources(source,directions,s):
    shifted=[a-s*sum(p[j] for p in directions) for j,a in enumerate(source[-1])]
    return source[:-1]+[[s*x for x in p] for p in directions]+[shifted]


def ratio(graph,e,source,directions,t,rho,sigma):
    h=len(source[0])
    u,v=graph.edges[e]
    child=graft(graph,e,h)
    x=[a-b for a,b in zip([F(1)]+sigma,sigma+[F(0)])]
    child_t=t[:]
    child_t[e]=t[e]/rho
    child_t.extend(t[e]*v/(1-rho) for v in x)
    full=child_sources(source,directions,1)
    reduced=[row[:] for row in source[:-1]]+[full[-1][:]]
    for p,s in zip(directions,sigma):
        for j,value in enumerate(p):
            reduced[u][j]+=s*value
            reduced[v][j]+=(1-s)*value
    current=currents_exact(graph,t,reduced[:-1])[e]
    parent=density_exact(graph,t,reduced[:-1])
    local=rho**(h-1)*fdet([current]+directions)
    J=t[e]**h/(rho**2*(1-rho)**(h+1))
    actual=density_exact(child,child_t,full[:-1])*J
    if not parent*local:
        assert not actual
        return 0
    sign=actual/(local*parent)
    assert sign in (-1,1)
    return sign


def sign_at_random_points(graph,e,h,seed):
    rng=Random(seed)
    signs=[]
    while len(signs)<3:
        lam=[[rng.randrange(-4,5) for _ in range(h)] for v in range(graph.root)]
        lam.append([-sum(row[j] for row in lam) for j in range(h)])
        p=[[rng.randrange(-3,4) for _ in range(h)] for _ in range(h-1)]
        t=[F(rng.randrange(1,8),rng.randrange(1,6)) for edge in graph.edges]
        rho=F(rng.randrange(1,11),11)
        sigma=sorted([F(i,19) for i in rng.sample(range(1,19),h-1)],reverse=True)
        value=ratio(graph,e,lam,p,t,rho,sigma)
        if value:signs.append(value)
    assert len(set(signs))==1
    return signs[0]


def vertex_test(h):
    sys.path.insert(0,str(VERIFY.parents[2]/'code'))
    from htk.graph import Graph
    from htk.spred import sp_reduce
    import sympy as sp
    edges=[(i,i+1) for i in range(h)]+[(0,h)]
    graph=G(h+1,edges)
    source=[[int(i==j) for j in range(h)] for i in range(h)]+[[-1]*h]
    directions=source[1:h]
    parent=sp_reduce(Graph(h+1,edges),h).integrate(source[:-1])
    value=0
    t=sp.Symbol('s')
    for e,pair in enumerate(edges):
        if 0 not in pair:
            continue
        child=graft(graph,e,h)
        reduce=sp_reduce(Graph(child.nv,child.edges),h)
        samples=[(s,reduce.integrate(child_sources(source,directions,s)[:-1])) for s in range(2*h+1)]
        coefficient=sp.Poly(sp.interpolate(samples,t),t).nth(h-1)
        value+=((1 if pair[0]==0 else -1)*sign_at_random_points(graph,e,h,477+e)*coefficient)
    assert factorial(h)*value==parent
    print('integrated vertex identity H',h,'parent',parent,flush=True)


def integrated_test():
    sys.path.insert(0,str(VERIFY/'prism_forcer'))
    from four_loop import evaluate
    import sympy as sp
    edges=[(0,1),(0,2),(0,3),(1,4),(4,2),(1,3),(2,3)]
    graph=G(5,edges)
    source=[[1,2],[2,-1],[-1,3],[-3,-4],[1,0]]
    directions=[[0,1]]
    child=graft(graph,0,2)
    vals=[evaluate(None,None,child.edges,child_sources(source,directions,s)) for s in range(9)]
    assert vals[0]==0
    first=sp.expand(sum(F((-1)**(s+1)*comb(8,s),s)*vals[s] for s in range(1,9)))
    sign=sign_at_random_points(graph,0,2,41429)
    expected=F(222317,25200)-F(3021,700)*sp.Symbol('z3')
    assert sp.expand(first-sign*expected/2)==0
    print('signed first coefficient on the 01 graft',first,'sign',sign,flush=True)
    # The remaining insertion on the direct edge brings one further factor rho.
    vals=[evaluate(None,None,child.edges,child_sources(source,directions,s),[(0,0)]) for s in range(10)]
    assert vals[0]==0
    first=sp.expand(sum(F((-1)**(s+1)*comb(9,s),s)*vals[s] for s in range(1,10)))
    expected=-F(4003781,302400)+F(3772,525)*sp.Symbol('z3')
    assert sp.expand(first-sign*expected/3)==0
    print('same graft retaining one 01 insertion',first,flush=True)


if __name__=='__main__':
    sys.path.insert(0,str(VERIFY.parents[0]))
    from graphs import K4SUB,PRISM
    examples=[(G(3,[(0,1),(1,2),(0,2)]),2),(G(6,PRISM),2),
              (G(4,[(0,1),(1,2),(2,3),(0,3)]),3),(G(8,K4SUB[0][1]),3),
              (G(5,[(0,1),(1,2),(2,3),(3,4),(0,4)]),4)]
    for i,(graph,h) in enumerate(examples):
        incident=[e for e,pair in enumerate(graph.edges) if 0 in pair]
        print(h,graph.n,[sign_at_random_points(graph,e,h,3100+i*13+e) for e in incident],flush=True)
    for h in (2,3):
        vertex_test(h)
    integrated_test()
