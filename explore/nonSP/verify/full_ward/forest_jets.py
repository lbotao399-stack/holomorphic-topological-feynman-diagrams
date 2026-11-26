"""Exact current insertions on cubic marked-product forests with root sources."""
from functools import cache
from fractions import Fraction as F
from math import factorial, prod
from pathlib import Path
import sys
from typing import NamedTuple

import networkx as nx
from flint import fmpq, fmpq_mpoly_ctx

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from four_mark_relations import sign


class Current(NamedTuple):
    x: object
    y: object

    def __add__(self,other):
        return Current(self.x+other.x,self.y+other.y)

    def __sub__(self,other):
        return Current(self.x-other.x,self.y-other.y)

    def __neg__(self):
        return Current(-self.x,-self.y)

    def __mul__(self,w):
        return Current(self.x*w,self.y*w)

    def bracket(self,other):
        return self.x*other.y-self.y*other.x


def edge(i,j):
    return (i,j) if i<j else (j,i)


def get(currents,i,j):
    return currents[i,j] if i<j else -currents[j,i]


def put(currents,i,j,current):
    currents[edge(i,j)]=current if i<j else -current


@cache
def moves(edges,a,b):
    graph=nx.Graph(edges)
    if len(graph)==2:
        assert set(edges)=={edge(a,b)}
        return []
    leaves={}
    for v in sorted(set(graph)-{a,b}):
        neighbors=set(graph[v])
        if neighbors=={a,b}:
            kind,vertices='series',(v,a,b)
            break
        if {a,b}<neighbors and len(neighbors)==3:
            u,=neighbors-{a,b}
            leaves.setdefault(u,[]).append(v)
    else:
        for u,group in leaves.items():
            spokes=set(graph[u])&{a,b}
            if len(spokes)==1:
                x,=spokes
                w,={a,b}-{x}
                kind,vertices='single',(group[0],u,w,x)
                break
            if len(group)==2 and not spokes:
                kind,vertices='paired',(*group,u,a,b)
                break
        else:
            raise ValueError('The graph admits no rooted cubic-forest reduction.')
    if kind=='series':
        v,u,w=vertices
        before=((v,u),(v,w))
        after=()
    elif kind=='single':
        v,u,w,x=vertices
        before=((v,u),(v,w),(v,x),(u,x),(w,x))
        after=((u,w),(u,x),(w,x))
    else:
        p,q,u,a,b=vertices
        before=((p,u),(p,a),(p,b),(q,u),(q,a),(q,b),(a,b))
        after=((u,a),(u,b),(a,b))
    before=[edge(*e) for e in before]
    after=[edge(*e) for e in after]
    rest=[e for e in edges if e not in before]
    orient=sign([edges.index(e) for e in before+rest])*sign(after+rest)
    reduced=tuple(sorted(after+rest))
    return [(kind,vertices,orient)]+moves(reduced,a,b)


def chart(edges,a,b,source,factored=False,context=None):
    """The returned density uses the catalog K order; source rows match vertex labels."""
    steps=moves(tuple(edges),a,b)
    loops=len(edges)-len(source)+1
    ctx=context if context is not None else fmpq_mpoly_ctx.get(tuple('r'+str(j) for j in range(2*loops)))
    def scalar(v):
        return v if hasattr(v,'context') and v.context()==ctx else ctx.constant(fmpq(str(v)))
    momenta={i:Current(*(scalar(v) for v in row)) for i,row in enumerate(source)}
    if any(sum((p[j] for p in momenta.values()),ctx.constant(0)) for j in (0,1)):
        raise ValueError('The source momenta must sum to zero.')
    variables=ctx.gens()
    replay=[]
    simplices=[]
    start=0
    phase=(-1)**(loops*(loops+1)//2)
    for kind,vertices,orient in steps:
        phase*=orient
        n=4 if kind=='paired' else 2
        r=variables[start:start+n]
        if kind=='series':
            v,u,w=vertices
            local=momenta.pop(v)
            momenta[u]+=local*r[0]
            momenta[w]+=local*(1-r[0])
        else:
            local=None
            for v in vertices[:2 if kind=='paired' else 1]:
                p=momenta.pop(v)
                if p.x or p.y:
                    raise ValueError('Star reduction requires zero source at the eliminated vertex.')
        replay.append((kind,vertices,r,local))
        simplices+=([(start,),(start+1,),(start+2,start+3)] if n==4 else [(start,),(start+1,)])
        start+=n
    assert start==2*loops
    current={}
    put(current,a,b,momenta[a])
    factors=[ctx.constant(phase)]
    for kind,vertices,r,local in reversed(replay):
        if kind=='series':
            v,u,w=vertices
            y=get(current,u,w)
            factors.extend((r[1],local.bracket(y)))
            first=local*r[0]-y*(1-r[1])
            put(current,v,u,first)
            put(current,v,w,local-first)
            put(current,u,w,y*r[1])
        elif kind=='single':
            v,u,w,x=vertices
            y=get(current,u,w)
            l,h=get(current,u,x),get(current,w,x)
            del current[edge(u,w)]
            factors.extend((r[0],r[1],l.bracket(h)))
            put(current,v,u,-y-l*(1-r[0]))
            put(current,v,w,y-h*(1-r[1]))
            put(current,v,x,l*(1-r[0])+h*(1-r[1]))
            put(current,u,x,l*r[0])
            put(current,w,x,h*r[1])
        else:
            p,q,u,a,b=vertices
            l,h,y=get(current,u,a),get(current,u,b),get(current,a,b)
            del current[edge(u,a)],current[edge(u,b)]
            factors.extend((r[0]-r[1],1-r[2]-r[3],l.bracket(y),h.bracket(y)))
            put(current,p,u,-l*r[0]-h*r[1])
            put(current,q,u,-l*(1-r[0])-h*(1-r[1]))
            put(current,p,a,l*r[0]-y*r[2])
            put(current,q,a,l*(1-r[0])-y*r[3])
            put(current,p,b,h*r[1]+y*r[2])
            put(current,q,b,h*(1-r[1])+y*r[3])
            put(current,a,b,y*(1-r[2]-r[3]))
    assert set(current)==set(edges)
    for v,row in enumerate(source):
        for j in (0,1):
            assert sum((get(current,v,t[1] if t[0]==v else t[0])[j] for t in edges if v in t),
                       ctx.constant(0))==scalar(row[j])
    return factors if factored else prod(factors),current,tuple(simplices)


@cache
def volume(exponents,simplices):
    return F(prod(factorial(k) for k in exponents),
             prod(factorial(len(s)+sum(exponents[j] for j in s)) for s in simplices))


def integrate(poly,simplices):
    return sum((F(int(c.p),int(c.q))*volume(tuple(map(int,exp)),simplices)
                for exp,c in poly.terms()),F())


def evaluate(edges,a,b,source,insertions=()):
    factors,current,simplices=chart(edges,a,b,source,factored=True)
    for e,j in insertions:
        factors.append(get(current,*edges[e])[j])
    return F(str(integrate_factors(factors,simplices)))


def integrate_factors(factors,simplices):
    """Contract each ratio simplex with its dependent factors before expanding others."""
    if not all(factors):
        return factors[0].context().constant(0)
    active=[(p,set(i for i,v in enumerate(p.degrees()) if v)) for p in factors]
    remaining=list(simplices)
    while remaining:
        group=min(remaining,key=lambda s:prod(len(p) for p,idx in active if idx.intersection(s)))
        remaining.remove(group)
        terms=[p for p,idx in active if idx.intersection(group)]
        if not terms:
            active.append((factors[0].context().constant(fmpq(1,factorial(len(group)))),set()))
            continue
        active=[item for item in active if not item[1].intersection(group)]
        polynomial=prod(terms)
        coeff={}
        for exp,c in polynomial.terms():
            powers=[int(exp[j]) for j in group]
            scale=fmpq(prod(factorial(k) for k in powers),factorial(len(group)+sum(powers)))
            key=tuple(0 if j in group else int(v) for j,v in enumerate(exp))
            coeff[key]=coeff.get(key,0)+c*scale
        reduced=polynomial.context().from_dict(coeff)
        if not reduced:
            return polynomial.context().constant(0)
        active.append((reduced,set(i for i,v in enumerate(reduced.degrees()) if v)))
    return prod(p for p,_ in active)
