"""Integration by parts of the convergent zero-time faces of the four-spoke wheel."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
from math import factorial, prod

from joint_wheels import CORE, DUAL, ORBIT
from cubature import checked_cut, k4orbit, reduced_powers, reduction
import networkx as nx


@cache
def representative(ex,q=6):
    return min([tuple(ex[i] for i in s) for s in ORBIT]+
               [tuple(q-2-ex[DUAL[i]] for i in s) for s in ORBIT])


@cache
def tree_terms():
    result=[]
    for removed in combinations(range(8),4):
        tree=nx.Graph(e for j,e in enumerate(CORE) if j not in removed)
        if len(tree)==5 and nx.is_tree(tree):
            result.append(tuple(int(i in removed) for i in range(8)))
    return result


def face_data(index,edge,q=6):
    a,b=CORE[edge]
    rest=[i for i in range(8) if i!=edge]
    contracted=[tuple(a if x==b else x for x in CORE[i]) for i in rest]
    plan,name,keep=reduction(contracted)
    zeros=[term for term in tree_terms() if term[edge]==0]
    derivative=[term for term in tree_terms() if term[edge]==1]
    bounds=[]
    for n in range(1,7):
        for subset in combinations(range(7),n):
            degree=min(sum(term[rest[j]] for j in subset) for term in zeros)
            if (q-1)*degree>=n:
                bounds.append((subset,(q-1)*degree-n))
    for powers in product(range(q-2),repeat=7):
        if sum(powers)!=4*(q-1)-7 or any(sum(powers[j] for j in s)<=k for s,k in bounds):
            continue
        beta=dict(zip(rest,powers))
        row=Counter()
        for c in derivative:
            ex=tuple(0 if i==edge else beta[i]+c[i] for i in range(8))
            row[index[representative(ex,q)]]+=q-1
        core,w=reduced_powers(powers,q-1,plan,keep)
        if name=='edge':
            assert core==(0,)
            value=(w,F())
        else:
            assert name=='kite'
            w*=F(prod(factorial(k-1) for k in core),factorial(q-2))
            value=tuple(w*c for c in checked_cut(k4orbit(core)))
        yield (edge,powers),row,value


@cache
def star_plans():
    plans=[]
    for vertices in combinations(range(5),3):
        opposite=[tuple(sorted(set(vertices)-{v})) for v in vertices]
        if not all(e in CORE for e in opposite):
            continue
        ids=tuple(CORE.index(e) for e in opposite)
        star=list(CORE)
        for j,v in zip(ids,vertices):
            star[j]=(v,5)
        plans.append((ids,reduction(star)))
    assert len(plans)==4
    return plans


def star_value(ex,q=6):
    """Polynomial triangle-to-star reduction in the resistance variables."""
    from collections import defaultdict
    results=[]
    for ids,(moves,name,keep) in star_plans():
        surplus=sum(ex[i]+1 for i in ids)-2*q
        if surplus<0:
            continue
        start=list(ex)
        for i in ids:
            start[i]=q-2-ex[i]
        terms={tuple(start):1}
        for _ in range(surplus):
            new=defaultdict(int)
            for term,c in terms.items():
                for i,j in combinations(ids,2):
                    vec=list(term)
                    vec[i]+=1;vec[j]+=1
                    new[tuple(vec)]+=c
            terms=new
        value=[F(),F()]
        for powers,c in terms.items():
            core,w=reduced_powers(powers,q,moves,keep)
            if name=='edge':
                assert core==(0,)
                value[0]+=c*w
            else:
                assert name=='kite'
                w*=F(prod(factorial(k-1) for k in core),factorial(q-1))
                for j,x in enumerate(checked_cut(k4orbit(core))):
                    value[j]+=c*w*x
        results.append(tuple(value))
    if results:
        assert len(set(results))==1
        return results[0]
