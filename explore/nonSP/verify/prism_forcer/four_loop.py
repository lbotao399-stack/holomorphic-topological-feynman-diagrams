"""Reduce a four-loop holomorphic current moment to its primitive graph core."""
from collections import defaultdict
from fractions import Fraction as F
from math import factorial,lcm,prod
from pathlib import Path
import sys

import networkx as nx
import sympy as sp
from flint import fmpz_mpoly_ctx

from periods import PRISM, PERMS, edge_order, form, read_vectors, tagged
from shift import K33,period_evaluator,shifted_polynomial
from check import AUTOMORPHISMS
from vk import G

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'k4_periods'))
from kernel import numerator,tree_flows
from core import EDGES as K4,checked_cut,orbit as k4orbit

# Contract the prism rung 25. The other eight propagators realize the wheel W4.
WHEEL = [tuple(sorted(2 if v==5 else v for v in edge)) for edge in PRISM]
WHEEL.pop(6)


def reduction(edges):
    working = {j:tuple(sorted(e)) for j,e in enumerate(edges)}
    moves = []
    while len(working)>1:
        groups = defaultdict(list)
        for e,pair in working.items():
            groups[pair].append(e)
        parallel = next((ids for ids in groups.values() if len(ids)>1),None)
        if parallel:
            keep,drop = parallel[:2]
            moves.append(('parallel',keep,drop))
            del working[drop]
            continue
        neighbors = defaultdict(list)
        for e,pair in working.items():
            for v in pair:
                neighbors[v].append(e)
        pair = next((ids for ids in neighbors.values() if len(ids)==2),None)
        if not pair:
            break
        keep,drop = pair
        contracted = set(working[keep])^set(working[drop])
        assert len(contracted)==2
        working[keep] = tuple(sorted(contracted))
        del working[drop]
        moves.append(('series',keep,drop))
    if len(working)==1:
        return moves,'edge',list(working)
    graph = nx.Graph(working.values())
    for name,target in [('kite',K4),('wheel',WHEEL),('prism',PRISM),('K33',K33)]:
        iso=nx.isomorphism.GraphMatcher(graph,nx.Graph(target))
        if iso.is_isomorphic():
            order={tuple(sorted((iso.mapping[u],iso.mapping[v]))):i for i,(u,v) in working.items()}
            return moves,name,[order[tuple(sorted(e))] for e in target]
    raise ValueError('No primitive core for this graph.')


def reduced_powers(ex,q,moves,keep):
    powers=[k+1 for k in ex]
    norm=F(1)
    for kind,i,j in moves:
        a,b=powers[i],powers[j]
        if kind=='parallel':
            powers[i]=a+b-q
            a,b=q-a,q-b
        else:
            powers[i]=a+b
        assert a>0 and b>0
        norm *= F(factorial(a-1)*factorial(b-1),factorial(a+b-1))
    return tuple(powers[i] for i in keep),norm


def current_terms(edges,source,insertions):
    graph=G(max(max(e) for e in edges)+1,edges)
    h,remainder=divmod(graph.N-1,4)
    if graph.n!=4 or remainder or len(source)!=graph.nv or any(len(v)!=h for v in source):
        raise ValueError('Expected a four-loop holomorphic graph with a momentum at each vertex.')
    source=[list(map(F,row)) for row in source]
    if any(map(sum,zip(*source))):
        raise ValueError('Sources must sum to zero.')
    den=lcm(*(x.denominator for v in source for x in v))
    integer=[[int(x*den) for x in v] for v in source]
    trees=tree_flows(graph,integer)
    raw=numerator(graph,integer,trees)
    if insertions:
        ctx=fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(graph.N)))
        poly=ctx.from_dict(raw)
        for e,j in insertions:
            if not 0<=e<graph.N or not 0<=j<h:
                raise ValueError('Unknown edge or current component.')
            poly*=ctx.from_dict({tuple(int(i in miss) for i in range(graph.N)):v[e][j]
                                 for miss,v in trees.items() if v[e][j]})
        raw={tuple(map(int,ex)):int(c) for ex,c in poly.terms()}
    q=h+4+len(insertions)
    moves,name,keep=reduction(edges)
    weights=defaultdict(F)
    for ex,c in raw.items():
        powers,norm=reduced_powers(ex,q,moves,keep)
        if name=='edge':
            assert powers==(0,)
        else:
            assert sum(powers)==(3 if name=='kite' else 4)*q and min(powers)>0 and max(powers)<q
            norm*=F(prod(factorial(a-1) for a in powers),factorial(q-1))
        weights[powers]+=c*norm
    scale=den**(4*h+len(insertions))
    return name,q,{p:c/scale for p,c in weights.items() if c}


def two_cuts(folder,executable,name,q,weights):
    if name=='K33':
        edges,perms=K33,AUTOMORPHISMS
    else:
        edges,perms=PRISM,PERMS
    if name=='wheel':
        # Zero propagator power on the contracted rung; cut neither at zero
        # nor on an image of this rung under a prism automorphism.
        weights={p[:6]+(0,)+p[6:]:c for p,c in weights.items()}
    vectors=[]
    for strategy in (0,1):
        cut = 0 if name=='K33' else strategy
        index=edge_order(cut)
        parts=defaultdict(F)
        for powers,c in weights.items():
            orbit={tuple(powers[p[j]] for j in index) for p in perms}
            valid=[v for v in orbit if v[0]]
            levels=sorted({v[0] for v in valid},reverse=True)
            # For a transitive graph choose the two largest distinct exponents.
            level=levels[min(strategy,len(levels)-1)] if name=='K33' else levels[0]
            chosen=max(v for v in valid if v[0]==level)
            parts[chosen[1:]]+=c
        terms=[]
        for powers,c in parts.items():
            if c:
                prop='*'.join(f'1/(p{j}.p{j})^{a}' for j,a in enumerate(powers,1) if a)
                terms.append(f'+({c})*'+prop)
        if terms:
            text=form(folder,executable,f'{name}{strategy}',tagged([(0,''.join(terms))]),cut,edges=edges)
            vectors.append(read_vectors(text).get(0,{}) if text.strip()!='0' else {})
        else:
            vectors.append({})
    if not any(vectors):
        return sp.Integer(0)
    shifts=q-2
    def valuation(p):
        return next(i for i,c in enumerate(shifted_polynomial(p,shifts)) if c)
    last=max([0]+[valuation(b)-valuation(a) for row in vectors for a,b in row.values()])
    get=period_evaluator(folder,executable,last,shifts,edges)
    pair=[get(row) for row in vectors]
    assert pair[0]==pair[1]
    return pair[0]


def evaluate(folder,executable,edges,source,insertions=()):
    name,q,weights=current_terms(edges,source,insertions)
    if name=='edge':
        return sp.Rational(weights.get((0,),0))
    if name=='kite':
        total=[F(0),F(0)]
        for powers,c in weights.items():
            for i,v in enumerate(checked_cut(k4orbit(powers))):
                total[i]+=c*v
        return total[0]+total[1]*sp.Symbol('z3')
    return two_cuts(folder,executable,name,q,weights)
