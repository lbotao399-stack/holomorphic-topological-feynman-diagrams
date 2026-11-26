"""Universal four-loop C2 wheel monomial table in the two subdivided-wheel graphs."""
from functools import cache
from itertools import product
import json
from pathlib import Path
import networkx as nx
import sympy as sp

from periods import BATCH, CUTS, PERMS, edge_order, form, read_vectors, tagged
from shift import period_evaluator
from four_loop import WHEEL

GRAPH=nx.Graph(WHEEL)
INDEX={e:i for i,e in enumerate(WHEEL)}
WHEEL_PERMS=[tuple(INDEX[tuple(sorted((p[u],p[v])))] for u,v in WHEEL)
             for p in nx.isomorphism.GraphMatcher(GRAPH,GRAPH).isomorphisms_iter()]
assert len(WHEEL_PERMS)==8


def orbit(powers):
    return min(tuple(powers[j] for j in p) for p in WHEEL_PERMS)


@cache
def support():
    bounds=[]
    for mask in range(1,255):
        edges=[j for j in range(8) if mask & (1<<j)]
        g=nx.Graph(WHEEL[j] for j in edges)
        loops=len(edges)-len(g)+nx.number_connected_components(g)
        if loops:
            bounds.append((edges,6*loops))
    seen=set()
    keep=[]
    for ex in product(range(1,6),repeat=8):
        if sum(ex)!=24 or ex.count(5)>1 or ex in seen:
            continue
        seen.update(tuple(ex[j] for j in p) for p in WHEEL_PERMS)
        if all(sum(ex[j] for j in edges)>degree for edges,degree in bounds):
            keep.append(ex)
    assert len(keep)==1421
    return keep


def generate(folder,executable,save):
    integral=period_evaluator(folder,executable,4)
    results=[]
    z3,z5=sp.symbols('z3 z5')
    for start in range(0,len(support()),BATCH):
        ids=range(start,min(start+BATCH,len(support())))
        pairs=[]
        for cut in CUTS:
            terms=[]
            for i in ids:
                powers=support()[i]
                lifted=powers[:6]+(0,)+powers[6:]
                # The contracted rung has zero power and must stay off the cut.
                patterns=(tuple(lifted[p[j]] for j in edge_order(cut)) for p in PERMS)
                rep=max(t for t in patterns if t[0])
                graph='*'.join(f'1/(p{j}.p{j})^{a}' for j,a in enumerate(rep) if j and a)
                terms.append((i,graph))
            vectors=read_vectors(form(folder,executable,f'wheel{start}_{cut}',tagged(terms),cut))
            assert set(vectors)==set(ids)
            pairs.append([integral(vectors[i]) for i in ids])
        assert pairs[0]==pairs[1]
        for i,raw in zip(ids,pairs[0]):
            c=[raw.subs({z3:0,z5:0}),raw.coeff(z3),raw.coeff(z5)]
            assert all(v.is_Rational for v in c) and raw==c[0]+c[1]*z3+c[2]*z5
            # Store the unit-momentum propagator, since current_terms has
            # already included the original simplex and cut factorials.
            results.append([support()[i],list(map(str,c))])
        print('wheel periods in both cuts',len(results),flush=True)
    save.write_text(json.dumps({'edges':WHEEL,'cut_integrals':results},separators=(',',':'))+'\n')


if __name__=='__main__':
    from replay import run
    def verify(folder,executable):
        output=folder/'wheel.json'
        generate(folder,executable,output)
        assert output.read_bytes()==Path(__file__).with_name('wheel.json').read_bytes()
    run(verify)
