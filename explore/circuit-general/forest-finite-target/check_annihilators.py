"""Independent exact audit of the new K33+edge target annihilator rows.

No full numerator or period table is generated.  The row expansion below
is independent of syz5.relation; archived low-q periods are read only.
"""
import argparse
import collections
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import itertools
import json
from math import prod
from pathlib import Path
import pickle
import sys
from time import perf_counter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'explore/circuit-general/symmetry'))
import engine
from family_target import FamilyTarget
sys.path.insert(0,str(ROOT/'explore/circuit-general/fiveloop'))
import per5
import syz5

original_table=per5.table
per5.table=lambda name,q,build=False,verbose=False: original_table(name,q,build=False,verbose=False)
certificate=json.loads((ROOT/'explore/circuit-general/forest-target-reduction/dimension_certificate.json').read_text())
E=tuple(map(tuple,certificate['scalar_core']['edges']))
oldE=tuple(per5.SKEL['K33+e'])
mapping=next(engine.nx.isomorphism.GraphMatcher(engine.nx.Graph(oldE),engine.nx.Graph(E)).isomorphisms_iter())
indices={frozenset(edge):i for i,edge in enumerate(E)}
permutation=tuple(indices[frozenset((mapping[u],mapping[v]))] for u,v in oldE)

def transport(monomial):
    result=[0]*len(E)
    for i,j in enumerate(permutation):result[j]=monomial[i]
    return tuple(result)

cache=ROOT/'explore/circuit-general/fiveloop/syz_cache/0-1_0-2_0-3_0-4_1-2_1-3_1-4_2-5_3-5_4-5_k2.pkl'
generators=[]
for degree,rows in pickle.loads(cache.read_bytes()).items():
    for aa,bb in rows:
        a=[{} for _ in E]
        for i,j in enumerate(permutation):a[j]={transport(m):F(c) for m,c in aa[i].items()}
        generators.append((degree,(a,{transport(m):F(c) for m,c in bb.items()})))

def circuit_rank(edges):
    vertices={v for edge in edges for v in edge}
    if not edges:return 0
    return len(edges)-len(vertices)+engine.nx.number_connected_components(engine.nx.Graph(edges))

@lru_cache(None)
def subsets(edges):
    answer=[]
    for mask in range(1,(1<<len(edges))-1):
        selected=tuple(i for i in range(len(edges)) if mask>>i&1)
        loops=circuit_rank(tuple(edges[i] for i in selected))
        if loops:answer.append((selected,loops))
    return tuple(answer)

@lru_cache(None)
def convergent(edges,powers,q):
    if min(powers)<=0 or sum(powers)!=q*circuit_rank(edges):return False
    return all(sum(powers[i] for i in selected)>q*loops for selected,loops in subsets(edges))

def face(edge):
    u,v=E[edge]
    return tuple((u if a==v else a,u if b==v else b) for j,(a,b) in enumerate(E) if j!=edge)

def raw_relation(syzygy,beta,q):
    """Divergence and every coordinate face, with a separate strict guard."""
    a,b=syzygy
    values=collections.defaultdict(F)
    for edge,polynomial in enumerate(a):
        for monomial,coefficient in polynomial.items():
            powers=[x+y for x,y in zip(beta,monomial)]
            derivative=powers[edge]-1
            if derivative:
                powers[edge]-=1
                if not convergent(E,tuple(powers),q):return None
                values[E,tuple(powers)]+=coefficient*derivative
            else:
                boundary=tuple(p for j,p in enumerate(powers) if j!=edge)
                contracted=face(edge)
                if not convergent(contracted,boundary,q):return None
                values[contracted,boundary]+=coefficient
    for monomial,coefficient in b.items():
        powers=tuple(x+y for x,y in zip(beta,monomial))
        values[E,powers]-=q*coefficient
    values={key:coefficient for key,coefficient in values.items() if coefficient}
    if any(not convergent(graph,powers,q) for graph,powers in values):return None
    return values

def collected(items):
    result=collections.defaultdict(F)
    for graph,powers,coefficient in items:result[tuple(graph),tuple(powers)]+=F(coefficient)
    return {key:coefficient for key,coefficient in result.items() if coefficient}

@lru_cache(None)
def value(graph,powers,q):
    result=per5.value(list(graph),list(powers),q)
    return result

def polynomials():
    # A separate spanning-tree enumeration, without syz5.upoly/gutil.cotrees.
    archive=json.loads((HERE/'annihilators.json').read_text())
    assert tuple(map(tuple,archive['edges']))==oldE
    transported=[]
    for row in archive['generators']:
        vector=[{} for _ in E]
        for old,new in enumerate(permutation):
            vector[new]={transport(tuple(powers)):F(coefficient)
                         for powers,coefficient in row['vector_fields'][old]}
        scalar={transport(tuple(powers)):F(coefficient)
                for powers,coefficient in row['scalar_multiplier']}
        transported.append((row['degree'],(vector,scalar)))
    assert transported==generators,'JSON archive differs from the independently loaded original generators'
    U={}
    for tree in itertools.combinations(range(len(E)),len({v for e in E for v in e})-1):
        if engine.nx.is_tree(engine.nx.Graph([E[i] for i in tree])):
            U[tuple(int(i not in tree) for i in range(len(E)))]=F(1)
    assert U==syz5.upoly(E)
    for degree,(a,b) in generators:
        residual=collections.defaultdict(F)
        assert all(sum(m)==degree for polynomial in a for m in polynomial)
        assert all(sum(m)==degree-1 for m in b)
        for i,polynomial in enumerate(a):
            for m,c in polynomial.items():
                for u,uc in U.items():
                    if u[i]:
                        exponent=[x+y for x,y in zip(m,u)];exponent[i]-=1
                        residual[tuple(exponent)]+=c*uc*u[i]
        for m,c in b.items():
            for u,uc in U.items():residual[tuple(x+y for x,y in zip(m,u))]-=c*uc
        assert not {m:c for m,c in residual.items() if c}
    return len(U)

def compositions(total,parts,maximum):
    if parts==0:
        if total==0:yield ()
        return
    for first in range(max(1,total-maximum*(parts-1)),min(maximum,total-parts+1)+1):
        for rest in compositions(total-first,parts-1,maximum):yield (first,)+rest

def level(q,maximum_rows=None):
    counts=collections.Counter();by_generator=collections.Counter();boundary_generators=set()
    classifier=FamilyTarget(E,[],q).system
    classifier.counts['independent_guard']=0
    for identifier,(degree,syzygy) in enumerate(generators):
        for beta in compositions(5*q-degree+1,len(E),q):
            counts['candidates']+=1
            reference=syz5.relation(E,syzygy,list(beta),q)
            if reference is None:
                counts['divergence_guard_rejected']+=1;continue
            row=raw_relation(syzygy,beta,q)
            if row is None:
                counts['strict_face_or_summed_guard_rejected']+=1
                assert not classifier.add_relation(reference,'independent_guard'),(q,identifier,beta,'classifier accepted a divergent face')
                counts['classifier_guard_agrees']+=1
                continue
            assert row==collected(reference),(q,identifier,beta,'row mismatch')
            total=collections.defaultdict(F);unknown=False
            for (graph,powers),coefficient in row.items():
                exact=value(graph,powers,q)
                if exact is None:unknown=True;break
                for label,entry in exact.items():total[label]+=coefficient*entry
            if unknown:
                counts['unavailable_low_period']+=1;continue
            assert not {label:c for label,c in total.items() if c},(q,identifier,beta,dict(total))
            counts['verified']+=1;by_generator[identifier]+=1
            if row:counts['nontrivial']+=1
            if any(len(graph)<len(E) for graph,powers in row):
                counts['with_coordinate_faces']+=1;boundary_generators.add(identifier)
            if maximum_rows is not None and by_generator[identifier]>=maximum_rows:break
    return {'q':q,'counts':dict(counts),'verified_by_generator':dict(by_generator),
            'boundary_generators':sorted(boundary_generators)}

def rational_q7_checks():
    U=syz5.upoly(E)
    def evaluate(poly,point):
        return sum((c*prod(x**p for x,p in zip(point,m)) for m,c in poly.items()),F())
    def derivative(poly,index,point):
        result=F()
        for monomial,coefficient in poly.items():
            if monomial[index]:
                result+=coefficient*monomial[index]*prod(
                    x**(p-int(j==index)) for j,(x,p) in enumerate(zip(point,monomial)))
        return result
    totals=collections.Counter();records=[]
    for identifier,(degree,syzygy) in enumerate(generators):
        a,b=syzygy
        shifts=set(b)
        for i,polynomial in enumerate(a):
            for monomial in polynomial:
                shift=list(monomial);shift[i]-=1;shifts.add(tuple(shift))
        selected=None
        for alpha,unused in certificate['all_residual_scalar_periods']:
            for shift in sorted(shifts):
                beta=tuple(x-y for x,y in zip(alpha,shift))
                if min(beta)>0:
                    row=raw_relation(syzygy,beta,7)
                    if row is not None:
                        selected=beta,row;break
            if selected is not None:break
        assert selected is not None,identifier
        beta,row=selected
        assert sum(beta)==35-degree+1
        for point in (tuple(F(i+2,i+1) for i in range(len(E))),
                      tuple(F(2*i+3,i+2) for i in range(len(E)))):
            u=evaluate(U,point)
            monomial=prod(x**(p-1) for x,p in zip(point,beta))
            direct=sum(derivative(polynomial,i,point)+(beta[i]-1)*evaluate(polynomial,point)/point[i]
                       -7*evaluate(polynomial,point)*derivative(U,i,point)/u
                       for i,polynomial in enumerate(a))*monomial/u**7
            expanded=sum(c*prod(x**(p-1) for x,p in zip(point,powers))/u**7
                         for (graph,powers),c in row.items() if graph==E)
            assert direct==expanded,(identifier,beta,'rational divergence')
            totals['rational_divergence_values']+=1
            for i,polynomial in enumerate(a):
                if beta[i]!=1:continue
                boundary_point=list(point);boundary_point[i]=F()
                retained_point=tuple(x for j,x in enumerate(point) if j!=i)
                contracted_u=evaluate(syz5.upoly(face(i)),retained_point)
                assert contracted_u==evaluate(U,boundary_point)
                boundary_monomial=prod(x**(p-1) for x,p in zip(boundary_point,beta))
                direct_face=evaluate(polynomial,boundary_point)*boundary_monomial/evaluate(U,boundary_point)**7
                expanded_face=sum(c*prod(x**(p-1) for j,(x,p) in enumerate(zip(boundary_point,beta)) if j!=i)
                                  *prod(x**p for j,(x,p) in enumerate(zip(boundary_point,m)) if j!=i)
                                  /contracted_u**7 for m,c in polynomial.items() if not m[i])
                assert direct_face==expanded_face,(identifier,beta,i,'face value')
                totals['rational_face_values']+=1
        records.append({'generator':identifier,'degree':degree,'beta':beta})
    totals['generators']=len(records)
    return {'counts':dict(totals),'records':records}

def target():
    terms=[(tuple(p),F(c)) for p,c in certificate['all_residual_scalar_periods']]
    prepared=FamilyTarget(E,terms,7)
    assert len(prepared.system.reps)==len(terms)==len(prepared.target)==17
    assert all(tuple(prepared.system.reps[i][1])==powers for i,(powers,coefficient) in enumerate(terms))
    assert prepared.target=={i:c for i,(_,c) in enumerate(terms)}
    result={label:F(coefficient) for label,coefficient in certificate['finite_coefficients'].items()}
    extra={'1':F(529509659453,3990767616000),'z3':F(272322049,769824000),'z5':-F(57527,106920)}
    for label,coefficient in extra.items():result[label]=result.get(label,F())+coefficient
    result={label:c for label,c in result.items() if c}
    assert result=={'1':F(1187,435456000),'z3':-F(816,435456000)}
    return {label:str(c) for label,c in result.items()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='refresh the complete independent audit record')
    arguments=parser.parse_args()
    start=perf_counter();out={'generators':len(generators),'U_terms':polynomials(),'target_sum':target()}
    print(out,flush=True)
    out['q3']=level(3);print('q3',out['q3'],flush=True)
    out['q4']=level(4,maximum_rows=50);print('q4',out['q4'],flush=True)
    out['q7']=rational_q7_checks();print('q7',out['q7']['counts'],flush=True)
    assert not out['q3']['counts'].get('unavailable_low_period')
    assert not out['q4']['counts'].get('unavailable_low_period')
    assert len(out['q4']['verified_by_generator'])==len(generators)
    assert len(out['q4']['boundary_generators'])==len(generators)
    sources=[Path(__file__).resolve(),cache,HERE/'annihilators.json',
             ROOT/'explore/circuit-general/forest-target-reduction/dimension_certificate.json',
             ROOT/'explore/circuit-general/symmetry/family_target.py',
             ROOT/'explore/circuit-general/cores/gutil.py']
    sources += [ROOT/'explore/circuit-general/fiveloop'/name for name in
                ('syz5.py','per5.py','msys5.py','tables/K33+e_q3.json','tables/K33+e_q4.json',
                 'tables/prism_q3.json','tables/prism_q4.json','tables/K33_q3.json','tables/K33_q4.json')]
    out['source_sha256']={str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest()
                          for path in sources}
    out['scope']='Independent polynomial, divergence, boundary, low-level period, target-order and final-arithmetic controls; the complete target witness is checked separately.'
    out['edges']=E
    # Normalize JSON tuple/key types before a deterministic whole-record replay.
    out=json.loads(json.dumps(out))
    destination=HERE/'annihilator_audit.json'
    if arguments.write:
        destination.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
        print('Wrote the independent annihilator audit record.',flush=True)
    else:
        assert out==json.loads(destination.read_text()),'independent audit record changed'
        print('Every independent annihilator audit field replays.',flush=True)
    print('All independent annihilator checks passed (%.3f seconds).'%(perf_counter()-start),flush=True)


if __name__=='__main__':
    main()
