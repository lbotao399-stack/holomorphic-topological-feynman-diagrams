"""Projective integration by parts on orbit-summed graph numerators."""
from collections import Counter
from argparse import ArgumentParser
from fractions import Fraction as F
from gzip import open as open_gzip
from hashlib import sha256
from itertools import product
from math import factorial, gcd, lcm, prod
from pathlib import Path
from random import Random
import json
import numpy as np

from flint import nmod_mat

from quartic_planar_periods import checked_cut, core_data, orbit, reduced_powers, reduction
from four_loop import AUTOMORPHISMS, K33, PERMS, PRISM


def fold(row,permutations):
    return min(tuple(row[j] for j in p) for p in permutations)


def primitive(row):
    content = gcd(*row.values())
    sign = 1 if row[min(row)]>0 else -1
    return {j:v//(sign*content) for j,v in row.items() if v}


def identities(monomial,cotrees,children,q):
    """Eliminate the derivative at fixed dimension using U times the lower form."""
    total = Counter(children)
    for i,k in enumerate(monomial):
        containing = Counter(c for tree,c in zip(cotrees,children) if tree[i])
        row = {c:(q-1)*containing[c]-(k+1)*v for c,v in total.items()}
        row = {c:v for c,v in row.items() if v}
        if row:
            yield primitive(row)


def triangular(rows):
    # Sparse exact row differences with the same leading monomial give
    # additional pivots without ever constructing a dense rational matrix.
    basis = {}
    for initial in rows:
        row = initial
        while row:
            lead = min(row)
            previous = basis.get(lead)
            if previous is None:
                basis[lead]=row
                break
            if len(row)<len(previous):
                row,basis[lead]=previous,row
                previous=basis[lead]
            a,b = row[lead],previous[lead]
            g = gcd(a,b)
            a,b = a//g,b//g
            new = {j:b*v for j,v in row.items()}
            for j,v in previous.items():
                new[j]=new.get(j,0)-a*v
            new = {j:v for j,v in new.items() if v}
            row = primitive(new) if new else {}
    return basis


def collapse(weights,basis):
    result = dict(weights)
    for pivot,row in sorted(basis.items()):
        weight = result.pop(pivot,0)
        if weight:
            scale = F(weight,row[pivot])
            for j,c in row.items():
                if j!=pivot:
                    result[j]=result.get(j,0)-scale*c
    return {j:v for j,v in result.items() if v}


def faces(graph,e):
    a,b=graph[e]
    delete=[edge for i,edge in enumerate(graph) if i!=e]
    contract=[tuple(sorted(a if v==b else v for v in edge)) for edge in delete]
    return delete,contract


def lower_period(moves,name,keep,powers,q):
    propagated,norm=reduced_powers(powers,q,moves,keep)
    if name=='edge':
        assert propagated==(0,)
        return norm,F(0),F(0)
    assert name=='kite' and sum(propagated)==3*q
    w=norm*F(prod(factorial(v-1) for v in propagated),factorial(q-1))
    return tuple(w*c for c in checked_cut(orbit(propagated)))+(F(0),)


def boundary_rows(graph,u,lookup,count,values,q,cap,edge_reps):
    ctx=u.context()
    for e in edge_reps:
        b=u.derivative(e)
        a=u-ctx.gen(e)*b
        delete,contract=faces(graph,e)
        for vertex,remaining in ((False,delete),(True,contract)):
            moves,name,keep=reduction(remaining)
            if name not in ('edge','kite'):
                continue
            total=0
            for power in range(1,min(cap+2,q)):
                side=power if vertex else q-power
                if side>cap: continue
                factor=(b if vertex else a)**side
                entries=[(tuple(map(int,ex)),int(c)) for ex,c in factor.terms()]
                # The multiplier of B^power (or A^(q-power)) must keep
                # every child in the strictly convergent supplied bank.
                limits=[cap-max(ex[i] for ex,c in entries) for i in range(9)]
                limits[e]=0
                degree=4*q-9-(power-1)-int(factor.total_degree())
                if degree<0 or degree>sum(limits):continue
                beta=F(factorial(power-1)*factorial(q-power-1),factorial(q-1))
                for source in product(*(range(k+1) for k in limits)):
                    if sum(source)!=degree:continue
                    m=list(source)
                    m[e]=power-1
                    children=Counter()
                    for ex,c in entries:
                        target=tuple(i+j for i,j in zip(m,ex))
                        index=lookup.get(target)
                        if index is None:break
                        children[index]+=c
                    else:
                        ex=tuple(v for i,v in enumerate(source) if i!=e)
                        qlower=q-power if vertex else power
                        rhs=tuple(beta*v for v in lower_period(moves,name,keep,ex,qlower))
                        assert all(sum(v*values[j][s] for j,v in children.items())==rhs[s] for s in range(3))
                        den=lcm(*(v.denominator for v in rhs))
                        row={j:v*den for j,v in children.items()}
                        row.update({count+s:-int(v*den) for s,v in enumerate(rhs) if v})
                        total+=1
                        yield primitive(row)
            print('integrated face',e,'contraction' if vertex else 'deletion',name,'identities',total,flush=True)


def six_power_bank():
    here=Path(__file__).resolve().parents[1]/'nonSP/verify/prism_forcer'
    data=json.loads((here/'periods.json').read_text())
    values=[tuple(map(F,row)) for row in data['coefficients']]
    exponents=[tuple(row) for row in data['exponents']]
    lookup={tuple(ex[j] for j in p):i for i,ex in enumerate(exponents) for p in PERMS}
    tree=core_data(tuple(PRISM))[1]
    cotrees=[tuple(map(int,ex)) for ex,c in tree.terms()]
    assert len(cotrees)==75 and len(values)==1481
    scales=[lcm(*(row[i].denominator for row in values)) for i in range(3)]
    scaled=[tuple(int(a*b) for a,b in zip(row,scales)) for row in values]
    rows=[]
    eligible=0
    for m in product(range(3),repeat=9):
        if sum(m)!=11 or m!=fold(m,PERMS):
            continue
        children=[lookup.get(tuple(a+b for a,b in zip(m,t))) for t in cotrees]
        if None in children:
            continue
        eligible+=1
        for row in identities(m,cotrees,children,6):
            assert not any(sum(v*scaled[j][s] for j,v in row.items()) for s in range(3))
            rows.append(row)
    print('eligible lower exponent orbits',eligible,'exact IBP rows',len(rows),flush=True)
    prime=1009
    matrix=nmod_mat([[row.get(j,0)%prime for j in range(len(values))] for row in rows],prime)
    rank=matrix.rank()
    basis=triangular(rows)
    assert (eligible,len(rows),len(basis),rank)==(147,1323,852,852)
    # Dense, asymmetric test coefficients exercise every pivot simultaneously.
    weights={i:(-1)**i*(i+1) for i in range(len(values))}
    reduced=collapse(weights,basis)
    def channels(w):
        return tuple(sum(v*values[j][i] for j,v in w.items()) for i in range(3))
    assert channels(weights)==channels(reduced)
    print('pivot rank',rank,'remaining moments',len(reduced),'largest row',max(map(len,basis.values())),
          'three period channels',channels(reduced),flush=True)
    del matrix
    face=list(boundary_rows(PRISM,tree,lookup,len(values),values,6,3,(0,1)))
    combined=rows+face
    rank=nmod_mat([[row.get(j,0)%prime for j in range(len(values))] for row in combined],prime).rank()
    basis=triangular(combined)
    assert len(face)==111 and len(basis)==rank==877 and max(basis)<len(values)
    reduced=collapse(weights,basis)
    components=channels({j:v for j,v in reduced.items() if j<len(values)})
    assert tuple(v+reduced.get(len(values)+i,0) for i,v in enumerate(components))==channels(weights)
    print('with exact contracted or deleted faces','pivot rank',rank,
          'remaining unknown moments',sum(j<len(values) for j in reduced),flush=True)
    assert sum(j<len(values) for j in reduced)==604
    return basis


def twelve_power_block(name,number):
    graph,group=(PRISM,PERMS) if name=='prism' else (K33,AUTOMORPHISMS)
    folder=Path(__file__).resolve().parent
    filename=folder/('quartic_prism.weights.gz' if name=='prism' else 'quartic_k33.weights.gz')
    tree=core_data(tuple(graph))[1]
    cotrees=[tuple(map(int,ex)) for ex,c in tree.terms()]
    # Each orbit key packs its nine propagator powers in lexicographic order.
    code=np.array([1<<(4*i) for i in reversed(range(9))],dtype=np.int64)
    transformations=np.array([code[np.argsort(p)] for p in group],dtype=np.int64)
    increments=transformations@np.array(cotrees,dtype=np.int64).T
    rng=Random(1739)
    lower={}
    needed=set()
    while len(lower)<number:
        m=[0]*9
        for _ in range(35):
            m[rng.randrange(9)]+=1
        if max(m)>9:continue
        starts=transformations@(np.array(m,dtype=np.int64)+1)
        key=int(starts.min())
        if key in lower:continue
        children=[int(v) for v in (starts[:,None]+increments).min(axis=0)]
        lower[key]=(m,children)
        needed.update(children)
    coefficients=[{},{}]
    with open_gzip(filename,'rt') as stream:
        for line in stream:
            *fields,k,j=line.split()
            p=tuple(map(int,fields))
            key=sum(v<<(32-4*i) for i,v in enumerate(p))
            if key in needed:
                factor=F(factorial(11),prod(factorial(v-1) for v in p))
                for table,entry in zip(coefficients,(k,j)):
                    if entry!='0': table[key]=F(entry)*factor
    seen=set.union(*(set(t) for t in coefficients))
    rows=[]
    eligible=0
    for m,children in lower.values():
        # When every positive tree summand is in the published convergent
        # bank, the lower monomial is convergent by dividing their sum by U.
        if all(k in seen for k in children):
            eligible+=1
            rows.extend(identities(m,cotrees,children,12))
    basis=triangular(rows)
    print(name,'lower candidates',number,'eligible',eligible,'distinct pivots',len(basis),flush=True)
    checks=[]
    for channel,original in zip(('K4','J4'),coefficients):
        reduced=collapse(original,basis)
        # Replay the exact correction in reverse to recover every original
        # weight. Record the pivot multipliers during an independent sweep.
        intermediate=dict(original)
        multipliers={}
        for lead,row in sorted(basis.items()):
            w=intermediate.get(lead,0)
            if w:
                multipliers[lead]=w/row[lead]
                for k,c in row.items():
                    intermediate[k]=intermediate.get(k,0)-multipliers[lead]*c
        assert {k:v for k,v in intermediate.items() if v}==reduced
        rebuilt=dict(reduced)
        for lead,w in multipliers.items():
            for k,c in basis[lead].items():
                rebuilt[k]=rebuilt.get(k,0)+w*c
        assert {k:v for k,v in rebuilt.items() if v}==original
        digest=sha256()
        for k,w in sorted(reduced.items()):
            digest.update(f'{k} {w}\n'.encode('ascii'))
        result=(len(original),len(reduced),digest.hexdigest())
        print(channel,'affected input/output',*result,flush=True)
        checks.append(result)
    if number==256:
        expected={
            'prism':(160,1277,[(13186,11871,'fbb082b0bbed9f550a9b7a0a77849e87f80585568186ef7724e809b5892afa45'),
                               (13183,11868,'2abc4c244213817ec85351bbd33faa8e5f64d0741c2f9f8ce363c53448bdb912')]),
            'K33':(231,1818,[(15232,13409,'11a2910d90579fae2fc1e5258c0a0dae79ba7f4f73b6b1056422c4eea049f3b3'),
                             (15232,13408,'9f915ab63e0259e432df0967af1261bca8010772de9091a6b4625f76eaea6f4e')]),
        }
        assert (eligible,len(basis),checks)==expected[name]
    return checks


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--physical-block',choices=['prism','K33'])
    parser.add_argument('--lower-orbits',type=int,default=256)
    args=parser.parse_args()
    if args.physical_block:
        twelve_power_block(args.physical_block,args.lower_orbits)
    else:
        six_power_bank()
