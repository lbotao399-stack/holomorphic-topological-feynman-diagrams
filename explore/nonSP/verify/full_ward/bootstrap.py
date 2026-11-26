"""Recover all thirteen bare kernels from lower faces and analytically evaluated radial wheels."""
from fractions import Fraction as F
import json
from math import lcm
from random import Random

import numpy as np
from flint import fmpz_mat, nmod_mat

from verify_identities import evaluate_basis, face_values
from contractions import HERE, basis, edge_order, substitution, types
from homogeneous import PRIME, blocks, nullspace, reconstruct, signed_basis
from lower_faces import tight_faces
from radial_wheel import PROBES, recover
from sampling import pivots


def pool():
    rng=Random(171406)
    def axis():
        while True:
            v=[rng.choice((-2,-1,-1,0,1,1,2)) for _ in range(6)]
            if sum(map(abs,v))<=9 and abs(sum(v))<=3:
                return v+[-sum(v)]
    return np.array([list(zip(axis(),axis())) for _ in range(320)],dtype=np.int64)


def choose_equations(invariant):
    graphs=types(6)
    matrices=list(blocks(graphs,invariant))
    source=pool()
    ev=evaluate_basis(source,basis(6))
    rng=np.random.default_rng(10251)
    picks=[rng.permutation(len(source)) for _ in matrices]
    def sampled(at,until):
        rows,keys=[],[]
        for g,(matrix,indices) in enumerate(zip(matrices,picks)):
            if matrix.nnz:
                ids=indices[at:until]
                rows.append(matrix.T.dot(ev[ids].T).T)
                keys.extend((g,j) for j in ids)
        return np.vstack(rows),keys
    rows,keys=sampled(0,72)
    chosen=pivots(rows.T)
    rows=rows[chosen]
    keys=[keys[j] for j in chosen]
    kernel,_=nullspace(rows%PRIME,PRIME)
    def append(part,labels):
        nonlocal rows,keys,kernel
        remainder=((part%PRIME) @ kernel)%PRIME
        new=pivots(remainder.T)
        if new:
            rows=np.vstack([rows,part[new]])
            keys.extend(labels[j] for j in new)
            change,_=nullspace(remainder[new],PRIME)
            kernel=(kernel @ change)%PRIME
    for pos in range(72,len(source),40):
        if kernel.shape[1]==15:
            break
        append(*sampled(pos,pos+40))
    assert kernel.shape[1]==15
    offsets=np.cumsum([0]+[m.shape[1] for m in invariant])
    forbidden=[j for j,s in enumerate(basis(5)) if sum(v in (3,4) for pair in s for v in pair)<2]
    zeros=np.zeros((len(forbidden),rows.shape[1]),dtype=np.int64)
    zeros[:,offsets[9]:offsets[10]]=invariant[9][forbidden].toarray()
    append(zeros,[('zero',j) for j in forbidden])
    assert kernel.shape[1]==10
    ibp=json.loads((HERE/'wheel_ibp.json').read_text())
    assert ibp['catalogue_index']==10 and ibp['edges']==list(map(list,edge_order(graphs[10])))
    functionals=np.array(ibp['relations'],dtype=np.int64)
    probe_basis=evaluate_basis(np.array(ibp['probes'],dtype=np.int64),basis(5))
    tests=np.zeros((len(functionals),rows.shape[1]),dtype=np.int64)
    tests[:,offsets[10]:offsets[11]]=invariant[10].T.dot((functionals @ probe_basis).T).T
    append(tests,[('ibp',j) for j in range(len(tests))])
    assert kernel.shape[1]==2
    tests=np.zeros((len(PROBES),rows.shape[1]),dtype=np.int64)
    wheel_ev=evaluate_basis(np.array(PROBES),basis(5))
    tests[:,offsets[10]:offsets[11]]=invariant[10].T.dot(wheel_ev.T).T
    append(tests,[('wheel',j) for j in range(len(tests))])
    assert rows.shape==(2359,2359) and kernel.shape[1]==0
    print('2344 Ward evaluations, 5 numerator zeros, 8 wheel IBPs, 2 radial evaluations',flush=True)
    return rows,keys,source


def lower_rhs(keys,source):
    radial=recover(PROBES)
    rhs=[[F(),F(),F()] for _ in keys]
    for i,(kind,j) in enumerate(keys):
        if kind=='wheel':
            rhs[i]=list(radial[j])
    for number,g in enumerate(types(7)):
        positions=[i for i,(k,_) in enumerate(keys) if k==number]
        if not positions:
            continue
        probes=source[[keys[i][1] for i in positions]]
        for face in tight_faces(edge_order(g)):
            out,norm=face_values(face,probes)
            for idx,row in zip(positions,out):
                for c in (0,1):
                    rhs[idx][c]+=F(int(row[c]),norm[c])
    den=lcm(*(q.denominator for row in rhs for q in row))
    return den,np.array([[int(q*den) for q in row] for row in rhs],dtype=object)


def solve(matrix,rhs):
    crt=np.zeros(rhs.shape,dtype=object)
    modulus=1
    for p in (1000000007,1000000009,998244353,65521,65519,65497):
        a=nmod_mat((matrix%p).tolist(),p)
        b=nmod_mat((rhs%p).tolist(),p)
        try:
            sol=a.solve(b)
        except ZeroDivisionError:
            continue
        vals=np.array([[int(x) for x in row] for row in sol.tolist()],dtype=object)
        crt+=modulus*(((vals-crt)*pow(modulus,-1,p))%p)
        modulus*=p
        try:
            lift=np.vectorize(lambda x:reconstruct(int(x),modulus),otypes=[object])(crt)
        except ValueError:
            continue
        d=lcm(*(v.denominator for v in lift.flat))
        integer=fmpz_mat([[int(q*d) for q in row] for row in lift])
        if fmpz_mat(matrix.tolist())*integer==d*fmpz_mat(rhs.tolist()):
            return integer,d
    raise ArithmeticError('Failed to reconstruct the rational solution.')


def compute():
    from product_forests import placements
    invariant=[signed_basis(g) for g in types(6)]
    matrix,keys,source=choose_equations(invariant)
    den,rhs=lower_rhs(keys,source)
    physical_rhs=rhs[:,:2].copy()
    physical_rhs[[i for i,(kind,_) in enumerate(keys) if kind=='wheel']]=0
    integer,d=solve(matrix,np.column_stack([rhs,physical_rhs]))
    coeff=[]
    pos=0
    for graph,projection in zip(types(6),invariant):
        width=projection.shape[1]
        sub=fmpz_mat([[integer[i,j] for j in range(5)] for i in range(pos,pos+width)])
        standard=fmpz_mat(projection.toarray().tolist())*sub
        for a,b,roots in placements(graph):
            images=tuple(tuple(int(v==j)-int(v==b) for j in [a]+roots) for v in range(5))
            physical=fmpz_mat(substitution(images).toarray().tolist())*standard
            for i in range(physical.nrows()):
                assert physical[i,0]==physical[i,3]
                assert all(physical[i,j]==0 for j in (1,2,4))
        coeff.append({str(i):[str(F(int(standard[i,j]),d*den)) for j in range(3)]
                      for i in range(490) if any(standard[i,j] for j in range(3))})
        pos+=width
    return coeff


if __name__=='__main__':
    coeff=compute()
    assert coeff==json.loads((HERE/'polynomials.json').read_text())['coefficients']
    print('All 3892 nonzero coefficient triples recovered with both wheel periods evaluated analytically.',flush=True)
    print('All 23 cubic product kernels recovered from the lower-loop terms with zero wheel data.',flush=True)
