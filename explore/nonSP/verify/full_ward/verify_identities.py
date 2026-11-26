"""Coefficientwise comparison of all outer edge sums with all lower-loop faces."""
from functools import cache
from fractions import Fraction as F
from itertools import combinations
from math import factorial, lcm, prod
from random import Random
import numpy as np
from flint import fmpz_mat

from contractions import (basis, contractions, edge_order, exact_product, pullback,
                          record, scaled_coefficients, types)
from lower_faces import (SOURCES, canonical_inner, inner_polynomial, sign,
                         standard_basis, tight_faces)
from sampling import pivots
from tree_moments import moment_rule


def evaluate_basis(source,basis_rows):
    pairs=tuple(sorted({pair for s in basis_rows for pair in s}))
    bracket=np.array([source[...,i,0]*source[...,j,1]-source[...,i,1]*source[...,j,0]
                      for i,j in pairs])
    inds=np.array([[pairs.index(p) for p in s] for s in basis_rows])
    return np.moveaxis(np.prod(bracket[inds],axis=1),0,-1)


def probes():
    rng=Random(2184013)
    def axis():
        while True:
            v=[rng.choice((-2,-1,-1,0,1,1,2)) for _ in range(6)]
            if sum(map(abs,v))<=8 and abs(sum(v))<=2:
                return v+[-sum(v)]
    p=np.array([list(zip(axis(),axis())) for _ in range(2030)])
    ev=evaluate_basis(p,basis(6))
    chosen=pivots(ev.T)
    assert len(chosen)==len(basis(6))==1764
    return p[chosen],ev[chosen]


def lagrange(p,axis):
    m=int(sum(axis[0]))
    total=p.sum(axis=1)
    values=[]
    for a in axis:
        col=prod((m*p[:,j]-r*total) for j,k in enumerate(a) for r in range(k))
        values.append(col*(factorial(m)//prod(factorial(k) for k in a)))
    return np.array(values).T,factorial(m)*m**m


def int_matmul(a,b):
    bound=int(np.max(np.abs(a)))*int(np.max(np.abs(b)))*a.shape[1]
    if bound<2**63:
        return a.astype(np.int64,copy=False) @ b.astype(np.int64,copy=False)
    exact=fmpz_mat(a.tolist())*fmpz_mat(b.tolist())
    return np.array([[int(exact[i,j]) for j in range(b.shape[1])] for i in range(a.shape[0])],
                    dtype=object)


def dot_rows(a,b):
    bound=int(np.max(np.abs(a)))*int(np.max(np.abs(b)))*a.shape[1]
    if bound<2**63:
        return np.sum(a.astype(np.int64,copy=False)*b.astype(np.int64,copy=False),axis=1)
    return np.sum(a.astype(object)*b.astype(object),axis=1)


@cache
def exact_inner(edges):
    return [(lcm(*(v.denominator for v in col)),col) for col in inner_polynomial(edges)]


def evaluate_inner(tau,n,edges):
    # tau has one source array per tree-weight lattice point per probe.
    output=[]
    for den,coeff in exact_inner(edges):
        columns=[(s,int(c*den)) for s,c in zip(standard_basis(n),coeff) if c]
        if not columns:
            output.append((1,None))
            continue
        pairs={p for s,_ in columns for p in s}
        bound=2*int(np.max(np.abs(tau)))**2
        if bound**n*sum(abs(c) for _,c in columns)>=2**63:
            tau=tau.astype(object)
        brackets={(i,j):tau[:,:,i,0]*tau[:,:,j,1]-tau[:,:,i,1]*tau[:,:,j,0] for i,j in pairs}
        output.append((den,sum(c*prod(brackets[p] for p in s) for s,c in columns)))
    return output


def face_values(face,source):
    cut,inner,rest,mapping,outer,parity=face
    n=len(cut)-2
    m=4-n
    canonical,vertex_map,outer_sign=canonical_inner(tuple(outer))
    target,inner_map,inner_sign=canonical_inner(tuple(inner))
    axis,tree,beta,rule=moment_rule(canonical)
    # All series and edge reversals use the original directions in the sliding graph.
    merged=np.zeros((len(source),m+2,2),dtype=np.int64)
    for v in range(7):
        merged[:,vertex_map[mapping[v]],:]+=source[:,v,:]
    incidence=np.zeros((len(cut),len(canonical)),dtype=np.int64)
    for (a,b),(i,j) in zip(rest,outer):
        oriented=(vertex_map[i],vertex_map[j])
        e=canonical.index(tuple(sorted(oriented)))
        orient=1 if oriented==canonical[e] else -1
        for h,v in enumerate(cut):
            incidence[inner_map[h],e]+=orient*(int(v==b)-int(v==a))
    shifted=np.einsum('bt,ce,teq->bcq',beta,incidence,tree,optimize=True)
    old=np.argsort(inner_map)
    inner_source=source[:,np.array(cut)[old],:]
    sizes=beta.sum(axis=1)
    prefactor=-(-1)**(m*n)*parity*inner_sign*outer_sign
    norms=[1,1]
    out=np.zeros((len(source),2),dtype=object)
    for start in range(0,len(source),96):
        stop=min(start+96,len(source))
        active=merged[start:stop,:-1,:]
        tx,den=lagrange(active[:,:,0],axis)
        ty,_=lagrange(active[:,:,1],axis)
        basis_products=(tx[:,:,None]*ty[:,None,:]).reshape(len(tx),-1)
        tau=(sizes[None,:,None,None]*inner_source[start:stop,None,:,:]+
             np.einsum('bcq,pqd->pbcd',shifted,active,optimize=True))
        assert not np.any(tau.sum(axis=2))
        inputs=evaluate_inner(tau,n,target)
        for j,(outer_den,weight) in enumerate(rule):
            if not np.any(weight):
                continue
            integrated=int_matmul(basis_products,weight)
            for k,(inner_den,input_values) in enumerate(inputs):
                if input_values is None:
                    continue
                assert j+k<=1
                scale=den**2*inner_den*outer_den
                norms[j+k]=lcm(norms[j+k],scale)
                # Only one split can produce zeta(3) at these loop orders.
                assert norms[j+k]==scale
                out[start:stop,j+k]+=prefactor*dot_rows(integrated,input_values).astype(object)
    return out,norms


def run(fixture=False):
    data,graphs=record()
    den,arrays=scaled_coefficients(data)
    if fixture:
        source=np.array([SOURCES])
        ev=evaluate_basis(source,basis(6))
    else:
        source,ev=probes()
    counts=[0]*3
    for i,g in enumerate(types(7)):
        lhs=np.zeros((1764,3),dtype=object)
        for j,w,image in contractions(g,graphs):
            lhs+=w*exact_product(pullback(tuple(image[v] for v in range(7))),arrays[j])
        assert not np.any(lhs[:,2])
        left=int_matmul(ev,lhs[:,:2]).astype(object)
        rhs=np.zeros_like(left)
        divisors=[1,1]
        for face in tight_faces(edge_order(g)):
            values,norm=face_values(face,source)
            counts[len(face[0])-3]+=1
            for k,d in enumerate(norm):
                new=lcm(divisors[k],d)
                rhs[:,k]=rhs[:,k]*(new//divisors[k])+values[:,k]*(new//d)
                divisors[k]=new
        assert np.array_equal(left*divisors,rhs*den),i
        print('verified',i,'sources',len(source),flush=True)
    assert counts==[41,21,22]
    print('all 60 Ward polynomials; nontrivial lower faces',counts,flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture',action='store_true')
    run(parser.parse_args().fixture)
