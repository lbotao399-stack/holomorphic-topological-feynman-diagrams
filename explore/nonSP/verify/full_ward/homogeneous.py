"""Rank and rational kernel of the complete homogeneous four-loop Ward system."""
from math import gcd, isqrt, lcm
from fractions import Fraction
import numpy as np
from flint import fmpz_mat, nmod_mat
from scipy.sparse import csc_matrix

from contractions import basis, contractions, edge_order, edge_sign, nx, pullback, record, types

PRIME = 65521


def signed_basis(graph):
    order = edge_order(graph)
    projector = csc_matrix((490,490), dtype=np.int64)
    automorphisms = list(nx.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter())
    for p in automorphisms:
        # The quotient map is also the map on labels for a bijection.
        projector += edge_sign(order,p)*pullback(tuple(p[i] for i in range(6)))
    assert (projector @ projector - len(automorphisms)*projector).nnz==0
    rank,rem = divmod(int(projector.diagonal().sum()),len(automorphisms))
    assert rem==0
    reduced,modrank = nmod_mat((projector.toarray()%PRIME).tolist(),PRIME).rref()
    assert modrank==rank
    pivots = []
    for row in range(rank):
        start = pivots[-1]+1 if pivots else 0
        pivots.append(next(i for i in range(start,490) if reduced[row,i]))
    return projector[:,pivots]


def blocks(graphs, invariant):
    from scipy.sparse import hstack
    for graph in types(7):
        terms = [csc_matrix((len(basis(6)),m.shape[1]),dtype=np.int64) for m in invariant]
        for index,weight,image in contractions(graph,graphs):
            terms[index] += weight*(pullback(tuple(image[v] for v in range(7))) @ invariant[index])
        yield hstack(terms,format='csr')


def nullspace(matrix,prime):
    reduced,rank = nmod_mat(matrix.tolist(),prime).rref()
    pivots = []
    for i in range(rank):
        start = pivots[-1]+1 if pivots else 0
        pivots.append(next(j for j in range(start,reduced.ncols()) if reduced[i,j]))
    free = np.setdiff1d(np.arange(reduced.ncols()),pivots).astype(int)
    vectors = np.zeros((reduced.ncols(),len(free)),dtype=np.int64)
    vectors[free,:] = np.eye(len(free),dtype=np.int64)
    for i,j in enumerate(pivots):
        vectors[j] = [(-int(reduced[i,int(k)]))%prime for k in free]
    return vectors,free


def modular_kernel(matrices,prime):
    from scipy.sparse import vstack
    rng = np.random.default_rng(481701)
    sample = []
    for m in matrices:
        nz = np.flatnonzero(np.diff(m.indptr))
        if len(nz):
            sample.append(m[rng.choice(nz,size=min(78,len(nz)),replace=False)])
    matrix = vstack(sample).toarray()
    vectors,free = nullspace(matrix,prime)
    print('sample nullity',len(free),flush=True)
    for m in matrices:
        rem = (m @ vectors)%prime
        if not rem.any():
            continue
        update,keep = nullspace(rem,prime)
        # Modular products of basis matrices may exceed int64 at larger primes.
        if len(free)*(prime-1)**2 < 2**63:
            vectors = (vectors @ update)%prime
        else:
            vectors = ((vectors.astype(object) @ update.astype(object))%prime).astype(np.int64)
        free = free[keep]
    assert np.array_equal(vectors[free],np.eye(len(free),dtype=np.int64))
    return free,vectors


def reconstruct(value,prime):
    bound = isqrt(prime//2)
    a,b,c,d = prime,value,0,1
    while abs(b)>bound:
        q,r=divmod(a,b)
        a,b,c,d = b,r,d,c-q*d
    if d==0 or abs(d)>bound or gcd(b,d)!=1:
        raise ValueError('Insufficient modulus for rational reconstruction')
    return Fraction(b,d)


def exact_kernel():
    _,graphs = record()
    invariant = [signed_basis(g) for g in graphs]
    dims = [x.shape[1] for x in invariant]
    print('signed source dimensions',dims,flush=True)
    matrices = list(blocks(graphs,invariant))
    modulus = 1
    crt = None
    for prime in (65521,65519,65497,65479,65449,65447):
        free,vectors = modular_kernel(matrices,prime)
        if crt is None:
            crt = np.zeros(vectors.shape,dtype=object)
            reference = free
        assert np.array_equal(free,reference)
        crt += modulus*(((vectors-crt)*pow(modulus,-1,prime))%prime)
        modulus *= prime
        print('full modular nullity',len(free),'modulus',modulus,flush=True)
        try:
            lift = np.vectorize(lambda v:reconstruct(int(v),modulus),otypes=[object])(crt)
        except ValueError:
            continue
        common = [lcm(*(q.denominator for q in col)) for col in lift.T]
        lift = np.array([[int(v*d) for v,d in zip(row,common)] for row in lift],dtype=object)
        bound = max(abs(int(x)) for x in lift.flat)*max(int(abs(m).sum(axis=1).max()) for m in matrices)
        if bound>=2**63:
            continue
        lift = lift.astype(np.int64)
        if all(not np.any(m @ lift) for m in matrices):
            assert np.array_equal(lift[free],np.diag(common))
            print('rational rank',sum(dims)-len(free),'nullity',len(free),'free',list(map(int,free)),flush=True)
            return invariant,matrices,free,lift
    raise ArithmeticError('The rational kernel could not be certified.')


def restrictions(invariant,lift):
    offsets=np.cumsum([0]+[m.shape[1] for m in invariant])
    slices=[lift[a:b] for a,b in zip(offsets,offsets[1:])]
    def rank(ids):
        return fmpz_mat(np.concatenate([slices[i] for i in ids]).tolist()).rank()
    from itertools import combinations
    ranks = [rank([i]) for i in range(13)]
    print('one-graph ranks',ranks,flush=True)
    sp = (0,1,2,3,5)
    kite = (4,6,7,9)
    wheel = (8,10)
    for group in (sp,kite,sp+kite,wheel,(11,12),wheel+(11,12),sp+kite+wheel):
        print('restricted rank',group,rank(group),flush=True)
    for size in range(1,5):
        full=[g for g in combinations(range(13),size) if rank(g)==lift.shape[1]]
        if full:
            print('smallest determining sets',full,flush=True)
            break
    return ranks


if __name__=='__main__':
    invariant,_,_,lift=exact_kernel()
    restrictions(invariant,lift)
