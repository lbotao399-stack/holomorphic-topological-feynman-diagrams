"""Transverse cubic consistency on a single tuple, retaining every distinct fusion."""
from itertools import combinations
from math import comb
from random import Random

import numpy as np
from flint import nmod_mat

from soft_consistency import cochain_row, transverse_vector
from soft_target_basis import det

PRIME = 1009
RANK_PATTERNS = (
    (4,3), (4,2,1), (4,1,1,1), (3,3,1), (3,2,2), (3,2,1,1),
    (3,1,1,1,1), (2,2,2,1), (2,2,1,1,1), (2,1,1,1,1,1), (1,)*7,
)


def fixture(n, seed):
    rng = Random(seed)
    waves = [tuple(rng.randrange(-7,8) for _ in range(5)) for _ in range(n-1)]
    waves.append(tuple(-sum(col) for col in zip(*waves)))
    return waves


def coordinates(waves):
    """Return a coordinate slice for potentials modulo the common and longitudinal shifts."""
    n = len(waves)
    gauges = np.zeros((5+n,n,5),dtype=np.int64)
    for i in range(5):
        gauges[i,:,i] = 1
    for i,q in enumerate(waves):
        gauges[5+i,i] = q
    rows,rank = nmod_mat(gauges.reshape(n+5,-1).tolist(),PRIME).rref()
    assert rank == n+5
    reduced = np.asarray(rows.tolist(),dtype=np.int64)
    pivots = np.argmax(reduced!=0,axis=1)
    free = np.setdiff1d(np.arange(5*n),pivots)
    return pivots,free,reduced


class Fiber:
    def __init__(self, p):
        self.p = np.array(p,dtype=np.int64)
        self.n = len(p)
        self.pairs = list(combinations(range(self.n),2))
        self.partition_ids = {}
        blocks = []
        reciprocity = []
        all_free = []
        rows_after = []
        for pair in self.pairs:
            others = [i for i in range(self.n) if i not in pair]
            group = [(i,) for i in others]+[pair]
            rows = [p[i] for i in others]+[tuple(self.p[list(pair)].sum(axis=0))]
            pivots,free,gauge = coordinates(rows)
            all_free.append(free)
            lookup = [others.index(i) if i in others else len(others) for i in range(self.n)]
            rows_after.append(lookup)
            columns = []
            values = []
            for i,j in combinations(range(len(rows)),2):
                partition = tuple(sorted([group[k] for k in range(len(group)) if k not in (i,j)]
                                         + [tuple(sorted(group[i]+group[j]))]))
                index = self.partition_ids.setdefault(partition,len(self.partition_ids))
                raw = np.zeros((len(rows),5),dtype=np.int64)
                raw[i] = rows[j]
                raw[j] = rows[i]
                flat = raw.flatten()
                values.append((flat[free] - flat[pivots]@gauge[:,free])%PRIME)
                columns.append(index)
            blocks.append((np.array(columns),np.array(values)))
            rec = []
            for i,j in combinations(range(5),2):
                raw = np.zeros((len(rows),5),dtype=np.int64)
                raw[:,j] = [q[i] for q in rows]
                raw[:,i] = [-q[j] for q in rows]
                rec.append(raw.flatten()[free])
            reciprocity.append(rec)
        self.free = np.array(all_free)
        self.blocks = blocks
        self.reciprocity = np.array(reciprocity)
        self.rows_after = np.array(rows_after)
        self.width = len(self.free[0])
        assert len(self.partition_ids) == comb(self.n,3)+3*comb(self.n,4)

    def trace_map(self):
        out = np.zeros((len(self.pairs),self.width,len(self.partition_ids)),dtype=np.int64)
        for i,(columns,block) in enumerate(self.blocks):
            out[i][:,columns] = block.T
        return out.reshape(-1,len(self.partition_ids))

    def test(self, positions, vectors):
        # Seven transverse vectors are wedged within their respective slots.
        # A pair fusion removes one vector by divergence; T_A removes a
        # second one. The other five give their signed volume.
        coeff = np.zeros((7,7),dtype=np.int64)
        for i,j in combinations(range(7),2):
            value = det([v for h,v in enumerate(vectors) if h not in (i,j)]) % PRIME
            coeff[i,j] = (-1)**(i+j-1)*value
            coeff[j,i] = -coeff[i,j]
        a,b = np.asarray(self.pairs).T
        positions = np.asarray(positions)
        own = (positions==a[:,None]) | (positions==b[:,None])
        first = (((self.p[a]+self.p[b])@np.asarray(vectors).T)%PRIME)*own
        weight = first@coeff % PRIME
        blocks = np.zeros((len(a),self.n-1,5),dtype=np.int64)
        for j,v in enumerate(vectors):
            blocks[np.arange(len(a)),self.rows_after[:,positions[j]]] += weight[:,j,None]*np.array(v)
        return np.take_along_axis(blocks.reshape(len(a),-1),self.free,axis=1)%PRIME

    def test_rows(self,count,seed):
        rng = Random(seed)
        result = []
        for trial in range(count):
            pattern = RANK_PATTERNS[trial%len(RANK_PATTERNS)]
            labels = rng.sample(range(self.n),len(pattern))
            positions = sorted(i for i,count_here in zip(labels,pattern) for _ in range(count_here))
            vectors = []
            for i in positions:
                q = self.p[i].tolist()
                k = next(j for j,x in enumerate(q) if x)
                z = [rng.randrange(-4,5) if j!=k else 0 for j in range(5)]
                v = [q[k]*x for x in z]
                v[k] = -sum(x*y for x,y in zip(q,z))
                vectors.append(v)
            result.append(self.test(positions,vectors).flatten())
        return np.array(result)


def modular_matrix(matrix):
    # Populate the sparse entries without a second, boxed copy of a dense matrix.
    n,m = matrix.shape
    out = nmod_mat(n,m,PRIME)
    for i,row in enumerate(matrix):
        cols = np.flatnonzero(row)
        for j,x in zip(cols,row[cols]):
            out[i,int(j)] = int(x)
    return out


def modular_rank(matrix):
    return modular_matrix(matrix).rank()


def kernel_basis(matrix):
    basis,dim = nmod_mat((matrix%PRIME).tolist(),PRIME).nullspace()
    return np.array(basis.tolist(),dtype=np.int64)[:,:dim]


def transverse_sum_zero(p):
    n = len(p)
    constraints = np.zeros((n+5,n,5),dtype=np.int64)
    for i,row in enumerate(p):
        constraints[i,i] = row
    for j in range(5):
        constraints[n+j,:,j] = 1
    return kernel_basis(constraints.reshape(n+5,-1)).reshape(n,5,-1)


def exterior_consistency(fiber):
    motions = transverse_sum_zero(fiber.p)
    n,_,dim = motions.shape
    # On the sum-zero transverse space, multiply the fusion covector and
    # the evaluation of the potential as an exterior two-form.
    a,b = np.triu_indices(dim,k=1)
    pieces = []
    for index,(i,j) in enumerate(fiber.pairs):
        linear = (fiber.p[i]@motions[j]+fiber.p[j]@motions[i]) % PRIME
        other = [h for h in range(n) if h not in (i,j)]
        fused = np.concatenate((motions[other],(motions[i]+motions[j])[None,:,:]))
        marked = fused.reshape(-1,dim)[fiber.free[index]]
        block = (linear[a,None]*marked.T[b] - linear[b,None]*marked.T[a]) % PRIME
        pieces.append(block)
    return np.concatenate(pieces,axis=1)


def check_exterior_sign():
    p = fixture(11,78965)
    fiber = Fiber(p)
    rng = Random(72)
    anchor = next(i for i,q in enumerate(p) if q[0])
    legs = [i for i in (1,3,5,7,8) if i!=anchor][:3]
    vectors = [transverse_vector(p[i],rng) for i in legs]
    columns = [(anchor,tuple(-p[anchor][j] if k==0 else p[anchor][0]*int(k==j)
                             for k in range(5))) for j in range(1,5)]
    columns = sorted(columns+list(zip(legs,vectors)), key=lambda field:field[0])
    positions,polarizations = zip(*columns)
    raw = np.zeros((55,50),dtype=np.int64)
    for index,value in cochain_row(p,anchor,legs,vectors):
        raw[index//50,index%50] += value
    expected = np.take_along_axis(raw,fiber.free,axis=1)*p[anchor][0]**3%PRIME
    assert np.any(expected) and np.array_equal(fiber.test(positions,polarizations),expected)


def check_gluing(fiber,trace,rank):
    # Up to relabeling, these exhaust the simple graphs with one, two, or
    # three edges. Prescribe all scalar values on the indicated fusion faces.
    graphs = [
        [(0,1)], [(0,1),(0,2)], [(0,1),(2,3)],
        [(0,1),(0,2),(1,2)], [(0,1),(1,2),(2,3)], [(0,1),(0,2),(0,3)],
        [(0,1),(0,2),(3,4)], [(0,1),(2,3),(4,5)],
        *[[(0,i) for i in range(1,size+1)] for size in range(4,10)],
    ]
    for pairs in graphs:
        ids = [fiber.pairs.index(pair) for pair in pairs]
        columns = np.unique(np.concatenate([fiber.blocks[i][0] for i in ids]))
        rest = np.ones(trace.shape[1],dtype=bool)
        rest[columns] = False
        local = trace.reshape(len(fiber.pairs),fiber.width,-1)[ids][:,:,columns]
        local_rank = modular_rank(local.reshape(-1,len(columns)))
        outside = modular_rank(trace[:,rest])
        assert local_rank==33*len(pairs) and outside==rank-local_rank
        print('gluing',pairs,'shared scalar values',len(columns),'free choices',len(columns)-local_rank,
              flush=True)


def check(n,seed,count,glue=False):
    fiber = Fiber(fixture(n,seed))
    trace = fiber.trace_map()
    rank = modular_rank(trace)
    print(n,'scalar partitions',trace.shape[1],'potential coordinates',trace.shape[0],
          'trace rank',rank,flush=True)
    ce = exterior_consistency(fiber)
    rows = ce.reshape(ce.shape[0],len(fiber.pairs),fiber.width)
    restricted = []
    boundary_rank = 0
    for i,(columns,values) in enumerate(fiber.blocks):
        rec = fiber.reciprocity[i]
        assert not np.any((rec@values.T)%PRIME)
        tangent = kernel_basis(rec)
        boundary_rank += fiber.width-tangent.shape[1]
        restricted.append((rows[:,i]@tangent)%PRIME)
    outgoing = modular_rank(np.concatenate(restricted,axis=1))
    if n >= 7:
        assert outgoing == comb(4*n-15,2)-comb(5,3)
    # Check every elementary composition, including an independent set
    # of decomposable test polyvectors with all eleven rank partitions.
    tests = np.concatenate((ce,fiber.test_rows(count,seed+1))) if count else ce
    for start in range(0,len(tests),100):
        batch = tests[start:start+100].reshape(-1,len(fiber.pairs),fiber.width)
        product = np.zeros((len(batch),trace.shape[1]),dtype=np.int64)
        for i,(columns,values) in enumerate(fiber.blocks):
            product[:,columns] += batch[:,i]@values.T
        assert not np.any(product%PRIME)
    print(n,'reciprocity',boundary_rank,'consistent two-form rank',outgoing,
          'kernel minus image',trace.shape[0]-rank-outgoing-boundary_rank,flush=True)
    expected = {9:(412,360,200), 11:(979,550,396), 13:(1918,780,656), 15:(3325,1050,980)}
    if n in expected:
        assert (rank,boundary_rank,outgoing)==expected[n]
        assert trace.shape[0]==rank+boundary_rank+outgoing
    if glue:
        assert n==13
        check_gluing(fiber,trace,rank)
    return rank,boundary_rank,outgoing


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n',type=int)
    parser.add_argument('--seed',type=int,default=78965)
    parser.add_argument('--tests',type=int,default=22)
    parser.add_argument('--glue',action='store_true')
    args = parser.parse_args()
    check_exterior_sign()
    check(args.n,args.seed,args.tests,args.glue)
