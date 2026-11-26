"""Edge contractions in the standard bracket basis, with all sources retained."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product
import json
from math import lcm
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy.sparse import csc_matrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from four_mark_relations import contractions, edge_order, edge_sign, types
from sampling import PAIRS, SEEDS


@cache
def basis(size):
    pairs = tuple(combinations(range(size), 2))
    return tuple(s for s in combinations_with_replacement(pairs, 4)
                 if all(i[1] <= j[1] for i, j in zip(s, s[1:])))


@cache
def straighten(columns):
    """For i<j<k<l, [il][jk]=[ik][jl]-[ij][kl]."""
    for pos, ((i,l),(j,k)) in enumerate(zip(columns, columns[1:])):
        if l <= k:
            continue
        assert i < j < k < l
        rest = columns[:pos] + columns[pos+2:]
        result = Counter(dict(straighten(tuple(sorted(rest+((i,k),(j,l)))))))
        result.subtract(dict(straighten(tuple(sorted(rest+((i,j),(k,l)))))))
        return tuple((s,v) for s,v in result.items() if v)
    return ((columns,1),)


def pullback(image):
    """image[v] is the quotient label of old vertex v. Eliminate the last source."""
    size = len(image)-1
    source = tuple(tuple(int(image[v]==j)-int(image[-1]==j) for v in range(size)) for j in range(5))
    return substitution(source)


@cache
def substitution(source):
    """Rows give the first five project sources in the independent input momenta."""
    size = len(source[0])
    target = {s:i for i,s in enumerate(basis(size))}
    brackets = {}
    for i,j in PAIRS:
        terms = []
        for u,v in combinations(range(size),2):
            w = source[i][u]*source[j][v] - source[i][v]*source[j][u]
            if w:
                terms.append(((u,v),w))
        brackets[i,j] = terms
    indices, data, ptr = [], [], [0]
    for monomial in basis(5):
        row = Counter()
        for (a,w),(b,x),(c,y),(d,z) in product(*(brackets[p] for p in monomial)):
            factor = w*x*y*z
            for s,val in straighten(tuple(sorted((a,b,c,d)))):
                row[target[s]] += val*factor
        for s,val in row.items():
            if val:
                indices.append(s)
                data.append(val)
        ptr.append(len(indices))
    return csc_matrix((data,indices,ptr),shape=(len(target),490),dtype=np.int64)


def record():
    data = json.loads((HERE/'polynomials.json').read_text())
    assert tuple(tuple(PAIRS[i] for i in s) for s in SEEDS)==basis(5)
    assert data['pairs']==list(map(list,PAIRS)) and data['seeds']==list(map(list,SEEDS))
    return data, [nx.Graph(edges) for edges in data['edges']]


def scaled_coefficients(data):
    rows = [{int(i):tuple(map(Fraction,coeff)) for i,coeff in graph.items()}
            for graph in data['coefficients']]
    den = lcm(*(v.denominator for row in rows for values in row.values() for v in values))
    # Apply small integer pullbacks to scaled coefficients with Python integers.
    arrays = []
    for row in rows:
        array = np.zeros((490,3),dtype=object)
        for i,vals in row.items():
            array[i] = [int(v*den) for v in vals]
        arrays.append(array)
    return den, arrays


def exact_product(matrix, array):
    # Choose native integers only after a global row-sum bound.
    bound = max(abs(int(x)) for x in array.flat)*int(abs(matrix).sum(axis=1).max())
    if bound <= np.iinfo(np.int64).max:
        return (matrix @ array.astype(np.int64)).astype(object)
    return matrix.toarray().astype(object) @ array


def lower_faces(graph):
    edges = edge_order(graph)
    return [cut for k in (3,4,5) for cut in combinations(sorted(graph),k)
            if sum(i in cut and j in cut for i,j in edges)==2*k-3]


def check():
    data, graphs = record()
    den, kernel = scaled_coefficients(data)
    pure, nonzero, edges = [], [], 0
    assert len(basis(6))==1764
    for i,graph in enumerate(types(7)):
        residual = np.zeros((1764,3),dtype=object)
        for number,sign,image in contractions(graph,graphs):
            mat = pullback(tuple(image[v] for v in range(7)))
            residual += sign*exact_product(mat,kernel[number])
            edges += 1
        assert not any(residual[:,2]),i
        if not lower_faces(graph):
            assert not np.any(residual),i
            pure.append(i)
        nonzero.append([np.count_nonzero(x) for x in residual.T])
    print('signed contractions',edges,'pure graphs',pure,'nonzero graph residuals per channel',
          np.count_nonzero(nonzero,axis=0).tolist(),'common denominator',den)


if __name__=='__main__':
    check()
