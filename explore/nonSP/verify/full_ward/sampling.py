"""Degree-eight bracket basis on the five independent two-component sources."""
from functools import cache
from itertools import combinations,combinations_with_replacement
from random import Random

import numpy as np
from flint import nmod_mat

PAIRS=tuple(combinations(range(5),2))
SEEDS=tuple(s for s in combinations_with_replacement(range(10),4)
            if all(PAIRS[i][1]<=PAIRS[j][1] for i,j in zip(s,s[1:])))
assert len(SEEDS)==490


def evaluation(probes):
    p=np.asarray(probes,dtype=np.int64)
    assert np.all(p.sum(axis=1)==0)
    b=np.array([p[:,i,0]*p[:,j,1]-p[:,i,1]*p[:,j,0] for i,j in PAIRS])
    return np.prod(b[np.array(SEEDS)],axis=1).T


def pivots(array):
    p=65521
    rows,rank=nmod_mat((array%p).tolist(),p).rref()
    result=[]
    for i in range(rank):
        start=result[-1]+1 if result else 0
        result.append(next(j for j in range(start,rows.ncols()) if rows[i,j]))
    return result


@cache
def sampling():
    rand=Random(617472)
    def axis():
        while True:
            v=[rand.choice((-2,-1,-1,0,1,1,2)) for _ in range(5)]
            if sum(abs(i) for i in v)<=7 and abs(sum(v))<=2:
                return v+[-sum(v)]
    probes=[list(zip(axis(),axis())) for _ in range(720)]
    matrix=evaluation(probes)
    selected=pivots(matrix.T)
    assert len(selected)==490
    extra=[i for i in range(len(probes)) if i not in selected][-2:]
    return [probes[i] for i in selected+extra],matrix[selected+extra,:]
