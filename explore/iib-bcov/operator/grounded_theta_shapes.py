"""Unisolvence and individual theta shapes with an arbitrary grounded charge."""
from itertools import combinations, permutations
from pathlib import Path
from random import Random
import sys

from flint import fmpz_mat, nmod_mat
import numpy as np

from genus_two_probes import orbit_data
from soft_trace_dimensions import multiply_by_slot


WORD_ROWS = (2,3,4,7,11,12,16,18,19,24,25,28,32,36,40,41,44,45,51,57,66,79,80,81,
             104,105,106,117,141,156,173,200,222,286,310,383,390,601)
PROBES = (*range(32),35,40,41,42,80,81)


def dimension():
    shape = (5,4,4,4,4)
    size = tuple(weight+5-i for i,weight in enumerate(shape))
    unit = np.zeros(size,dtype=object)
    unit[(0,)*5] = 1
    sym = [unit]
    for n in range(1,8):
        total = sum(multiply_by_slot(sym[n-j],j,2) for j in range(1,n+1))
        sym.append(total//n)
    # A distinguished output slot, seven indistinguishable slots, all double soft.
    # The eight independent momenta entail no total-momentum Koszul factors.
    a = multiply_by_slot(sym[7],1,2)
    return sum((-1)**sum(i>j for i,j in combinations(perm,2))
               *a[tuple(w-i+j for i,(w,j) in enumerate(zip(shape,perm)))]
               for perm in permutations(range(5)))


def momenta(seed):
    rng = Random(seed)
    return tuple(tuple(rng.randint(-3,3) for _ in range(5)) for _ in range(8))


def words(subsets):
    rng = Random(981433)
    cover = [sum(1<<i for i in s) for s in subsets]
    while True:
        free = rng.randrange(2)
        kind = rng.randrange(4)
        first,second = rng.randrange(56),rng.randrange(56)
        if kind==0:
            word = (first,first,second,second)
        elif kind==1:
            word = (first,first,second,rng.randrange(56))
        else:
            word = (first,second,rng.randrange(56),rng.randrange(56))
        word = tuple(sorted(word))
        if all((i==free)+sum(cover[k]>>i&1 for k in word)>=2 for i in range(8)):
            yield free,word


def orbit_rows(points):
    subsets,images,pair = orbit_data()
    order = np.asarray(pair,dtype=np.int64)
    mask = np.array(images,dtype=np.int64)
    measures = []
    for p in points:
        values = np.array([0]+[int(fmpz_mat([p[i] for i in s]).det()) for s in subsets])
        signs = np.concatenate((values,-values[:0:-1]))
        measures.append((signs[mask]%1009,np.asarray(p,dtype=np.int64)))
    rows = []
    for index,(free,word) in enumerate(words(subsets)):
        if index not in WORD_ROWS:
            continue
        evals = []
        for d,p in measures:
            weights = (d[word[0]]*d[word[1]]%1009)*(d[word[2]]*d[word[3]]%1009)%1009
            if free:
                coeff = np.bincount(order,weights=weights,minlength=64).astype(np.int64)
                value = coeff.reshape(8,8)@p
            else:
                coeff = np.bincount(order//8,weights=weights,minlength=8).astype(np.int64)
                value = coeff[:,None]*p
            evals.extend((value%1009).flatten().tolist())
        rows.append([evals[j] for j in PROBES])
        if index==WORD_ROWS[-1]:
            return rows


def check():
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
    from bulk_eight import grounded_potentials,SHAPES
    points = [momenta(seed) for seed in (31307,41063,39222)]
    assert dimension()==38
    rows = orbit_rows(points)
    assert int(nmod_mat(rows,1009).det())==397
    print('Independent raw polynomial potentials: 38. Evaluation determinant 397 mod 1009.')
    for p in points:
        assert any(map(sum,zip(*p)))
        columns = [grounded_potentials(p,shapes=(s,)) for s in SHAPES]
        for entries in zip(*columns):
            for a,b,c,d in zip(*entries):
                assert 3*c == -a+2*b and 9*d == -20*a+13*b
    print('All raw off-conservation slot components obey both theta identities.')


if __name__=='__main__':
    check()
