"""Bare restricted and projected theta graphs in the defining trace pairing."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path
import sys

import numpy as np
import sympy as sp

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"o3/code"))
from ribbon import local_orders, orientation_sign, perm_sign


def units(n,i,j):
    result = np.zeros((n,n),dtype=np.int64)
    result[i,j] = 1
    return result


def colours(kind):
    n = 5 if kind == "so" else 4
    e = lambda i,j: units(n,i,j)
    if kind == "so":
        basis = [e(i,j)-e(j,i) for i,j in combinations(range(n),2)]
        sharp = lambda x: x.T
        commuting = (e(0,1)-e(1,0),e(2,3)-e(3,2))
        jacobi = (commuting[0],e(1,2)-e(2,1))
    else:
        r = n//2
        J = np.block([[np.zeros((r,r),dtype=int),np.eye(r,dtype=int)],
                      [-np.eye(r,dtype=int),np.zeros((r,r),dtype=int)]])
        sharp = lambda x: -J@x.T@J
        basis = [e(i,j)-e(r+j,r+i) for i,j in product(range(r), repeat=2)]
        for i in range(r):
            for j in range(i,r):
                basis.extend((e(i,r+j)+(e(j,r+i) if i!=j else 0),
                              e(r+i,j)+(e(r+j,i) if i!=j else 0)))
        commuting = (e(0,0)-e(2,2),e(1,1)-e(3,3))
        jacobi = (commuting[0],e(0,2))
    assert all(np.array_equal(x,-sharp(x)) for x in basis)
    return basis,sharp,commuting,jacobi


def comm(x,y):
    return x@y-y@x


def det(p,q):
    return p[0]*q[1]-p[1]*q[0]


def theta_colour(shape,basis,commuting):
    n = len(basis[0])
    dim = len(basis)
    gram = sp.Matrix([[int(np.trace(x@y)) for y in basis] for x in basis])
    twice_inverse = 2*gram.inv()
    assert all(c.is_Integer for c in twice_inverse)
    dual = np.einsum("ij,jmn->imn",np.array(twice_inverse,dtype=np.int64),basis)
    assert np.array_equal(np.einsum("iab,jba->ij",basis,dual),2*np.eye(dim,dtype=int))
    ends = (np.array(basis),dual)
    ext = (1,2,3,4)
    pos = 0
    edges = []
    for length in shape:
        path = (0,*ext[pos:pos+length-1],5)
        edges.extend(zip(path,path[1:]))
        pos += length-1
    # Each defining matrix has operator norm at most one before doubling its dual.
    # This bound controls every partial integer contraction for both test ranks.
    assert dim**7 * (2*n)**6 * 2**7 < 2**63
    tensors = []
    for v,slots in enumerate(local_orders(6,edges,ext)):
        inputs = [np.array(commuting) if slot[0]=="x" else ends[slot[2]] for slot in slots]
        xyz = np.einsum("iab,jbc,kca->ijk",*inputs,optimize=True)
        xzy = np.einsum("iab,kbc,jca->ijk",inputs[0],inputs[2],inputs[1],optimize=True)
        labels = [7+ext.index(v) if h[0]=="x" else h[1] for h in slots]
        tensors.extend([xyz-xzy,labels])
    signed = orientation_sign(6,edges,ext)
    # Four external axes select either of the two commuting matrices.
    return np.einsum(*tensors,[7,8,9,10],optimize=True),F(signed,2**7)


def bare(basis,commuting,momenta):
    external = (0,0,1,1)
    terms = []
    for shape,factor in (((1,3,3),F(-1,20)),((2,2,3),F(1,30))):
        tensor,scale = theta_colour(shape,basis,commuting)
        result = 0
        for i,j,k,l in permutations(range(4)):
            signs = perm_sign((i,j,k,l))*det(momenta[i],momenta[j])*det(momenta[k],momenta[l])
            result += signs*int(tensor[external[i],external[j],external[k],external[l]])
        terms.append(factor*scale*result/4)
    return terms


def pole_exchange(sharp,commuting,momenta):
    xs = (*([commuting[0]]*2),*([commuting[1]]*2))
    n = len(xs[0])
    total = 0
    for i,j in combinations(range(4),2):
        k,l = [s for s in range(4) if s not in (i,j)]
        first = xs[i]@xs[j]+xs[j]@xs[i]
        second = xs[k]@xs[l]+xs[l]@xs[k]
        assert np.array_equal(first,sharp(first)) and np.array_equal(second,sharp(second))
        pairing = F(int(np.trace(first@second)),2)-F(int(np.trace(first))*int(np.trace(second)),2*n)
        total += (-1)**(i+j-1)*det(momenta[i],momenta[j])*det(momenta[k],momenta[l])*pairing
    return total


def check():
    momenta = ((1,0),(0,1),(1,0),(-2,-1))
    for kind in ("so","sp"):
        basis,sharp,commuting,jacobi = colours(kind)
        x,y = jacobi
        assert np.array_equal(comm(comm(x,y),x), y if kind=="so" else -4*y)
        n = len(x)
        restricted = bare([units(n,i,j) for i,j in product(range(n),repeat=2)],commuting,momenta)
        projected = bare(basis,commuting,momenta)
        exchange = pole_exchange(sharp,commuting,momenta)
        assert exchange == F(16,n)
        assert sum(restricted) == n*n*exchange/4
        expected = (F(3,5),F(2,5)) if kind=="so" else (F(24,5),F(16,5))
        assert tuple(projected) == expected
        expected_full = (F(12),F(8)) if kind=="so" else (F(48,5),F(32,5))
        assert tuple(restricted) == expected_full
        print(kind,n,"pure propagator",projected,"sum",sum(projected),
              "full matrix propagator",restricted,"sum",sum(restricted),"pole bracket",exchange)


if __name__ == "__main__":
    check()
