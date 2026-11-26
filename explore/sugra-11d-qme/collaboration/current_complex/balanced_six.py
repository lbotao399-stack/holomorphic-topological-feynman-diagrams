"""Degree-six scalar divergence cokernel and its CE differential on balanced waves."""
from itertools import permutations, product
from random import Random

import numpy as np

PRIME = 251
PERM = np.array(list(permutations(range(6))), dtype=np.int16)
SLOTS = np.arange(6)


def orbits():
    # The six oriented five-minors coincide when the momenta sum to zero.
    # Their cube already contributes the alternating sign. Pairing monomials
    # therefore use a symmetric orbit of functions with no fixed inputs.
    unseen = set(product(*(tuple(j for j in range(6) if j != i) for i in range(6))))
    groups = []
    while unseen:
        f = min(unseen)
        targets = PERM[:,f]
        terms = np.empty_like(targets)
        np.put_along_axis(terms, PERM, targets, axis=1)
        orbit = sorted(set(map(tuple,terms)))
        unseen.difference_update(orbit)
        groups.append(np.array(orbit,dtype=np.int16))
    return groups


def det_mod(matrix):
    a = matrix.copy()%PRIME
    d = 1
    for i in range(len(a)):
        candidates = np.flatnonzero(a[i:,i])
        if not len(candidates):
            return 0
        k = i+int(candidates[0])
        if k != i:
            a[[i,k]] = a[[k,i]]
            d = -d
        pivot = int(a[i,i])
        d = d*pivot%PRIME
        a[i+1:] = (a[i+1:]-a[i+1:,i,None]*pow(pivot,-1,PRIME)*a[i])%PRIME
    return d


def evaluate(ps, us, groups):
    assert np.all(ps.sum(axis=0)%PRIME == 0)
    pairing = us@ps.T%PRIME
    assert not np.any(np.diag(pairing))
    volume = pow(det_mod(ps[:5]),3,PRIME)
    if not volume:
        return np.zeros(len(groups),dtype=np.int64)
    return np.array([int(np.prod(pairing[SLOTS,terms],axis=1).sum())%PRIME*volume%PRIME
                     for terms in groups],dtype=np.int64)


def waves(n,rng):
    ps = np.array([[rng.randrange(PRIME) for j in range(5)] for i in range(n)],dtype=np.int64)
    ps[-1] = -np.sum(ps[:-1],axis=0)%PRIME
    us = np.array([[rng.randrange(PRIME) for j in range(5)] for i in range(n)],dtype=np.int64)
    for p,u in zip(ps,us):
        j = next(j for j in range(5) if p[j])
        u[j] -= int(u@p)*pow(int(p[j]),-1,PRIME)
    return ps,us%PRIME


def coboundary(ps,us,groups):
    answer = np.zeros(len(groups),dtype=np.int64)
    for i in range(7):
        for j in range(i+1,7):
            keep = [k for k in range(7) if k not in (i,j)]
            p = (ps[i]+ps[j])%PRIME
            u = ((us[i]@ps[j])*us[j]-(us[j]@ps[i])*us[i])%PRIME
            answer += (-1)**(i+j)*evaluate(np.vstack((p,ps[keep])),np.vstack((u,us[keep])),groups)
    return answer%PRIME


def check():
    groups = orbits()
    rng = Random(20261121)
    pivots = {}
    raw_rows, selected = [], []
    for n in range(20):
        ps,us = waves(7,rng)
        row = coboundary(ps,us,groups)
        raw_rows.append(row.copy())
        for i in sorted(pivots):
            row = (row-row[i]*pivots[i])%PRIME
        indices = np.flatnonzero(row)
        if len(indices):
            k = int(indices[0])
            pivots[k] = row*pow(int(row[k]),-1,PRIME)%PRIME
            selected.append(n)
    minor = det_mod(np.array(raw_rows)[np.ix_(selected,sorted(pivots))])
    assert len(groups)==40 and len(pivots)==9 and minor==198
    print("graph orbits",len(groups),"pivots",sorted(pivots),"rank",len(pivots),"minor",minor)


if __name__=='__main__':
    check()
