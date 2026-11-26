"""First relative CE differential on scalar and vector-output cochains of weight det^-3."""
from itertools import combinations, permutations, product
from random import Random

import numpy as np
from sympy import Matrix

PRIME = 1000003


def orbit_columns():
    """The unmarked and marked functional-graph bases on five arguments."""
    perms = list(permutations(range(5)))
    maps = [f for f in product(range(5), repeat=5) if all(i != j for i,j in enumerate(f))]
    columns = {}
    for marks in [(-1,), range(5)]:
        unused = {(f,j) for f in maps for j in marks}
        index = 0
        while unused:
            f,j = min(unused)
            orbit = set()
            for p in perms:
                g = [None]*5
                for i in range(5):
                    g[p[i]] = p[f[i]]
                orbit.add((tuple(g),p[j] if j >= 0 else -1))
            for key in orbit:
                columns[key] = index
            unused -= orbit
            index += 1
        assert index == (13 if len(marks) == 1 else 52)
    return maps, columns


def data(seed):
    rng = Random(seed)
    momenta, vectors = [], []
    for _ in range(6):
        p = [rng.randint(-2,2) for _ in range(5)]
        if not any(p):
            p[0] = 1
        k = next(i for i,x in enumerate(p) if x)
        u = [rng.randint(-2,2) for _ in range(5)]
        c = sum(x*y for x,y in zip(p,u))
        u = [p[k]*x for x in u]
        u[k] -= c
        assert sum(x*y for x,y in zip(p,u)) == 0
        momenta.append(p)
        vectors.append(u)
    return momenta,vectors,[rng.randint(-2,2) for _ in range(5)]


def binary_factor(poly, tcoef, scoef):
    """Multiply homogeneous binary form by tcoef*t+scoef*s (batch of evaluations)."""
    out = [tcoef*c % PRIME for c in poly] + [0]
    for k,c in enumerate(poly):
        out[k+1] = (out[k+1] + scoef*c) % PRIME
    return out


def det(rows):
    return int(Matrix(rows).det(method="domain-ge")) % PRIME


def batch(size):
    """Each row evaluates on six transverse vector fields u*(exp(p.z)-1-p.z)."""
    inputs = [data(k) for k in range(size)]
    maps, columns = orbit_columns()
    scalar = np.zeros((size,13),dtype=np.int64)
    current = np.zeros((size,52),dtype=np.int64)
    dot = np.array([[[sum(x*y for x,y in zip(us[i],ps[j])) for ps,us,_ in inputs]
                     for j in range(6)] for i in range(6)], dtype=np.int64)
    for i,j in combinations(range(6),2):
        rest = [k for k in range(6) if k not in (i,j)]
        sign = (-1)**(i+j)
        linear_dets = []
        for mark in range(-1,5):
            pair = []
            for source in (j,i):
                values = []
                for ps,_,gamma in inputs:
                    rows = [ps[k] for k in [source]+rest]
                    if mark >= 0:
                        rows[mark] = gamma
                    values.append(det(rows))
                pair.append(np.array(values,dtype=np.int64))
            linear_dets.append(pair)
        tj,si = linear_dets[0]
        det2 = [tj*tj%PRIME,2*tj*si%PRIME,si*si%PRIME]
        detpolys = [binary_factor(det2,*linear_dets[0])]
        for mark in range(5):
            v,w = linear_dets[mark+1]
            detpolys.append([c*v%PRIME for c in det2] if mark==0 else binary_factor(det2,v,w))
        # For given determinant and homogeneous input monomial s^k*t^(d-k)
        # retain the terms of total power >=2 in each of s,t.
        cut = {(mark,d,k): sum(c for m,c in enumerate(poly)
                                if k+m >= 2 and d-k+len(poly)-1-m >= 2)
               for mark,poly in enumerate(detpolys,-1) for d in range(1,6) for k in range(d+1)}
        for f in maps:
            target = rest[f[0]-1]
            parts = [dot[i,j]*dot[j,target]%PRIME, -dot[j,i]*dot[i,target]%PRIME]
            spectator = 1
            for s,to in enumerate(f[1:]):
                u = rest[s]
                if to:
                    spectator = spectator * dot[u,rest[to-1]] % PRIME
                else:
                    parts = binary_factor(parts,dot[u,j],dot[u,i])
            d = len(parts)-1
            for mark in range(-1,5):
                value = sum(c*cut[mark,d,k] % PRIME for k,c in enumerate(parts))
                dest = scalar if mark == -1 else current
                dest[:,columns[f,mark]] += sign*(value*spectator%PRIME)
    return scalar%PRIME,current%PRIME


def determinant(matrix):
    aa = matrix.copy()
    sign=1
    for k in range(len(aa)):
        nonzero=np.flatnonzero(aa[k:,k])
        if len(nonzero)==0:
            return 0
        pivot=k+int(nonzero[0])
        if k!=pivot:
            aa[[k,pivot]]=aa[[pivot,k]]
            sign=-sign
        factors=aa[k+1:,k]*pow(int(aa[k,k]),-1,PRIME)%PRIME
        aa[k+1:]=(aa[k+1:]-factors[:,None]*aa[k])%PRIME
    answer=sign
    for d in aa.diagonal():
        answer=answer*int(d)%PRIME
    return answer


def independent_rows(matrix):
    pivots={}
    selected=[]
    for i,vals in enumerate(matrix):
        row=vals.copy()
        for col,pivot in sorted(pivots.items()):
            row=(row-row[col]*pivot)%PRIME
        support=np.flatnonzero(row)
        if len(support):
            j=int(support[0])
            pivots[j]=row*pow(int(row[j]),-1,PRIME)%PRIME
            selected.append(i)
    return selected,sorted(pivots)


if __name__=='__main__':
    scalar,current=batch(72)
    for matrix,expected in [(scalar,736898),(current,746682)]:
        rows,cols=independent_rows(matrix)
        value=determinant(matrix[np.ix_(rows,cols)])
        assert len(rows)==matrix.shape[1] and value==expected
        print(f'prime {PRIME}, rank {len(rows)}, determinant {value}, rows {rows}')
