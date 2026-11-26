"""Exact image and cokernel of relative divergence from its flat Koszul complex."""
from itertools import combinations
from math import comb

import numpy as np

from relative_dimensions import bounds, exterior, nodes, torus_data


def counts(n, order, qmax, degree, determinant, prime, generator):
    factors, _, characters = torus_data(n,order,qmax,degree,determinant,prime,generator)
    coefficients = exterior(characters,degree,qmax,
                            lambda value,q:value*pow(q,-1,prime)%prime,
                            lambda value:value%prime)
    positions = np.array([p for p,size in nodes(n,order)])
    primitive = pow(generator,(prime-1)//order,prime)
    roots = np.array([pow(primitive,i,prime) for i in range(order)],dtype=np.int64)
    result = []
    for j in range(n+1):
        wedge = np.zeros(len(factors),dtype=np.int64)
        for indices in combinations(range(n),j):
            eigenvalue = roots[-positions[:,list(indices)].sum(axis=1,dtype=np.int64)%order]
            wedge = (wedge + eigenvalue) % prime
        measure = factors*wedge%prime
        result.append([int((measure*coefficients[q][degree-j]%prime).sum()%prime)
                       for q in range(qmax+1)])
    return result


def exact_dimensions():
    n,degree,qmax,order = 5,15,15,32
    determinant = -3
    ordinary = bounds(n,degree,qmax)
    upper = max(comb(n,j)*row[degree-j] for j in range(n+1) for row in ordinary)
    mod = 1
    answer = [[0]*(qmax+1) for j in range(n+1)]
    for prime,generator in [(998244353,3),(469762049,3),(167772161,3),(754974721,11),
                            (2013265921,31),(1811939329,13)]:
        values = counts(n,order,qmax,degree,determinant,prime,generator)
        inv = pow(mod,-1,prime)
        for j in range(n+1):
            for q,v in enumerate(values[j]):
                answer[j][q] += mod*((v-answer[j][q])*inv%prime)
        mod *= prime
        if mod > upper:
            break
    assert mod > upper
    for j,row in enumerate(answer):
        assert all(v<=comb(n,j)*ordinary[q][degree-j] for q,v in enumerate(row))
    return list(zip(*answer))


def main():
    rows = exact_dimensions()
    assert rows[8][:2] == (55495,109050)
    assert rows[9][:2] == (197563,302166)
    output = []
    print('q   scalar, current, wedge^2,...,wedge^5   div-rank   ker   coker')
    for q in range(5,16):
        row = rows[q]
        coker = sum((-1)**j*v for j,v in enumerate(row))
        image = row[0]-coker
        kernel = row[1]-image
        assert min(coker,image,kernel)>=0
        output.append((q,kernel,coker))
        print(q,*row,image,kernel,coker)
    reduced = output[9-5][1]+output[8-5][2]
    print('degree 9 dimension after the divergence differential',reduced)
    print('degree 9 relative divergence-free current space',output[9-5][1])


if __name__ == '__main__':
    main()
