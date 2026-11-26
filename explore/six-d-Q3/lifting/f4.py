"""Two-scalar F4 restriction and the full-adjoint Hessian contraction."""
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod
from collections import Counter
from random import Random

from quadratic import exact_rank


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def orbit_vectors():
    units = [tuple(s * 2 * (i == j) for i in range(4))
             for j in range(4) for s in (-1, 1)]
    short = units + list(product((-1, 1), repeat=4))
    long = []
    for i, j in combinations(range(4), 2):
        for a, b in product((-1, 1), repeat=2):
            long.append(tuple(a * (k == i) + b * (k == j) for k in range(4)))
    return short, short + long


SHORT, ROOTS = orbit_vectors()


def group_matrices():
    # Every row is twice a unit vector in the short-root orbit.
    def rows(chosen):
        if len(chosen) == 4:
            yield tuple(chosen)
            return
        for u in SHORT:
            if all(dot(u, v) == 0 for v in chosen):
                yield from rows([*chosen, u])
    # A signed Hadamard matrix divided by two takes all half-sign roots
    # to half-sign or unit roots; signed permutation matrices do likewise.
    return rows([])


def molien(max_degree):
    # Newton sums for Tr(Sym^d(w)), with w represented as 2w throughout.
    profiles = {}
    for matrix in group_matrices():
        powers = matrix
        traces = []
        for j in range(1, max_degree + 1):
            traces.append(sum(powers[i][i] for i in range(4)) // (2**j))
            assert sum(powers[i][i] for i in range(4)) % (2**j) == 0
            cols = tuple(zip(*matrix))
            powers = tuple(tuple(dot(row, col) for col in cols) for row in powers)
        profiles[tuple(traces)] = profiles.get(tuple(traces), 0) + 1
    assert sum(profiles.values()) == 1152
    answer = {(p,q): 0 for p in range(max_degree+1) for q in range(max_degree+1)}
    for traces, count in profiles.items():
        sym = [1]
        for d in range(1, max_degree + 1):
            numerator = sum(traces[j-1] * sym[d-j] for j in range(1,d+1))
            assert numerator % d == 0
            sym.append(numerator // d)
        for p, q in answer:
            answer[p,q] += count * sym[p] * sym[q]
    assert all(v % 1152 == 0 for v in answer.values())
    return {key: v // 1152 for key,v in answer.items()}


def polarized(d, k, x, y):
    return comb(d, k) * sum(dot(a,x)**(d-k) * dot(a,y)**k for a in SHORT)


def hessian(x):
    gradient = tuple(6 * sum(a[i] * dot(a,x)**5 for a in SHORT) for i in range(4))
    diagonal = tuple(30 * sum(a[i] * a[j] * dot(a,x)**4 for a in SHORT)
                     for i,j in product(range(4), repeat=2))
    roots = tuple(Fraction(dot(a,gradient), dot(a,x)) for a in ROOTS)
    return diagonal, roots


def triple_invariants():
    # The Weyl numerator applied to the third symmetric power of the adjoint.
    # All weights and rho use doubled orthonormal coordinates.
    roots = SHORT + [tuple(2*x for x in a) for a in ROOTS[len(SHORT):]]
    positive = [a for a in roots if dot(a,(9,4,2,1)) > 0]
    rho = tuple(sum(a[i] for a in positive)//2 for i in range(4))
    simple = [(0,2,-2,0), (0,0,2,-2), (0,0,0,2), (1,-1,-1,-1)]
    zero = (0,0,0,0)
    def add(a,b):
        return tuple(i+j for i,j in zip(a,b))
    def sub(a,b):
        return tuple(i-j for i,j in zip(a,b))
    def act(matrix,a):
        return tuple(dot(row,a)//2 for row in matrix)
    def determinant(matrix):
        if len(matrix)==1:
            return matrix[0][0]
        return sum((-1)**j*x*determinant(tuple(row[:j]+row[j+1:] for row in matrix[1:]))
                   for j,x in enumerate(matrix[0]) if x)
    group = [(matrix,determinant(matrix)//16) for matrix in group_matrices()]
    assert all(abs(s)==1 for matrix,s in group)
    sym = [Counter({zero:1}),Counter(),Counter(),Counter()]
    for a in roots+[zero]*4:
        for d in range(1,4):
            sym[d].update({add(k,a):v for k,v in sym[d-1].items()})
    char = sym[3]
    delta = [(sub(rho,act(matrix,rho)),s) for matrix,s in group]
    summands = {}
    for k in char:
        if min(dot(k,a) for a in simple)<0:
            continue
        multiplicity = sum(s*char[add(k,v)] for v,s in delta)
        if multiplicity:
            summands[k]=multiplicity
    answer = 0
    for left,m in summands.items():
        for right,n in summands.items():
            answer += m*n*sum(s*char[sub(act(matrix,add(right,rho)),add(left,rho))]
                              for matrix,s in group)
    assert answer == 107
    return answer


def wedge_dimensions():
    true_roots = [tuple(Fraction(v,2) for v in root) for root in SHORT] + ROOTS[len(SHORT):]
    positive = [root for root in true_roots if dot(root,(9,4,2,1)) > 0]
    rho = tuple(sum(root[i] for root in positive) / 2 for i in range(4))
    def dimension(weight):
        shifted = tuple(a+b for a,b in zip(rho,weight))
        return prod(dot(root,shifted) / dot(root,rho) for root in positive)
    assert dimension((1,1,0,0)) == 52
    assert dimension((2,1,1,0)) == comb(52,2) - 52 == 1274


def checks():
    wedge_dimensions()
    dimensions = molien(8)
    assert dimensions[8,0] == 3
    rng = Random(2871)
    def point():
        while True:
            x = tuple(rng.randint(-10,10) for _ in range(4))
            if all(dot(a,x) for a in ROOTS):
                return x
    samples = [(point(),point()) for _ in range(23)]
    products = {(0,0): [[1] * len(samples)]}
    for d in (2,6,8):
        for k in range(d+1):
            dp, dq = d-k, k
            if max(dp,dq) > 4:
                continue
            values = [polarized(d,k,x,y) for x,y in samples]
            for p,q in product(range(5-dp),range(5-dq)):
                for old in products.get((p,q), ()):
                    products.setdefault((p+dp,q+dq), []).append([a*b for a,b in zip(old,values)])
    def rank(rows):
        return exact_rank([dict(enumerate(v)) for v in rows])
    for p,q in product(range(5),repeat=2):
        assert rank(products.get((p,q), [])) == dimensions[p,q] - ((p,q) == (4,4))
    rows = products[4,4]
    cartan, full = [], []
    for x,y in samples:
        x0, xr = hessian(x)
        y0, yr = hessian(y)
        cartan.append(dot(x0,y0))
        full.append(dot(x0,y0) + dot(xr,yr))
    assert rank([*rows,cartan]) == rank([*rows,full]) == 8
    return len(rows), rank(rows), dimensions[4,4]


if __name__ == '__main__':
    print('F4 bidegree (4,4) polarized rows, their rank, full dimension:', checks())
