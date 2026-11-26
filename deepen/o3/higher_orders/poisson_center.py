"""Exact rank bounds for the scalar Poisson CE consistency condition on four central waves."""
from itertools import combinations
from math import prod
from random import Random
import sys


def bracket(p, q):
    return p[0]*q[1] - p[1]*q[0]


def invariants(p):
    u = []
    for i in range(4):
        q, r, s = [p[j] for j in range(4) if j != i]
        u.append((-1)**i * (bracket(q, r) + bracket(r, s) + bracket(s, q)))
    assert sum(u) == 0
    elementary = tuple(sum(prod(xs) for xs in combinations(u, j)) for j in (2,3,4))
    vandermonde = prod(x-y for x, y in combinations(u, 2))
    return elementary, vandermonde


def columns(w):
    degree = w if w%2 else w-6
    return [(i,j,k) for k in range(degree//4+1) for j in range((degree-4*k)//3+1)
            for i in [(degree-4*k-3*j)//2] if 2*i+3*j+4*k == degree]


def merge_terms(p):
    for i,j in combinations(range(len(p)), 2):
        merged = tuple(x+y for x,y in zip(p[i], p[j]))
        yield ((-1)**(i+j)*bracket(p[i],p[j]),
               [merged] + [a for t,a in enumerate(p) if t not in (i,j)])


def incoming(w, p):
    return sum(sign*bracket(q[0],q[1])**(w-1) for sign,q in merge_terms(p))


def row(w, p, triples, prime):
    values = [0]*len(triples)
    for sign, q in merge_terms(p):
        (x,y,z), v = invariants(q)
        if w%2 == 0:
            sign *= v
        for t,(i,j,k) in enumerate(triples):
            values[t] += sign * pow(x,i,prime)*pow(y,j,prime)*pow(z,k,prime)
    return [a % prime for a in values]


def check(w):
    triples = columns(w)
    prime = 1000003
    # An integer minor nonzero modulo p certifies the lower bound over Q.
    # The Jacobi identity supplies the upper bound from the preceding cubic.
    probe = [(1,0), (0,1), (-3,2), (2,-3)]
    incoming_rank = int(w%2 == 0 and incoming(w, probe) != 0)
    target = len(triples)-incoming_rank
    basis = {}
    rng = Random(317)
    for _ in range(3*len(triples)+30):
        if len(basis) == target:
            break
        points = [(rng.randrange(-40,41), rng.randrange(-40,41)) for _ in range(4)]
        p = points + [tuple(-sum(x) for x in zip(*points))]
        v = row(w, p, triples, prime)
        for k in sorted(basis):
            if v[k]:
                mult = v[k]
                v = [(x - mult*y) % prime for x,y in zip(v,basis[k])]
        pivot = next((k for k,x in enumerate(v) if x), None)
        if pivot is not None:
            basis[pivot] = [x*pow(v[pivot],-1,prime)%prime for x in v]
    return len(triples), incoming_rank, len(basis)


if __name__ == '__main__':
    end = int(sys.argv[1]) if len(sys.argv)>1 else 32
    for w in range(2, end+1):
        count, incoming_rank, lower = check(w)
        print(w, count, incoming_rank, lower, count-incoming_rank-lower, flush=True)
