"""Replay the joint-soft transfer identity and the exceptional six-slot kernel."""
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import prod
from random import Random

from flint import fmpq_mat, fmpz_mat


def plus(p,q,scale=1):
    return tuple(a+scale*b for a,b in zip(p,q))


def balanced(p):
    return [*p,tuple(-sum(col) for col in zip(*p))]


def fuse(p,i,j):
    return [q for k,q in enumerate(p) if k not in (i,j)]+[plus(p[i],p[j])]


SUBSETS = list(combinations(range(7),5))
COVERING = [(i,j) for i,j in combinations(range(len(SUBSETS)),2)
            if len(set(SUBSETS[i]+SUBSETS[j]))==7]


@cache
def seven(p):
    # A symmetric genus-two symbol of four determinants and seven double-soft legs.
    sq = [int(fmpz_mat([p[i] for i in s]).det())**2 for s in SUBSETS]
    return sum(sq[i]*sq[j] for i,j in COVERING)


def potential(p,slot,trace):
    weights = [(p[j],trace(tuple(fuse(p,slot,j)))) for j in range(len(p)) if j!=slot]
    return tuple(sum(q[i]*w for q,w in weights) for i in range(5))


@cache
def weights(degree,order):
    # A dual Vandermonde row extracts the chosen power exactly.
    matrix = fmpq_mat([[t**k for k in range(degree+1)] for t in range(degree+1)])
    return tuple(Fraction(str(x)) for x in matrix.inv().tolist()[order])


def soft_points(p,soft,times):
    waves = [tuple(t*x for x in q) for t,q in zip(times,soft)]
    return [*waves,*p[:-1],tuple(p[-1][i]-sum(q[i] for q in waves) for i in range(5))]


def extract(p,soft,degree,fn):
    row = weights(degree,2)
    result = None
    for times in product(range(degree+1),repeat=len(soft)):
        factor = prod(row[t] for t in times)
        if not factor:
            continue
        values = fn(soft_points(p,soft,times))
        if result is None:
            result = [Fraction(0)]*len(values)
        for i,value in enumerate(values):
            result[i] += factor*value
    return tuple(result)


def jet(p,soft):
    return extract(p,soft,4,lambda q:(seven(tuple(q)),))[0]


def check_transfer(p,soft):
    # Compare the full row difference, including its longitudinal terms.
    i,j = 0,1
    def difference(q):
        a,b = [potential(q,k+len(soft),seven) for k in (i,j)]
        return plus(a,b,-1)
    actual = extract(p,soft,5,difference)
    hard = [potential(p,k,lambda q:jet(q,soft)) for k in (i,j)]
    rhs = plus(hard[0],hard[1],-1)
    for s,q in enumerate(soft):
        other = soft[:s]+soft[s+1:]
        derivative = Fraction(0)
        for t,factor in enumerate(weights(4,1)):
            shifted = list(p)
            shifted[i] = plus(shifted[i],q,t)
            shifted[j] = plus(shifted[j],q,-t)
            derivative += factor*jet(shifted,other)
        rhs = plus(rhs,q,derivative)
    assert actual==rhs and any(actual)


def check_exception(p):
    def six(q):
        return int(fmpz_mat(q[:-1]).det())**2
    common = potential(p,0,six)
    assert any(common) and all(potential(p,j,six)==common for j in range(1,len(p)))


def length_section(n):
    # Two independent Witt pairs realize arbitrary lengths on conserved waves.
    a = balanced([(1,i) for i in range(1,n)])
    matrix = []
    for i,(ai,bi) in enumerate(a):
        row = [0]*(2*n)
        row[2*i:2*i+2] = [2*ai,2*bi]
        matrix.append(row)
    matrix += [[int(j%2==k) for j in range(2*n)] for k in range(2)]
    assert fmpq_mat(matrix).rank()==n+2


if __name__=='__main__':
    rng = Random(131904)
    def waves(n):
        return [tuple(rng.randint(-2,2) for _ in range(5)) for j in range(n)]
    for k in range(3):
        check_transfer(balanced(waves(6-k)),waves(k+1))
        seven.cache_clear()
    check_exception(balanced(waves(6)))
    for g in range(1,16):
        m = 2*g
        for k in range(m-1):
            n = m+3-k
            length_section(n)
            assert 5*m-2*k < (1<<(n-1))*(2*m-k)
    print('One-, two-, and three-soft transfer coefficients agree, including hard fusions.')
    print('Witt length sections and degree bounds passed. Six-slot exception has common nonzero rows.')
