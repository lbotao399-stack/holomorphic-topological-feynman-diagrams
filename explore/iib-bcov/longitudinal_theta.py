"""Recover the longitudinal part of the eight-slot two-loop theta potential."""
from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial, lcm

from flint import fmpq_mat, fmpz_mat
from bulk_eight import BULK, grounded_potentials, potentials, primitives, residual
from marked_boundaries import quartics


# Oriented paths and cycles on omitted labels. Each word has four edges.
WORDS = (
    ((0,1),(1,2),(2,3),(3,0)),
    ((0,1),(0,1),(2,3),(2,3)),
    ((0,1),(1,2),(2,3),(3,4)),
    ((0,1),(0,1),(2,3),(3,4)),
    ((0,1),(1,2),(3,4),(4,5)),
)


def complementary_minors(p):
    minors = [[0]*7 for _ in p]
    for i,j in combinations(range(7),2):
        value = (-1)**(i+j-1)*int(fmpz_mat([row for k,row in enumerate(p)
                                              if k!=i and k!=j]).det())
        minors[i][j] = value
        minors[j][i] = -value
    return minors


def orbit_values(p):
    minors = complementary_minors(p)
    values = []
    for word in WORDS:
        size = 1+max(j for pair in word for j in pair)
        total = 0
        for slots in permutations(range(7),size):
            a,b,c,d = (minors[slots[i]][slots[j]] for i,j in word)
            total += a*b*c*d
        values.append(total*factorial(7-size))
    return values


def longitudinal(p):
    """The scalar multiplying the distinguished momentum, on seven integer waves."""
    minors = complementary_minors(p)
    pairs = {frozenset((i,j)): minors[i][j]**2 for i,j in combinations(range(7),2)}
    triples = {}
    for i,j,k in combinations(range(7),3):
        a,b,c = minors[i][j],minors[j][k],minors[k][i]
        triples[frozenset((i,j,k))] = a*b+b*c+c*a
    first = sum(pairs[s]*pairs[t] for s,t in combinations(pairs,2) if s.isdisjoint(t))
    mixed = sum(pairs[s]*triples[t] for s in pairs for t in triples if s.isdisjoint(t))
    last = sum(triples[s]*triples[t] for s,t in combinations(triples,2) if s.isdisjoint(t))
    return F(first,4860) + F(mixed,10368) + F(last,17280)


def closed_potentials(p):
    """All eight grounded potentials for conserved rational fivefold momenta."""
    if len(p)!=8 or any(len(q)!=5 for q in p) or any(map(sum,zip(*p))):
        raise ValueError("Expected eight conserved five-component momenta.")
    scale = lcm(*(F(x).denominator for q in p for x in q))
    p = tuple(tuple(int(F(x)*scale) for x in q) for q in p)
    trace = primitives(p)
    factors = [longitudinal(p[:i]+p[i+1:]) for i in range(8)]
    common = [sum(w*wave[k] for w,wave in zip(factors,p))/2 for k in range(5)]
    return [tuple((sum(c*a[i][k] for c,a in zip(BULK,trace)) + q[k]*factors[i]+common[k])
                  /scale**21 for k in range(5)) for i,q in enumerate(p)]


def scalar_square(scalars,vector,polarization):
    """Fully polarized second-order point symbol in reduced units."""
    if len(scalars)!=6 or len(polarization)!=5:
        raise ValueError("Expected six scalars and a five-component polarization.")
    fields = (*scalars,vector)
    root = tuple(-sum(col) for col in zip(*fields))
    phi = closed_potentials((*fields,root))
    return sum((v-r)*u for v,r,u in zip(phi[6],phi[7],polarization))


def primary_coefficient(scalars,vector,polarization):
    """Coefficient of f', f'', or f''' for Q2 f(phi) on primary waves."""
    fields=(*scalars,vector)
    if len(scalars) not in (5,6,7) or any(len(q)!=5 for q in (*fields,polarization)):
        raise ValueError("Expected five to seven scalar waves and a fivefold vector wave.")
    if len(scalars)==6:
        return scalar_square(scalars,vector,polarization)
    scale=lcm(*(F(x).denominator for q in fields for x in q))
    p=tuple(tuple(int(F(x)*scale) for x in q) for q in fields)
    if len(scalars)==5:
        root=tuple(-sum(col) for col in zip(*p))
        trace=sum(c*x for c,x in zip(BULK,quartics((*p,root))))
        return trace*sum(x*v for x,v in zip(root,polarization))/scale**21
    phi=grounded_potentials(p)
    return 2*sum(x*v for x,v in zip(phi[-1],polarization))/scale**21


def sample(p):
    fields = residual(p,potentials(p))
    common = tuple(sum(col)/10 for col in zip(*fields))
    rows, values = [], []
    for i,wave in enumerate(p):
        pivot = next(k for k,v in enumerate(wave) if v)
        factor = (fields[i][pivot]-common[pivot])/wave[pivot]
        assert all(factor*x==y-c for x,y,c in zip(wave,fields[i],common))
        rows.append(orbit_values([row for j,row in enumerate(p) if i!=j]))
        values.append(factor)
    return rows, values


def fixture(a,b):
    p = [tuple(int(i==j) for j in range(5)) for i in range(5)]+[a,b]
    return tuple([*p,tuple(-sum(col) for col in zip(*p))])


if __name__ == "__main__":
    rows,values = [],[]
    for a,b in (((1,-1,1,0,-1),(0,1,1,-1,-1)), ((1,2,-1,3,2),(2,-1,2,1,-3))):
        r,v = sample(fixture(a,b))
        rows.extend(r)
        values.extend(v)
    matrix = fmpq_mat([[str(x) for x in row] for row in rows])
    assert matrix.rank()==3
    for word,y in zip(rows,values):
        assert word[1]==2*word[0] and word[3]==-2*word[2]
        assert F(word[1],233280)+F(word[3],82944)+F(word[4],138240)==y
    fixtures = [fixture((1,2,-1,3,2),(2,-1,2,1,-3)),
                fixture((2,3,-2,-1,4),(-2,0,1,4,-1)),
                fixture((0,0,0,0,0),(0,0,0,0,0))]
    for p in fixtures:
        assert closed_potentials(p)==potentials(p)
    p = fixtures[1]
    half = tuple(tuple(F(x,2) for x in q) for q in p)
    assert closed_potentials(half)==[tuple(x/2**21 for x in q) for q in potentials(p)]
    basis = [tuple(int(i==j) for j in range(5)) for i in range(5)]
    scalars = [basis[k] for k in (0,0,1,1,2,3)]
    assert scalar_square(scalars,basis[4],(1,1,1,1,0))==F(1,405)
    assert scalar_square(scalars,basis[4],(0,0,1,1,0))==F(1,1215)
    assert primary_coefficient([basis[i] for i in (0,0,1,2,3)],basis[4],(0,0,0,-1,0))==-F(1,31104)
    assert primary_coefficient([basis[i] for i in (0,0,1,1,2,2,3)],basis[4],(1,1,1,1,0))==F(7,720)
    print("Three independent longitudinal coefficients; all 40 potentials agree with theta integration on every fixture.")
