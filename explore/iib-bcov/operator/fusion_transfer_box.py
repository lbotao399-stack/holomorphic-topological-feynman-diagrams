"""Split-independent fusion tests on waves in a hyperplane and its complement."""
from functools import cache
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
from random import Random
import sys

import numpy as np
from flint import fmpz_mat


def frame():
    rng = Random(4722)
    p = [tuple(rng.randrange(-3,4) for _ in range(4))+(0,) for _ in range(9)]
    p.append(tuple(rng.randrange(-3,4) for _ in range(4))+(1,))
    p.append(tuple(rng.randrange(-3,4) for _ in range(4))+(2,))
    p.append(tuple(-sum(column) for column in zip(*p)))
    return tuple(p)


def shift(p, i, j, q):
    result = list(p)
    result[i] = tuple(a+b for a,b in zip(p[i],q))
    result[j] = tuple(a-b for a,b in zip(p[j],q))
    return tuple(result)


@cache
def split_control(p):
    # Each six-slot block contributes a soft scalar of determinant weight four.
    five = {s: int(fmpz_mat([p[i] for i in s]).det())**2
            for s in combinations(range(12),5)}
    blocks = {}
    for s in combinations(range(12),6):
        omissions = [five[tuple(j for j in s if j!=i)] for i in s]
        blocks[s] = sum(a*b for a,b in combinations(omissions,2))
    value = 0
    for s,h in blocks.items():
        if 0 not in s or 1 in s:
            continue
        complement = tuple(j for j in range(12) if j not in s)
        value += -2*sum(p[j][4] for j in s)*h*blocks[complement]
    return value


@cache
def trace_control(p):
    # A single eleven-slot symmetric trace of degree forty.
    def fused(i,j):
        merged = [row for pos,row in enumerate(p) if pos not in (i,j)]
        merged.append(tuple(x+y for x,y in zip(p[i],p[j])))
        matrix = fmpz_mat(merged)
        return int((matrix.transpose()*matrix).det())**4
    return sum(p[k][4]*(fused(0,k)-fused(1,k)) for k in (9,10,11))


EDGES = ((0,9),(0,10),(0,11),(1,9),(1,10),(1,11))


def box(measure,p,directions):
    total = 0
    for subset in product(range(2),repeat=len(EDGES)):
        point = p
        for active,(i,j),q in zip(subset,EDGES,directions):
            if active:
                point = shift(point,i,j,q)
        total += (-1)**sum(subset)*measure(point)
    return total


def jet_product(a,b):
    # Three labeled nilpotent directions; array values remain exact integers.
    return np.array([sum(a[k]*b[mask^k] for k in range(mask+1) if not k&~mask)
                     for mask in range(8)],dtype=object)


def third_derivative(p,directions):
    edges = (EDGES[0],EDGES[1],EDGES[3])
    five = list(combinations(range(12),5))
    six = list(combinations(range(12),6))
    minors = []
    for mask in range(8):
        point = p
        for n,((i,j),q) in enumerate(zip(edges,directions)):
            if mask>>n&1:
                point = shift(point,i,j,q)
        minors.append([int(fmpz_mat([point[i] for i in s]).det()) for s in five])
    coeff = np.array(minors,dtype=object)
    # Each determinant is affine separately in each rank-one update.
    for bit in (1,2,4):
        for mask in range(8):
            if mask&bit:
                coeff[mask] -= coeff[mask^bit]
    squares = jet_product(coeff,coeff)
    subset_index = {s:i for i,s in enumerate(five)}
    omitted = [[subset_index[tuple(j for j in s if i!=j)] for i in s] for s in six]
    omit_squares = squares[:,omitted]
    total = omit_squares.sum(axis=2)
    block = (jet_product(total,total)-jet_product(omit_squares,omit_squares).sum(axis=2))//2
    # Lexicographic complement reverses the list of six-subsets.
    both = jet_product(block,block[:,::-1])[7]
    return sum(-2*sum(p[j][4] for j in s)*both[k] for k,s in enumerate(six)
               if 0 in s and 1 not in s)


def third_symbol():
    # On the physical three-outside slice the six transfer directions are X+tZ,
    # Y+tZ at three distinct t. All cubic minors of (X,Y,Z) annihilate each face.
    triples = list(combinations(range(4),3))
    for omit in range(6):
        rows = []
        rng = Random(921+omit)
        for _ in range(16):
            x,y,z = [[rng.randrange(-5,6) for _ in range(4)] for _ in range(3)]
            t = (0,3,2)[omit%3]
            if omit<3:
                x = [-t*a for a in z]
            else:
                y = [-t*a for a in z]
            rows.append([int(fmpz_mat([[v[j] for j in s] for v in (x,y,z)]).det())
                         for s in triples])
        assert not any(any(row) for row in rows)
    # Three independent edge blocks leave trilinear tensors in degree three.
    # The remaining faces impose skewness under the generating transpositions.
    words = list(product(range(4),repeat=3))
    matrix = []
    for word in words:
        for a,b in ((0,1),(0,2)):
            permuted = list(word)
            permuted[a],permuted[b] = permuted[b],permuted[a]
            row = [0]*len(words)
            row[words.index(word)] += 1
            row[words.index(tuple(permuted))] += 1
            matrix.append(row)
    assert len(words)-fmpz_mat(matrix).rank()==len(triples)==4


def check():
    p = frame()
    axes = tuple(tuple(int(j==i) for j in range(5)) for i in range(4))
    for directions in (axes+axes[:2],(axes[0],)*6):
        assert box(trace_control,p,directions)==0
    measured = box(split_control,p,(axes[0],)*6)
    print('Sixth difference of reciprocal first-jet control:',measured,'mod 1009:',measured%1009)
    assert measured == -140342433697054124544 and measured%1009==776
    third = sum((-1)**sum(i>j for i,j in combinations(perm,2))
                *third_derivative(p,[axes[i] for i in perm]) for perm in permutations(range(3)))
    print('Alternating physical third derivative:',third,'mod 1009:',third%1009)
    assert third == -4726500469225893262080 and third%1009==142
    third_symbol()
    print('All split trace differences vanish. Four independent cubic minors.')


def check_theta():
    sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
    from bulk_eight import SHAPES, potentials
    from genus_two_probes import MOMENTA, evaluate_potential

    # These fourteen rows are unisolvent even without the soft restrictions.
    columns = [evaluate_potential(potentials(MOMENTA,shapes=(shape,))) for shape in SHAPES]
    assert all(3*c==-a+2*b and 9*d==-20*a+13*b for a,b,c,d in zip(*columns))
    first,second = columns[:2]
    assert first[0]*second[1]-first[1]*second[0]

    original = frame()
    waves = [*original[:5],*original[9:11]]
    waves.append(tuple(-sum(col) for col in zip(*waves)))
    directions = [tuple(int(j==k) for j in range(5)) for k in (0,1,2,3,0,1)]
    result = [Fraction(0)]*4
    for subset in product(range(2),repeat=6):
        point = waves
        for active,(i,j),q in zip(subset,product((0,1),(5,6,7)),directions):
            if active:
                point = shift(point,i,j,q)
        for k,shape in enumerate(SHAPES):
            a = potentials(point,shapes=(shape,))
            result[k] += (-1)**sum(subset)*(a[0][4]-a[1][4])
    assert result == [Fraction(102820139,136080)*x for x in (2,1,0,-3)]
    print('Four integrated theta shapes:',*(str(x) for x in result),'sum',sum(result))


if __name__=='__main__':
    check()
    check_theta()
