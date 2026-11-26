"""Constant and linear soft bivector consistency on twelve primary potentials."""
from fractions import Fraction
from itertools import combinations
from math import comb
from random import Random

from flint import fmpz_mat
from soft_consistency import interior, vector, wedge


SLOTS = tuple(range(12))
FIVE = tuple(combinations(SLOTS, 5))
SIX = tuple(combinations(SLOTS, 6))
COMPLEMENT = {s: tuple(i for i in SLOTS if i not in s) for s in SIX}
SOFT = (1, 0, 0, 0, 0)
BIVECTOR = ((6, 1),)  # e_2 wedge e_3
SCALAR = ((0, 1),)


def fixture():
    rng = Random(90132)
    p = [tuple(rng.randint(-2, 2) for _ in range(5)) for _ in range(11)]
    return p + [tuple(-sum(x) for x in zip(*p))]


def bracket(p, q):
    return p[1]*q[2] - p[2]*q[1]


def array(p, derivative=False):
    # Each six-slot block contains two distinct squared five-row minors.
    minors = {}
    jets = {}
    for subset in FIVE:
        rows = [p[i] for i in subset]
        value = int(fmpz_mat(rows).det())
        minors[subset] = value**2
        if derivative:
            jets[subset] = {i: 2*value*int(fmpz_mat(rows[:j]+[SOFT]+rows[j+1:]).det())
                            for j, i in enumerate(subset)}
    omitted = {s: {i: tuple(j for j in s if j != i) for i in s} for s in SIX}
    block = {s: sum(minors[omitted[s][i]]*minors[omitted[s][j]]
                    for i, j in combinations(s, 2)) for s in SIX}
    weights = [[0]*12 for _ in SLOTS]
    jet = 0
    for s in SIX:
        rest = block[COMPLEMENT[s]]
        for i, j in combinations(s, 2):
            a, b = omitted[s][i], omitted[s][j]
            weights[i][j] += minors[a]*minors[b]*rest
            if derivative:
                # The derivatives common to both slots, including the balancing
                # derivative, cancel in the antisymmetric pair sum.
                jet += bracket(p[i], p[j])*rest*(minors[a]*jets[b][i]-minors[b]*jets[a][j])
    result = [tuple(sum(p[j][r]*weights[min(i,j)][max(i,j)] for j in SLOTS)
                    for r in range(5)) for i in SLOTS]
    return result, jet


def polarized(p):
    rng = Random(315)
    while True:
        vectors = []
        for row in p[:5]:
            axis = next(i for i, x in enumerate(row) if x)
            v = [rng.randint(-3,3) for _ in range(5)]
            v[axis] = 0
            v[axis] = -sum(a*b for a,b in zip(row,v))
            vectors.append(tuple(x if j==axis else x*row[axis] for j,x in enumerate(v)))
        volume = int(fmpz_mat(vectors).det())
        if volume:
            return [vector(v) for v in vectors] + [SCALAR]*7, volume


def cochain(a, fields):
    total = 0
    for i, covector in enumerate(a):
        term = interior(covector, fields[i])
        if not term:
            continue
        value = SCALAR
        for j, f in enumerate(fields):
            value = wedge(value, term if i==j else f)
        parity = sum(fields[j][0][0].bit_count() for j in range(i))
        total += (-1)**parity*dict(value).get(31,0)
    return total


def hard_soft_fusions(p, fields, t):
    # The omitted hard-hard fusions have at least two powers of t.
    hard = [*p[:-1], (p[-1][0]-t, *p[-1][1:])]
    total = 0
    for i in SLOTS:
        fused = list(hard)
        fused[i] = (fused[i][0]+t, *fused[i][1:])
        field = interior(fused[i], wedge(fields[i], BIVECTOR))
        if field:
            values, _ = array(fused)
            parity = sum(fields[j][0][0].bit_count() for j in range(i))
            total += (-1)**parity*cochain(values, fields[:i]+[field]+fields[i+1:])
    return total


def check():
    p = fixture()
    values, jet = array(p, derivative=True)
    assert all(sum(row[a]*value[b]-row[b]*value[a] for row,value in zip(p,values))==0
               for a,b in combinations(range(5),2))
    assert jet == -8729370490059451776 and jet % 1009 == 782
    fields, volume = polarized(p)
    # Eight determinants, one free covector and a cubic bracket have t-degree at
    # most ten, since every determinant changes by a rank-one row update.
    at_zero = hard_soft_fusions(p,fields,0)
    assert at_zero == 0
    first = sum(Fraction((-1)**(t-1)*comb(10,t),t)
                * (hard_soft_fusions(p,fields,t)-at_zero) for t in range(1,11))
    assert first == volume*jet
    print('Zero bivector projection vanishes. First jet:', jet, 'mod 1009:', jet % 1009)
    print('Direct transverse CE derivative:', first, 'with input volume', volume)


if __name__=='__main__':
    check()
