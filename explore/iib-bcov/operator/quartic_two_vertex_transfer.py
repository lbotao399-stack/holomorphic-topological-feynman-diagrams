"""Momentum-transfer tests of the physical {I_2,I_2}/2 primary tree."""
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import sys

import numpy as np
from flint import fmpz_mat

from fusion_transfer_box import EDGES, box, frame, jet_product, shift

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bulk_eight import BULK, primitives
from marked_boundaries import quartics


FIVE = list(combinations(range(12),5))
SIX = list(combinations(range(12),6))
INDEX = {s:i for i,s in enumerate(FIVE)}
OMIT = [[INDEX[s[:j]+s[j+1:]] for j in range(6)] for s in SIX]
SIGNS = np.array([1,-1,1,-1,1,-1],dtype=object)
DEN = 2987228160
assert BULK == (Fraction(-31,DEN), Fraction(42,DEN))


def minors(p):
    return [int(fmpz_mat([p[i] for i in s]).det()) for s in FIVE]


def numbers(p):
    # Each unordered six|six partition separating 0 and 1 appears once.
    return [-sum(p[j][4] for j in s) if 0 in s and 1 not in s else 0 for s in SIX]


def trace_numerators(signed, multiply):
    # The seventh cofactor (the internal wave) is zero in this basis.
    center = signed.sum(axis=-1,keepdims=True)
    six = 7*signed-center
    squares = multiply(six,six)
    last = multiply(center,center)
    norm = squares.sum(axis=-1,keepdims=True)+last
    fourth = multiply(squares,squares).sum(axis=-1,keepdims=True)+multiply(last,last)
    return (-31*multiply(norm,norm)+42*fourth)[...,0]


def measure(p):
    c = np.array(minors(p),dtype=object)[OMIT]*SIGNS
    values = trace_numerators(c, lambda a,b:a*b)
    return Fraction(int(np.dot(numbers(p),values*values[::-1])),DEN**2)


def literal(p):
    # Directly evaluate the two conserved seven-slot symbols on either side.
    value = Fraction()
    for s,weight in zip(SIX,numbers(p)):
        if not weight:
            continue
        factors = []
        for group in (s, tuple(j for j in range(12) if j not in s)):
            six = tuple(p[j] for j in group)
            recoil = tuple(-sum(c) for c in zip(*six))
            factors.append(sum(z*x for z,x in zip(BULK,quartics(six+(recoil,)))))
        value += weight*factors[0]*factors[1]
    return value


def derivative(p,directions):
    chosen = (EDGES[0], EDGES[1], EDGES[3])
    tables = []
    for mask in range(8):
        point = p
        for bit,((i,j),q) in enumerate(zip(chosen,directions)):
            if mask & (1<<bit):
                point = shift(point,i,j,q)
        tables.append(minors(point))
    jets = np.array(tables,dtype=object)
    # Each determinant is affine in each individual rank-one row update.
    for bit in (1,2,4):
        for mask in range(8):
            if mask & bit:
                jets[mask] -= jets[mask^bit]
    symbols = trace_numerators(jets[:,OMIT]*SIGNS,jet_product)
    product = jet_product(symbols,symbols[:,::-1])[7]
    return Fraction(int(np.dot(numbers(p),product)), DEN**2)


def cubic_orientation(p):
    waves = list(p[:5]) + list(p[9:11])
    waves.append(tuple(-sum(v) for v in zip(*waves)))
    actual = primitives(waves)
    expected = sum(c*(a[0][4]-a[1][4]) for c,a in zip(BULK,actual))
    test = Fraction()
    for i,j in combinations(range(8),2):
        sign = int(0 in (i,j)) - int(1 in (i,j))
        fused = tuple(x+y for x,y in zip(waves[i],waves[j]))
        others = tuple(waves[k] for k in range(8) if k not in (i,j))
        symbol = sum(c*b for c,b in zip(BULK,quartics(others+(fused,))))
        test += sign * fused[4] * symbol
    assert test == expected and expected


def check():
    p = frame()
    cubic_orientation(p)
    assert measure(p) == literal(p) == Fraction(7698906551931026804317,403107840)
    axes = tuple(tuple(int(i==j) for j in range(5)) for i in range(3))
    difference = box(measure,p,(axes[0],)*6)
    alternation = sum((-1)**sum(i>j for i,j in combinations(order,2))
                      * derivative(p,[axes[j] for j in order]) for order in permutations(range(3)))
    print("reference projection",measure(p),"sixth difference",difference,"alternating third derivative",alternation)
    assert difference == Fraction(-219328932226387589,2519424)
    assert alternation == Fraction(1460732123756158891,466560)


if __name__ == "__main__":
    check()
