"""Exact interpolation of the unsplit-source C² three-loop kernel at zero displacement."""
from itertools import combinations, combinations_with_replacement
from math import prod
from random import Random
import sympy as s

from three_loop_current import build, integrate, SOURCES, PAIRED_SOURCES

PAIRS = tuple(combinations(range(4),2))


def basis():
    result = []
    for factors in combinations_with_replacement(range(6),3):
        # Straighten [02][13] with the Pluecker relation.
        if 1 in factors and 4 in factors:
            continue
        indices = [i for pair in factors for i in PAIRS[pair]]
        # Graph symmetry is odd under (12), even under (03). A source at
        # the bivalent vertex is needed for its Schwinger split differential.
        if indices.count(0)%2 == 1 and indices.count(1)%2 == 0 and 3 in indices:
            result.append(factors)
    assert len(result) == 11
    return result


def row(sources,monomials):
    directions = []
    for i,j,sign in ((1,2,-1),(0,3,-1),(0,3,1)):
        directions.append(tuple(x+sign*y for x,y in zip(sources[i],sources[j])))
    directions.append(sources[4])
    brackets = [directions[i][0]*directions[j][1]-directions[i][1]*directions[j][0] for i,j in PAIRS]
    return [prod(brackets[i] for i in factors) for factors in monomials]


def main():
    rng = Random(16795)
    monomials = basis()
    rows,values = [],[]
    while len(rows)<len(monomials):
        sources = [(rng.randint(-2,2),rng.randint(-2,2)) for _ in range(4)]
        sources.append(tuple(-sum(entries) for entries in zip(*sources)))
        vector = row(sources,monomials)
        if s.Matrix(rows+[vector]).rank() == len(rows):
            continue
        det,num,sign = build(0,sources)
        value = integrate(det,num,sign,[])
        assert set(value) <= {0}
        rows.append(vector)
        values.append(value.get(0,0))
    matrix = s.Matrix(rows)
    assert matrix.det() != 0
    coefficients = matrix.inv()*s.Matrix(values)
    assert {key:value for key,value in zip(monomials,coefficients) if value} == {
        (0,0,2):-s.Rational(1,768),(2,3,3):s.Rational(1,384),(2,4,4):-s.Rational(1,2304)}
    for sources,expected in ((SOURCES,-s.Rational(71,12)),(PAIRED_SOURCES,0)):
        assert (s.Matrix([row(sources,monomials)])*coefficients)[0] == expected
    print('determinant',matrix.det(),'coefficients',list(zip(monomials,coefficients)),flush=True)


if __name__ == '__main__':
    main()
