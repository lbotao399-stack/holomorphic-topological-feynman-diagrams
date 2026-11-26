"""Leading normal variation at four collinear, nonzero momenta of full total rank."""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from flint import fmpq_mat
from elementary_scalar import ROOT
from worldsheet_primary import fourfold_primary
from three_loop_soft_jet import orbit_integers


def waves(strengths,variations,shears,parameter):
    first = [(a,*(parameter*v for v in variation)) for a,variation in zip(strengths,variations)]
    first.extend((shear,*(int(i==j) for i in range(4))) for j,shear in enumerate(shears))
    return [*first,tuple(-sum(column) for column in zip(*first))]


def normal_norm(strengths,variations):
    slopes = [tuple(F(x,a) for x in w) for a,w in zip(strengths,variations)]
    value = F(0)
    for left,right in combinations(slopes,2):
        difference = [x-y for x,y in zip(left,right)]
        value += sum(x*x for x in difference)-sum(difference)**2/5
    return prod(a*a for a in strengths)*value


def quadratic_row(exe,strengths,variations,shears):
    # Each minor has parameter degree at most the rank of all moving directions.
    degree = 6*fmpq_mat(variations).rank()
    grid = list(range(-degree//2,degree//2+1))
    points = [waves(strengths,variations,shears,t) for t in grid]
    values = fmpq_mat(orbit_integers(exe,points))
    powers = fmpq_mat([[t**i for i in range(degree+1)] for t in grid])
    rows = powers.solve(values).tolist()
    assert all(not x for row in rows[:2] for x in row)
    return rows[2]


def check():
    zero = (0,0,0,0)
    fixtures = [((1,1,1,1),(zero,zero,(1,0,0,0),(1,0,0,0)),(0,0,0,0)),
                ((1,-2,2,3),((1,0,0,0),(0,1,0,0),(1,-2,0,0),(-1,1,0,0)),(2,-1,1,0))]
    universal = [0]*18
    universal[2],universal[8] = -240,240
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for strengths,variations,shears in fixtures:
            row = quadratic_row(exe,strengths,variations,shears)
            norm = normal_norm(strengths,variations)
            assert [F(str(x)) for x in row]==[norm*x for x in universal]
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    a,b = (sum(value*F(x) for value,x in zip(universal,col)) for col in zip(*data['rational_and_zeta']))
    assert (a,b)==(F(-740466719,169474598092800),F(469153,130767436800))
    # The exact fourfold soft coefficient on the projected transverse simplex.
    # A single moving direction makes every four-row minor affine in t.
    grid = range(-3, 4)
    vandermonde = fmpq_mat([[t**j for j in range(7)] for t in grid])
    basis = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    for strengths, variations, _ in fixtures:
        total = [F(0), F(0)]
        slopes = [tuple(F(x, n) for x in row) for n, row in zip(strengths, variations)]
        for left, right in combinations(slopes, 2):
            delta = tuple(x-y for x, y in zip(left, right))
            values = []
            for t in grid:
                first = [*basis, tuple(t*x for x in delta)]
                last = tuple(-sum(col) for col in zip(*first))
                values.append([str(c) for c in fourfold_primary(3, [*first, last])])
            coefficient = vandermonde.solve(fmpq_mat(values)).tolist()[2]
            for i, c in enumerate(coefficient):
                total[i] += 2 * prod(x*x for x in strengths) * F(str(c))
        assert total == [-c*normal_norm(strengths, variations) for c in (a,b)]
    assert b>0 and a+b*F(2089,1728)<0
    print('Universal normal quadratic coefficient:',a,b)


if __name__=='__main__':
    check()
