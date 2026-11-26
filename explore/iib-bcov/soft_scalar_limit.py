"""Leading one-, two- and three-soft limits of the elementary nine-primary trace."""
from fractions import Fraction
from math import lcm
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from flint import fmpq_mat
from elementary_scalar import ROOT
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers


def soft_limit(exe,hard,velocities):
    k = len(velocities)
    if k not in (1,2,3) or len(hard)+k!=9:
        raise ValueError("Expected nine slots and one, two or three soft velocities.")
    points = [tuple(map(Fraction,p)) for p in (*hard,*velocities)]
    if any(len(p)!=5 for p in points) or any(sum(column) for column in zip(*points[:-k])):
        raise ValueError("Expected five-component covectors and conserved hard momenta.")
    scale = lcm(*(x.denominator for p in points for x in p))
    integers = [tuple(int(scale*x) for x in p) for p in points]
    row, = orbit_integers(exe,[integers],leading_soft=k)
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words']))==DOUBLE_SOFT
    return tuple(sum(v*Fraction(x) for v,x in zip(row,column))/scale**30
                 for column in zip(*data['rational_and_zeta']))


def check(exe):
    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    expected = [('-740466719/70614415872000','469153/54486432000'),
                ('-103647659/3851695411200','164077/7429968000'),
                ('-841479311/35307207936000','53359/2724321600')]
    vandermonde = fmpq_mat([[t**j for j in range(7)] for t in range(-3,4)])
    for k in (1,2,3):
        first = [*frame,*[frame[0]]*(3-k)]
        last = tuple(-sum(column) for column in zip(*first))
        hard = [*first,last]
        for direction,target in [(frame[0],('0','0')),((1,1,0,0,0),expected[k-1])]:
            velocities = [direction]*k
            pair = soft_limit(exe,hard,velocities)
            assert pair==tuple(map(Fraction,target))
            # Only one moving direction, so each of the six minors is affine in t.
            curves = []
            for t in range(-3,4):
                moved = tuple(p-k*t*x for p,x in zip(last,direction))
                curves.append([*first,moved,*[tuple(t*x for x in direction)]*k])
            polynomial = vandermonde.solve(fmpq_mat(orbit_integers(exe,curves))).tolist()
            assert all(not x for row in polynomial[:2*k] for x in row)
            direct, = orbit_integers(exe,[[*hard,*velocities]],leading_soft=k)
            assert polynomial[2*k]==direct
    print('All eighteen coefficients verified at soft orders two, four and six.')


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        check(exe)
