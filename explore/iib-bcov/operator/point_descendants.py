"""Point classes obtained by differentiating the minimal BCOV action in a scalar direction."""
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

from flint import fmpz_mat
from integrated_descendants import translated, derivative, multiply
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT))
from bulk_eight import BULK
from elementary_scalar import elementary_three_loop
from marked_boundaries import quartics


def check_classical(order=10):
    # The field series have two spare orders before both derivatives.
    potential = translated(0,order+2)
    for r in range(7):
        current = translated(r,order+2)
        successor = translated(r+1,order+2)
        first = derivative(current,0)
        for i in range(5):
            gradient = derivative(potential,i)
            if i<4:
                difference = Counter(multiply(first,gradient,order))
                difference.subtract(derivative(current,i+1))
                difference.subtract(derivative(successor,i))
                assert not any(c for powers,c in difference.items() if sum(powers)<=order)
            for j in range(5):
                difference = Counter(multiply(first,derivative(gradient,j),order))
                difference.subtract(derivative(derivative(successor,i),j))
                assert not any(c for powers,c in difference.items() if sum(powers)<=order)


def with_root(waves):
    return [*waves,tuple(-sum(component) for component in zip(*waves))]


def check_minimal_ranks():
    for r in range(5):
        for k in range(6):
            rank = 5-k
            if k<=r:
                # Complementary field in the action trace.
                assert 0<=k<=4-(r-k) and rank>=r-k+1
            else:
                level = k-r-1
                assert 0<=level<=4 and rank+level==4-r


def check_quantum():
    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    # On the first slot take the transverse four-vector e2 e3 e4 e5,
    # and on the last external slot take e1. The other fields are scalar.
    polyvectors = [*frame[1:],frame[0]]
    assert fmpz_mat(polyvectors).det()==1
    assert all(not sum(x*y for x,y in zip(frame[0],v)) for v in polyvectors[:4])
    assert not sum(x*y for x,y in zip(frame[-1],polyvectors[-1]))
    one = Fraction(int(fmpz_mat(frame).det())**2,24)
    two = -sum(c*q for c,q in zip(BULK,quartics(with_root([frame[0],*frame]))))
    waves = [frame[0]]*3+[frame[1]]*2+frame[2:]
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words']))==DOUBLE_SOFT
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-O3','-std=c++17',str(ROOT/'operator'/'nine_sextics.cpp'),
             '-o',str(exe)],check=True)
        row, = orbit_integers(exe,[with_root(waves)])
        terminal = elementary_three_loop(waves[:-1],waves[-1],frame[0],exe=exe)
    three = tuple(-sum(c*Fraction(x) for c,x in zip(row,column))
                  for column in zip(*data['rational_and_zeta']))
    assert one==Fraction(1,24) and two==Fraction(1,31104)
    assert three==(Fraction(49188779,2715939072000),Fraction(-7783,523908000))
    # Total external momentum evaluated on e1 is +3. The reduced
    # trace correction is -C, while the elementary differential is (-3) C.
    assert tuple(3*x for x in three)==terminal
    print('Point insertions, reduced periods at genera one, two and three:',one,two,*three)


if __name__=='__main__':
    check_minimal_ranks()
    check_classical()
    check_quantum()
