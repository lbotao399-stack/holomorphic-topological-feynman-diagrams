"""One invisible trace direction on the single-collinear-pair locus."""
from itertools import combinations
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

from flint import fmpq_mat, fmpz_mat

from three_loop_soft_jet import DOUBLE_SOFT, SUBSETS, native, orbit_integers, random_points


# Each of the 36 pairs is covered by one or more of the six five-row minors.
WORD = (10,22,27,41,104,108)


def collinear_points():
    rng = Random(89)
    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    points = []
    for _ in range(17):
        scale = rng.randint(1,8)
        first = [*frame,tuple(scale*x for x in frame[0])]
        first.extend(tuple(rng.randint(-3,3) for i in range(5)) for j in range(2))
        points.append([*first,tuple(-sum(column) for column in zip(*first))])
    return points


def certificate(exe):
    minors = [SUBSETS[j] for j in WORD]
    assert all(any(i in minor and j in minor for minor in minors)
               for i,j in combinations(range(9),2))
    assert all(sum(i in minor for minor in minors)>=2 for i in range(9))
    generic = orbit_integers(exe,random_points(18),(*DOUBLE_SOFT,WORD))
    matrix = fmpq_mat([r[:-1] for r in generic])
    rhs = fmpq_mat([[r[-1]] for r in generic])
    assert int(matrix.det())%1009==840 and rhs[0,0]==-98484480
    coordinates = matrix.solve(rhs)
    assert matrix*coordinates==rhs
    single = native(exe,collinear_points(),(*DOUBLE_SOFT,WORD),False)
    assert all(r[-1]==0 for r in single)
    minor = int(fmpz_mat([r[1:-1] for r in single]).det())%1009
    assert minor
    return minor,[str(x) for x in coordinates.entries()]


if __name__=='__main__':
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(exe)],check=True)
        minor,coordinates = certificate(exe)
    print('Single-pair rank 17, minor modulo 1009',minor)
    print('Unique invisible direction in the eighteen trace coordinates',*coordinates)
