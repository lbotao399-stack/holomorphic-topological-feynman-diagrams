"""Required individual bare four-loop potentials from the saturated fused trace."""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

from flint import fmpz_mat
from four_loop_planar import restricted_symbols
from quartic_seven_plane_hexagon import normal_component as hexagon_normal, planar_factors, restrict
from quartic_seven_plane_theta import NORMAL, normal_component as theta_normal, planar_factor
from quartic_scalar_response import BIPARTITE, total_channels, interval

CONTACT_CODE = Path(__file__).resolve().parents[2]/'twisted-holography-N4'/'code'
sys.path.insert(0,str(CONTACT_CODE))
from flux_maximal_minors import words

PLANES = (
    ((1,0),(1,0),(2,0),(0,2),(0,2),(1,2),(-2,-1)),
    ((2,0),(2,0),(0,-2),(0,1),(1,1),(2,2),(1,-1)),
)


def normal_boundary(normal, i, j):
    assert len(normal) == 5 and normal[i][2] == normal[j][2] == 0
    assert not any(map(sum,zip(*normal)))
    def minor(r,s):
        return int(fmpz_mat([row for k,row in enumerate(normal) if k not in (r,s)]).det())
    return sum(p[2]*(minor(i,s)**8-minor(j,s)**8)
               for s,p in enumerate(normal) if s not in (i,j))


def literal_boundary(plane,normal,lifts):
    """Fuse actual five-dimensional momenta, then sum both global determinant words."""
    all_rows = [(*r,0,0,0) for r in plane]+[(*xy,*q) for xy,q in zip(lifts,normal)]
    assert not any(map(sum,zip(*all_rows)))
    subsets, first, second = words()
    value = [0,0]
    for r,sign in ((7,1),(8,-1)):
        for s in (9,10,11):
            fused = [p for i,p in enumerate(all_rows) if i not in (r,s)]
            fused.append(tuple(x+y for x,y in zip(all_rows[r],all_rows[s])))
            minors = [int(fmpz_mat([fused[i] for i in ids]).det())**2 for ids in subsets]
            for ch,pool in enumerate((first,second)):
                trace = sum(minors[i]*minors[j]*minors[k]*minors[l] for i,j,k,l in pool)
                value[ch] += sign*all_rows[s][4]*trace
    return tuple(value)


def check():
    outside = normal_boundary(NORMAL,0,1)
    assert outside == -4548608064
    with TemporaryDirectory() as tmp:
        native = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(native)],check=True)
        hexagon_coefficients = restrict(native)
    one_loop = hexagon_normal(NORMAL,0,1)
    two_loop = theta_normal(NORMAL)
    assert one_loop == Fraction(653198,15) and two_loop == Fraction(8129517449,170100000)
    for plane, floor in zip(PLANES,(18550,36651),strict=True):
        a,b,d = restricted_symbols(plane)
        assert 22*a-5*b+18*d == 0
        for slopes in (((0,0),)*4,((2,1),(-1,0),(3,2),(0,-4))):
            balancing = tuple(-sum(c) for c in zip(*plane,*slopes))
            assert literal_boundary(plane,NORMAL,(*slopes,balancing)) == (24*a*outside,24*b*outside)
        k,j = Fraction(outside*b,576),Fraction(outside*(8*a-b),9)
        needed = tuple(k*x+j*y for x,y in total_channels(BIPARTITE))
        h0,h3 = (v*one_loop for v in planar_factors(plane,hexagon_coefficients))
        theta = two_loop*planar_factor(plane)
        bare = (needed[0]-h0-theta,needed[1]-h3,needed[2])
        print('eleven-field boundary times (K4,J4)', k,j)
        print('theta',theta,'hexagon',h0,h3)
        low,high = interval(bare)
        print('required bare (1,zeta(3),zeta(5))', *bare,'interval',low,high)
        assert floor < low <= high < floor+1


if __name__ == '__main__':
    check()
