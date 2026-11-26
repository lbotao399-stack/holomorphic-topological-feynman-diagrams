"""Exact order at t=0 of the real fivefold worldsheet primary on conserved polynomial arcs."""
from functools import cache
from itertools import combinations
from math import inf
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from flint import fmpq, fmpq_poly

from fourfold_contact_order import valuation
from elementary_scalar import ROOT
from three_loop_soft_jet import DOUBLE_SOFT, SUBSETS


def polynomial(entry):
    if isinstance(entry, fmpq_poly):
        return entry
    coefficients = entry if isinstance(entry, (tuple,list)) else [entry]
    return fmpq_poly([fmpq(str(c)) for c in coefficients])


class MomentumMinors:
    def __init__(self, points):
        if len(points) < 7 or len(points)%2 != 1 or any(len(p) != 5 for p in points):
            raise ValueError("Expected 2g+3 five-component momenta, g at least two")
        self.rows = [tuple(map(polynomial, p)) for p in points]
        if any(map(sum, zip(*self.rows))):
            raise ValueError("Momentum must be conserved identically")
        # Per-arc cache, released together with this computation.
        rows = self.rows
        @cache
        def det(slots, columns):
            if not slots:
                return fmpq_poly([1])
            return sum(((-1)**i * rows[slots[0]][c] *
                        det(slots[1:],columns[:i]+columns[i+1:])
                        for i,c in enumerate(columns) if rows[slots[0]][c]), fmpq_poly())
        self.det = det

    def cramer(self, basis, i):
        coordinates = []
        for j in range(5):
            indices = (*basis[:j],i,*basis[j+1:])
            inversions = sum(x>y for k,x in enumerate(indices) for y in indices[k+1:])
            coordinates.append((-1)**inversions*self.det(tuple(sorted(indices)),tuple(range(5)))
                               if len(set(indices)) == 5 else fmpq_poly())
        return coordinates

    def best_basis(self):
        basis = ()
        for i in range(len(self.rows)):
            extended = (*basis,i)
            if any(self.det(extended,c) for c in combinations(range(5),len(extended))):
                basis = extended
            if len(basis)==5:
                break
        if len(basis)<5:
            return None
        while True:
            order = valuation(self.det(basis,tuple(range(5))))
            coordinates = [self.cramer(basis,i) for i in range(len(self.rows))]
            best,i,j = min((valuation(x),i,j) for i,p in enumerate(coordinates) for j,x in enumerate(p))
            if best == order:
                return order, coordinates
            # Each exchange lowers an integer. At termination every row is
            # integral over the valuation ring in this basis, proving global minimality.
            basis = tuple(sorted((*basis[:j],i,*basis[j+1:])))


def contact_order(points):
    minors = MomentumMinors(points)
    basis = minors.best_basis()
    if basis is None:
        return inf
    minimum, coordinates = basis
    contents = [min(map(valuation,p)) for p in coordinates]
    if inf in contents:
        return inf
    genus = (len(points)-3)//2
    half_order = genus*minimum+sum(v-minimum for v in contents)
    rows = [[x.right_shift(v) for x in p] for p,v in zip(coordinates,contents)]
    # Order of projective agreement, measured by the two-row minors of unit rows.
    distances = {(i,j): min(valuation(p[a]*q[b]-p[b]*q[a]) for a,b in combinations(range(5),2))
                 for (i,p),(j,q) in combinations(enumerate(rows),2)}
    def cluster_size(level):
        return max(1+sum(distances[min(i,j),max(i,j)] >= level for j in range(len(rows)) if j!=i)
                   for i in range(len(rows)))
    if cluster_size(inf) > genus:
        return inf
    previous = 0
    for depth in sorted(set(distances.values())-{0,inf}):
        half_order += (depth-previous)*max(0,cluster_size(depth)-genus)
        previous = depth
    return 2*half_order


def covering_order(points):
    """Separate, exponential small-frame test via minimum-weight covering bracket products."""
    m = MomentumMinors(points)
    bases = [(sum(1<<i for i in b),valuation(m.det(b,tuple(range(5)))))
             for b in combinations(range(len(points)),5)]
    best = {0:0}
    for _ in range((len(points)-3)//2):
        following = {}
        for mask,value in bases:
            if value == inf:
                continue
            for covered,cost in best.items():
                union = covered|mask
                following[union] = min(following.get(union,inf),cost+value)
        best = following
    return 2*best.get((1<<len(points))-1,inf)


def three_loop_order(points, exe):
    """Evaluate every polynomial in the independent integrated trace, with exact coefficients."""
    m = MomentumMinors(points)
    factors = [fmpq_poly([1])]
    by_coefficients = {}
    ids = []
    for slots in SUBSETS:
        p = m.det(tuple(slots),tuple(range(5)))
        if not p:
            ids.append(0)
            continue
        sign = 1 if p.leading_coefficient()>0 else -1
        key = tuple(map(str,sign*p))
        if key not in by_coefficients:
            by_coefficients[key] = len(factors)
            factors.append(sign*p)
        ids.append(sign*by_coefficients[key])
    records = [(len(DOUBLE_SOFT),1,-1),*DOUBLE_SOFT,ids]
    out = run([exe],input='\n'.join(' '.join(map(str,row)) for row in records),
              text=True,capture_output=True,check=True).stdout
    orbit = []
    p = fmpq_poly()
    for line in out.splitlines():
        if line=='end':
            orbit.append(p)
            p = fmpq_poly()
        else:
            coefficient,*indices = map(int,line.split())
            term = polynomial(coefficient)
            for i in indices:
                term *= factors[i]
            p += term
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words'])) == DOUBLE_SOFT and len(orbit)==len(DOUBLE_SOFT)
    return min(valuation(sum(q*fmpq(x) for q,x in zip(orbit,col)))
               for col in zip(*data['rational_and_zeta']))


def clustered_fixture(genus):
    t = fmpq_poly([0,1])
    line_count = 2*genus-2
    points = [[1,0,0,0,0] for _ in range(line_count)]
    points[-1][1] = t
    points[-2][2] = t**3
    for j in range(max(0,genus-4)):
        points[-3-j][1+j%4] = t**(4+j)
    points[-2] = [t**2*x for x in points[-2]]
    points.extend([int(i==j) for j in range(5)] for i in range(1,5))
    points.append([-sum(col) for col in zip(*points)])
    return points


def check():
    t = fmpq_poly([0,1])
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for genus,expected in ((3,6),(4,12),(5,20)):
            points = clustered_fixture(genus)
            assert contact_order(points)==expected==covering_order(points)
            if genus==3:
                assert three_loop_order(points,exe)==expected
                scaled = [[a*f for a,f in zip(p,(1,t,1,1,t**2))] for p in points]
                assert contact_order(scaled)==expected+6*genus==three_loop_order(scaled,exe)
            # General-position constant changes of frame, slot order and longitudinal shear.
            sheared = [[p[0],p[1]+(1+t)*p[0],p[2]+p[1],p[3],p[4]] for p in points[::-1]]
            assert contact_order(sheared)==expected
            print("genus",genus,"nested fivefold contact",expected)
    identical = clustered_fixture(4)
    identical.pop()
    for p in identical[:6]:
        p[1:] = [0,0,0,0]
    identical.append([-sum(col) for col in zip(*identical)])
    assert contact_order(identical)==inf
    print("Independent integrated three-loop channels: orders 6 and 24.")


if __name__=='__main__':
    check()
