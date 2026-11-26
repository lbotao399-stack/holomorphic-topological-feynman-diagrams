"""Single numerator chord on a four-wheel in the mixed Ward chamber."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, prod
from pathlib import Path
import gzip
import json
import sys
from subprocess import run
from tempfile import TemporaryDirectory

from flint import fmpq as Q

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from quartic_bipartite_wheel import connection, geometry, jets_on_pole, rational_sum
from quartic_prism_completion import completed
from quartic_prism_physical_ibp import compile_program


def scalar_converges(spokes, ring, vertex, q, degree=1):
    graph = [[0]*6 for _ in range(6)]
    for i in range(4):
        graph[4][i] = graph[i][4] = spokes[i] - degree*int(i == vertex)
        graph[i][(i-1)%4] = graph[(i-1)%4][i] = ring[i]
    for i in range(5):
        graph[i][5] = graph[5][i] = 2*q - sum(graph[i])
    assert sum(graph[5]) == 2*q
    return (all(w < q for r in graph for w in r) and
            all(graph[i][j]+graph[j][k]+graph[k][i] < 2*q for i,j,k in combinations(range(6),3)))


def charts(graph, q, degree=1):
    for i,j in combinations(range(6),2):
        if graph[i][j] != 0:
            continue
        remaining = set(range(6))-{i,j}
        for k,l in combinations(sorted(remaining),2):
            if graph[k][l] != -degree:
                continue
            rim = (k,i,l,j)
            ring = tuple(graph[rim[h-1]][v] for h,v in enumerate(rim))
            for center in sorted(remaining-{k,l}):
                spokes = tuple(graph[center][v] for v in rim)
                if min(spokes+ring)>0 and all(scalar_converges(spokes,ring,h,q,degree) for h in (0,2)):
                    yield spokes,ring


def period(spokes, ring, q):
    """Normalized cut, with |x_0-x_2|^2 on the angular diagonal."""
    assert sum(spokes+ring)==4*q+1 and all(scalar_converges(spokes,ring,h,q) for h in (0,2))
    lam = q-1
    contour,cuts = geometry()
    total = defaultdict(Q)
    # First two terms come from r_0^2 + r_2^2. The cross term -2 r_0 r_2
    # pairs degree n on one half-rim with n+1 on the other half.
    sectors = (((0,0,0,0),(2,0,0,0),False), ((0,0,0,0),(0,0,2,0),False),
               ((0,1,1,0),(1,0,1,0),True), ((1,0,0,1),(1,0,1,0),True))
    for higher, radial_numerator, cross in sectors:
        angular = Counter({i:-1 for i in range(1,2*lam)})
        angular[lam] += 3
        if cross:
            # Ratio to the scalar angular measure, including the factor -2:
            # -(n+lambda)(n+2 lambda)/(n+lambda+1)^2.
            angular[lam] -= 1
            angular[2*lam] -= 1
            angular[lam+1] += 2
        radial = {}
        for key,(vertices,crossing,internal) in cuts.items():
            shift = sum(radial_numerator[i] for i in vertices)+sum(higher[i] for i in crossing)
            assert shift%2 == 0
            radial[key] = sum(q-spokes[i] for i in vertices)-sum(ring[i] for i in internal)+shift//2
        assert min(radial.values())>0
        for ks in product(*(range(q-x) for x in ring)):
            base = angular.copy()
            weight = Q(-1 if cross else 1)
            for power,k,delta in zip(ring,ks,higher):
                coeff,roots = connection(q,power,k)
                base.update(root+delta for root in roots)
                weight *= Q(str(coeff))
            rad = {cut:h+sum(ks[i] for i in cuts[cut][1]) for cut,h in radial.items()}
            for roots,multiplicity in Counter(tuple(sorted(rad[i] for i in row)) for row in contour).items():
                factors = base.copy()
                factors.update(roots)
                profile = tuple(sorted((a,p) for a,p in factors.items() if p))
                for k,value in rational_sum(profile).items():
                    total[k] += weight*multiplicity*value
        jets_on_pole.cache_clear()
    assert all(not c or k in (0,3,5) for k,c in total.items())
    norm = Q(lam**3, factorial(lam)**3*factorial(2*lam-1))
    return tuple(F(str(total[k]*norm)) for k in (0,3,5))


def controls():
    ref = json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    known = {}
    for ex,cols in zip(ref['exponents'],ref['coefficients'],strict=True):
        p = tuple(a+1 for a in ex)
        expected = tuple(F(str(c))*factorial(5)/prod(map(factorial,ex)) for c in cols)
        for wheel in charts(completed(p,6),6):
            known.setdefault(wheel,expected)
            assert known[wheel]==expected
    for wheel,val in known.items():
        assert period(*wheel,6)==val
    print('q=6 independent scalar, weight-three and weight-five controls:',len(known))


def physical(check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as f:
        rows = [r for line in f if (r:=line.split())[11]=='0']
    with TemporaryDirectory() as scratch:
        exe = Path(scratch)/'complete'
        compile_program('quartic_prism_completion.cpp',exe)
        data = ''.join(' '.join(row[:11])+'\n' for row in rows)
        out = run([exe,'13','--wheel-charts'],input=data,text=True,capture_output=True,check=True)
    old = {tuple(map(int,s.split()[:9])) for s in out.stdout.splitlines()}
    sums = [F(),F(),F()]
    count = 0
    for r in rows:
        seed = tuple(map(int,r[:9]))
        if seed in old:
            continue
        wheels = list(charts(completed(seed,13),13))
        selected = [w for w in wheels if sorted(w[1])==[6,9,11,12]]
        if not selected:
            continue
        val = period(*selected[0],13)
        if check:
            for w in set(selected[1:]):
                assert period(*w,13)==val
        for i,c in enumerate(val):
            sums[i] += F(r[9])*c
        count += 1
    assert count == 10 and tuple(sums) == (
        -F(15469511507765307278569,146309488733719525294080000000),
        -F(5575670452211,19748068342203801600),F(579197,1349265677760))
    print('q=13 chamber components with rim powers 6,9,11,12:',count,*sums)
    return count,tuple(sums)


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Verify both origins and the independent q=6 periods')
    args = parser.parse_args()
    if args.check:
        controls()
    physical(args.check)
