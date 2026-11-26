"""One empty star of the mixed K3,3 chamber on paths (0,0,0,5,5,5)."""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import hashlib

from flint import fmpq, fmpz_mpoly_ctx
from fourfold_quartic_prism import project_squares
from mixed_bare_primitive_pole import METRIC, PRIME, PrimitiveChamber, residue
from quartic_seven_plane_double import numerator
from quartic_seven_plane_prism_check import kirchhoff
from quartic_seven_plane_single import orbit, tetrahedron_bank
from quartic_prism_physical_ibp import compile_program


ARMS = (4,6,8)
# Core path indices on the dual prism in quartic_prism_completion.FACES order.
PRISM = (8,0,1,5,3,2,6,4,7)
HERE = Path(__file__).resolve().parent
TABLE = HERE/'mixed_bipartite_one_star.components.gz'
DIGEST = 'a935b0347332e14a15200bac11f097f83eeb981fa091c3f399683586a6ee4624'


def transformed_chamber():
    probe = PrimitiveChamber(4)
    ctx = fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(15)))
    times = ctx.gens()
    paths = []
    for group in probe.groups:
        start = sum(map(len,paths))
        paths.append(tuple(range(start,start+len(group)+1)))
    lengths = [sum(times[i] for i in path) for path in paths]
    minors = {(i,j):numerator(probe.graph,probe.sheared,probe.flow,[(3,i),(3,j)])
              for i,j in combinations(range(4),2)}
    squares = []
    for a,b in combinations(range(9),2):
        det = sum((probe.cycles[i][a]*probe.cycles[j][b]-probe.cycles[i][b]*probe.cycles[j][a])*poly
                  for (i,j),poly in minors.items())
        if not det:
            continue
        unit = int(det.content())*(1 if det.leading_coefficient()>0 else -1)
        det = det//unit
        for term in squares:
            if term[0]==det:
                term[1] += lengths[a]*lengths[b]*unit**2
                break
        else:
            squares.append([det,lengths[a]*lengths[b]*unit**2])
    assert [(len(p),len(w)) for p,w in squares]==[(183,15)]
    voltage = sum(k*prod(times[i] for i in ex) for ex,k in probe.voltage_terms)
    poly = project_squares((probe.graph,ctx,None,lengths,None,paths,None),
                           squares,voltage,1,(ARMS,))
    assert len(poly)==105218 and probe.weight==8
    edges = list(probe.core)
    edges[4],edges[6],edges[8] = (2,3),(1,3),(1,2)
    for metric in (METRIC,(11,5,2,7,9,6,8,3,4)):
        original = probe.numerator(metric,check=True)
        new = list(map(F,metric))
        delta = sum(new[i]*new[j] for i,j in combinations(ARMS,2))
        for i in ARMS:
            new[i] = delta/new[i]
        jacobian = prod(new[i] for i in ARMS)/sum(new[i] for i in ARMS)**3
        tree = kirchhoff(edges,new)
        value = F(str(poly(*(fmpq(str(x)) for x in new))))
        assert residue(probe.weight/2*value/tree**13)==(
            original*pow(probe.kirchhoff(metric),-13,PRIME)*residue(jacobian)%PRIME)
    return poly,probe.weight/2


def rebuild():
    poly,weight = transformed_chamber()
    with TemporaryDirectory() as scratch:
        tmp = Path(scratch)
        program = tmp/'complete'
        compile_program('quartic_prism_completion.cpp',program)
        with (tmp/'weights').open('w') as sink:
            for ex,c in poly.terms():
                powers = tuple(12-int(ex[i]) for i in PRISM)
                assert sum(powers)==52 and 0<min(powers) and max(powers)<13
                factor = F(prod(factorial(a-1) for a in powers),factorial(12))
                print(*powers,c*fmpq(str(weight*factor)),0,file=sink)
        with (tmp/'weights').open('rb') as source:
            result = run([program,'13','--full-components'],stdin=source,capture_output=True,check=True)
        assert result.stderr.decode().splitlines()==[
            'rows 105218 completed 47678 nonzero 47676 minimum margin 1',
            '4262 37403 37391 0 12 26788']
        return result.stdout


def subtotal(raw):
    assert hashlib.sha256(raw).hexdigest()==DIGEST
    counts = Counter()
    powers = set()
    total = [F(),F()]
    for line in raw.decode().splitlines():
        row = line.split()
        assert len(row) in (12,19) and row[10]=='0'
        counts[len(row),int(row[11])] += 1
        if len(row)==19:
            chart = orbit(tuple(map(int,row[12:18])))
            weight = F(row[9])*F(row[18])
            powers.add(chart)
            for j,v in enumerate(tetrahedron_bank()[chart]):
                total[j] += weight*v
    assert counts=={(12,0):35563,(12,1):1242,(19,1):586} and len(powers)==461
    assert total==[F(13780571003319104322144109,29035118039206639794610176000000000),
                   -F(1842108151731299,4665481145845648128000000)]
    return total


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true',help='Regenerate the chamber and all completed cuts')
    args = parser.parse_args()
    raw = gzip.decompress(TABLE.read_bytes())
    if args.rebuild:
        assert rebuild()==raw
    print('586 scalar-star components, rational and zeta(3) subtotal:',*subtotal(raw))
    print('1242 more components are integrated by mixed_bipartite_tensor.py; 35563 have no reachable bubble.')
