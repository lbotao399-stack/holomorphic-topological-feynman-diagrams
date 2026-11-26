"""One empty star and one parallel-edge contraction in the mixed genus-four prism."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations
from math import factorial, prod
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import json

from flint import fmpq
from fourfold_quartic_prism import project_slices
from mixed_bare_primitive_pole import METRIC, PRIME, PrimitiveChamber, image, numerator, residue
from mixed_prism_empty_stars import prism_data, square_roots
from quartic_bipartite_wheel import dual_powers
from quartic_seven_plane_prism import integer_parallel
from quartic_seven_plane_prism_check import kirchhoff
from quartic_wheel_exact import fields, field_table, reconstruct
from quartic_prism_physical_ibp import compile_program, zeta_interval

MARKS = (0,0,0,6,6,6)
STAR = (1,3,5)
HERE = Path(__file__).resolve().parent


@cache
def reduction(marks,star):
    probe = PrimitiveChamber(3,marks)
    edges = list(probe.core)
    center, = set.intersection(*(set(edges[i]) for i in star))
    neighbors = [next(v for v in edges[i] if v!=center) for i in star]
    for j,i in enumerate(star):
        edges[i] = tuple(sorted(v for k,v in enumerate(neighbors) if k!=j))
    pair = [i for i,e in enumerate(edges) if edges.count(e)>1]
    assert len(pair)==2
    wheel = list(dict.fromkeys(edges))
    representatives = [edges.index(e) for e in wheel]
    # Positive placeholders retain the original path index across planar duality.
    spokes,ring = dual_powers(wheel,[i+1 for i in representatives])
    spokes = tuple(i-1 for i in spokes)
    ring = tuple(i-1 for i in ring)
    # ring[i] joins dual spokes i-1 and i.
    packings = []
    for start in range(4):
        for step in (1,-1):
            s = [spokes[(start+step*i)%4] for i in range(4)]
            r = [ring[(start+step*i+(step<0))%4] for i in range(4)]
            packings.append((r[1],r[3],r[0],s[0],s[3],r[2],s[1],s[2]))
    return tuple(edges),pair,representatives,packings


def chamber_rows(marks=MARKS,core_index=3):
    probe = PrimitiveChamber(core_index,marks)
    stabilizer = [g for g in probe.actions if tuple(sorted(g[i][0] for i in marks))==marks]
    layouts = len(probe.actions)//len(stabilizer)
    colors = set(permutations((0,0,1,1,2,2)))
    while colors:
        color = min(colors)
        config = tuple(tuple(color[i] for i,k in enumerate(marks) if k==e) for e in range(9))
        orbit = {sum(image(config,g),()) for g in stabilizer}
        colors -= orbit
        yield color,8*layouts*len(orbit)


def transformed(color,check,marks=MARKS,star=STAR):
    data,soft = prism_data(marks)
    probe = PrimitiveChamber(3,marks,color)
    groups = []
    for d,rays,scale in square_roots(probe,soft):
        for row in groups:
            if d==row[0]:
                row[1] += scale*soft[rays]
                break
        else:
            groups.append([d,scale*soft[rays]])
    if not groups:
        return
    t = data[1].gens()
    voltage = sum(k*prod(t[i] for i in term) for term,k in probe.voltage_terms)
    metrics = []
    for m in ((METRIC,(11,5,2,7,9,6,8,3,4)) if check else ()):
        new = list(map(F,m))
        delta = sum(new[i]*new[j] for i,j in combinations(star,2))
        for i in star:
            new[i] = delta/new[i]
        jac = prod(new[i] for i in star)/sum(new[i] for i in star)**3
        tree = kirchhoff(reduction(marks,star)[0],new)
        original = probe.numerator(m,check=True)
        expected = original*pow(probe.kirchhoff(m),-13,PRIME)*residue(jac*tree**13)%PRIME
        metrics.append([tuple(fmpq(str(t)) for t in new),expected,F()])
    for poly in project_slices(data,groups,voltage,1,(star,)):
        poly *= fmpq(str(probe.weight/2))
        for item in metrics:
            item[2] += F(str(poly(*item[0])))
        yield poly
    assert all(residue(literal)==expected for _,expected,literal in metrics)


def wheel_orbit(p,packings):
    return min(tuple(p[i] for i in order) for order in packings)


def wheel_bank(polys,marks=MARKS,star=STAR):
    _,(keep,drop),reps,packings = reduction(marks,star)
    weights = defaultdict(int)
    count = 0
    for poly in polys:
        for ex,c in poly.terms():
            p = [int(j)+1 for j in ex]
            beta = integer_parallel(p[keep],p[drop])
            p[keep] += p[drop]-13
            assert sum(p[i] for i in reps)==52 and all(0<p[i]<13 for i in reps)
            key = wheel_orbit(p,packings)
            weights[key,int(c.denominator)] += int(c.numerator)*beta
            count += 1
    bank = defaultdict(F)
    for (key,den),numer in weights.items():
        if numer:
            scale = F(prod(factorial(i-1) for i in key),factorial(12)**2)
            bank[key] += scale*F(numer,den)
    return {sum(i<<(4*j) for j,i in enumerate(reversed(k))):v for k,v in bank.items() if v},count


def project_word(word,directory,check,marks=MARKS,star=STAR):
    bank,count = wheel_bank(transformed(tuple(map(int,word)),check,marks,star),marks,star)
    with gzip.open(directory/(word+'.gz'),'wt') as f:
        for key,value in sorted(bank.items()):
            print(key,value,file=f)
    print(word,count,len(bank),flush=True)



def census():
    from mixed_bare_symbolic import classified_orbits
    core = PrimitiveChamber(3).core
    by_layout = Counter()
    for marks,color,n,weight in classified_orbits(3)[0]:
        empty = [v for v in range(6) if all(v not in core[e] for e in marks)]
        if empty and all((i,j) in core for i,j in combinations(empty,2)):
            by_layout[marks] += 1
    assert len(by_layout)==89 and sum(by_layout.values())==5573 and by_layout[MARKS]==17
    return len(by_layout),sum(by_layout.values())


def read_bank(path):
    with gzip.open(path,'rt') as f:
        return {int(key):fmpq(value) for key,value in (row.split() for row in f)}


def capture(bank,values):
    return tuple(sum((weight*values[key][j] for key,weight in bank.items()),fmpq(0)) for j in range(3))


def wheel_values(keys,directory):
    manifest = directory/'keys'
    manifest.write_text(''.join(f'{k}\n' for k in sorted(keys)))
    wheel = directory/'wheel'
    compile_program('quartic_wheel_bank.cpp',wheel)
    primes = fields(13)
    tables = [field_table(wheel,manifest,directory,None,13,p) for p in primes]
    output = directory/'periods'
    reconstruct(13,primes,tables,manifest,output)
    with output.open() as f:
        return {int(k):tuple(map(fmpq,row)) for k,*row in (line.split() for line in f)}


def rebuild(record,directory):
    for word in record['chambers']:
        project_word(word,directory,True)
    marks = (3,3,3,7,7,7)
    control_word = '002112'
    star = (0,1,2)
    bank,_ = wheel_bank(transformed(tuple(map(int,control_word)),True,marks,star),marks,star)
    control = {k:fmpq(str(v)) for k,v in bank.items()}
    keys = set(control)
    combined = defaultdict(fmpq)
    for word in record['chambers']:
        for k,v in read_bank(directory/(word+'.gz')).items():
            combined[k] += v
            keys.add(k)
    values = wheel_values(keys,directory)
    assert sum(bool(w) for w in combined.values())==record['wheel_orbits']
    for word,row in record['chambers'].items():
        assert capture(read_bank(directory/(word+'.gz')),values)==tuple(map(fmpq,row))
    assert capture(combined,values)==tuple(map(fmpq,record['sum']))
    known = json.loads(HERE.joinpath('mixed_bare_prism_colors.json').read_text())
    expected, = [(fmpq(r),fmpq(s),fmpq(0)) for w,n,r,s in known['chambers'] if w==control_word]
    assert capture(control,values)==expected
    print('Independent two-star tetrahedron control:',*expected)


def validate(record):
    probe = PrimitiveChamber(3,MARKS)
    assert tuple(record['hard_paths'])==MARKS and tuple(record['star'])==STAR
    assert {v for v in range(6) if all(v not in probe.core[h] for h in MARKS)}=={2,5}
    assert (2,5) in probe.core
    layouts = {tuple(sorted(g[h][0] for h in MARKS)) for g in probe.actions}
    rows = list(chamber_rows())
    assert len(layouts)==3 and len(rows)==27 and sum(n for _,n in rows)==3*720
    signed = set()
    for color,n in rows:
        p = PrimitiveChamber(3,MARKS,color)
        assert p.weight*len(probe.actions)==n
        if any(numerator(p.graph,p.sheared,p.flow,[(3,i),(3,j)]) for i,j in combinations(range(4),2)):
            signed.add((''.join(map(str,color)),n))
    assert len(signed)==17 and sum(n for _,n in signed)==1296
    assert {w for w,n in signed}==record['chambers'].keys()
    total = tuple(sum(F(row[j]) for row in record['chambers'].values()) for j in range(3))
    assert total==tuple(map(F,record['sum']))
    intervals = [zeta_interval(n,terms=80) for n in (3,5)]
    def bounds(row):
        lo=hi=F(row[0])
        for c,interval in zip(row[1:],intervals,strict=True):
            a,b=sorted(F(c)*x for x in interval)
            lo+=a
            hi+=b
        return lo,hi
    signs = Counter()
    for row in record['chambers'].values():
        lo,hi = bounds(row)
        assert lo>0 or hi<0
        signs[lo>0]+=1
    assert signs=={True:9,False:8}
    lo,hi = bounds(total)
    assert 58244<10**20*lo<=10**20*hi<58245
    assert 10920<10**20*F(3,16)*lo<=10**20*F(3,16)*hi<10921
    print('17 nonzero chambers, 1296 labeled orders, coefficients:',*total)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true',help='Recompute all color periods and the independent control')
    parser.add_argument('--census',action='store_true',help='Recount the complete one-star prism sector')
    options=parser.parse_args()
    record=json.loads(Path(__file__).with_suffix('.json').read_text())
    validate(record)
    if options.census:
        print('one-star prism placements and nonzero color chambers:',*census())
    if options.rebuild:
        with TemporaryDirectory() as scratch:
            rebuild(record,Path(scratch))
