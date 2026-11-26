"""Three and two hard insertions on one triangle, one on the opposite triangle."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import json

from flint import fmpq
from mixed_bare_ranks import witness_counts
from mixed_prism_one_star import (
    PrimitiveChamber, capture, chamber_rows, image, numerator, read_bank,
    transformed, wheel_bank, wheel_values, zeta_interval,
)


def geometry():
    probe = PrimitiveChamber(3)
    triangles = [tuple(i for i,(u,v) in enumerate(probe.core) if u in side and v in side)
                 for side in ((0,1,2),(3,4,5))]
    by_orbit = defaultdict(set)
    for first,second in permutations(triangles):
        for triple,double in permutations(first,2):
            for single in second:
                marks = tuple(sorted((triple,)*3+(double,)*2+(single,)))
                orbit = min(tuple(sorted(g[i][0] for i in marks)) for g in probe.actions)
                by_orbit[orbit].add(marks)
    assert sorted(map(len,by_orbit.values()))==[12,12,12]
    assert sorted(by_orbit)==[(0,0,0,1,1,i) for i in (6,7,8)]
    return by_orbit


def unique_star(marks):
    core = PrimitiveChamber(3,marks).core
    center, = (v for v in range(6) if all(v not in core[i] for i in marks))
    return tuple(i for i,e in enumerate(core) if center in e)


def check_placement(marks,record):
    rows = [(marks,color,n,F(n,12)) for color,n in chamber_rows(marks)]
    assert len(rows)==90 and sum(row[2] for row in rows)==12*720
    missing = []
    assert sum(witness_counts(PrimitiveChamber(3,marks).core,rows,missing))==54
    for _,color,n,weight in missing:
        probe = PrimitiveChamber(3,marks,color)
        assert probe.weight==weight and n==96
        assert all(not numerator(probe.graph,probe.sheared,probe.flow,[(3,i),(3,j)])
                   for i,j in combinations(range(4),2))
    signed = {''.join(map(str,color)):n for m,color,n,w in set(rows)-set(missing)}
    assert len(signed)==54 and sum(signed.values())==5184 and signed.keys()==record['chambers'].keys()
    columns = tuple(sum(F(row[j]) for row in record['chambers'].values()) for j in range(3))
    assert columns==tuple(map(F,record['sum']))
    return columns


def rotate_to_zero(marks,color,star):
    probe = PrimitiveChamber(3,marks,color)
    config = tuple(tuple(color[i] for i in group) for group in probe.groups)
    move = next(g for g in probe.actions if {g[i][0] for i in star}=={0,1,2})
    new_config = image(config,move)
    paths = tuple(i for i,group in enumerate(new_config) for _ in group)
    word = sum(new_config,())
    return paths,word,(0,1,2)


def relabeled_control(marks,star,saved):
    color = (0,0,2,1,1,2)
    paths,word,arms = rotate_to_zero(marks,color,star)
    bank,count = wheel_bank(transformed(word,True,paths,arms),paths,arms)
    assert bank.keys()==saved.keys() and all(fmpq(str(v))==saved[k] for k,v in bank.items())
    return count,len(bank)


def integrate(marks,record,directory):
    star = unique_star(marks)
    for word in record['chambers']:
        paths,color,arms = rotate_to_zero(marks,tuple(map(int,word)),star)
        bank,_ = wheel_bank(transformed(color,True,paths,arms),paths,arms)
        with gzip.open(directory/(word+'.gz'),'wt') as f:
            for k,v in sorted(bank.items()):
                print(k,v,file=f)
    combined = defaultdict(fmpq)
    for word in record['chambers']:
        for k,v in read_bank(directory/(word+'.gz')).items():
            combined[k] += v
    assert sum(bool(v) for v in combined.values())==record['wheel_orbits']
    values = wheel_values(set(combined),directory)
    for word,row in record['chambers'].items():
        assert capture(read_bank(directory/(word+'.gz')),values)==tuple(map(fmpq,row))
    assert capture(combined,values)==tuple(map(fmpq,record['sum']))
    color = (0,0,2,1,1,2)
    direct,_ = wheel_bank(transformed(color,True,marks,star),marks,star)
    saved = read_bank(directory/'002112.gz')
    assert direct.keys()==saved.keys() and all(fmpq(str(v))==saved[k] for k,v in direct.items())
    print('Direct and relabeled chamber agree on',len(saved),'wheel powers.')


def bounds(columns):
    lo = hi = columns[0]
    for weight,n in zip(columns[1:],(3,5),strict=True):
        a,b = sorted(weight*z for z in zeta_interval(n,terms=80))
        lo += a
        hi += b
    return lo,hi


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true',help='Reintegrate every recorded color')
    parser.add_argument('--relabel',action='store_true',help='Compare every wheel weight of one color per placement under a prism symmetry')
    parser.add_argument('--edge',type=int,choices=(6,7,8),help='Select the singleton edge')
    opts=parser.parse_args()
    records=json.loads(Path(__file__).with_suffix('.json').read_text())
    sums=[]
    for marks in sorted(geometry()):
        if opts.edge is not None and marks[-1]!=opts.edge:
            continue
        record=records['placements'][str(marks[-1])]
        cols=check_placement(marks,record)
        assert bounds(cols)[0]>F(record['bound'],10**20)
        assert bounds(cols)[1]<F(record['bound']+1,10**20)
        sums.append(cols)
        print(marks,record['wheel_orbits'],*cols,flush=True)
        if opts.relabel:
            star=unique_star(marks)
            direct,_=wheel_bank(transformed((0,0,2,1,1,2),True,marks,star),marks,star)
            saved={k:fmpq(str(v)) for k,v in direct.items()}
            print('Relabeled control:',*relabeled_control(marks,star,saved),flush=True)
        if opts.rebuild:
            with TemporaryDirectory() as work:
                integrate(marks,record,Path(work))
    if opts.edge is None:
        total=tuple(map(sum,zip(*sums)))
        assert total==tuple(map(F,records['sum']))
        assert bounds(total)[0]>F(records['bound'],10**20)
        assert bounds(total)[1]<F(records['bound']+1,10**20)
        print('36 placements, 162 color orbits, 15552 labeled orders:',*total)
