"""The mixed bipartite occupancy with two triples on disjoint paths."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import combinations
from math import factorial,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip

from flint import fmpq
from mixed_bipartite_local_ibp import compile_program,reduced,zeta_interval
from mixed_bipartite_one_star import ARMS,HERE,METRIC,PRIME,PRISM,TABLE,PrimitiveChamber,kirchhoff,residue
from mixed_prism_one_star import chamber_rows
from mixed_prism_rung_sum import ColorMoments

PERIOD=(F(102163573657672606876837,40763513236013568000000000),
                   F(43027855856219197,6419046553920000000),
                   -F(6178645771,606485880000))
MARKS=(0,0,0,5,5,5)


def all_core_periods(sector,unit=False):
    from collections import Counter
    import json
    from mixed_bare_symbolic import classified_orbits,cubic_cores
    from mixed_bipartite_three_paths import bounds

    folder=Path(__file__).resolve().parent
    def read(name):
        return json.loads((folder/name).read_text())

    colors=Counter()
    orders=Counter()
    for core,edges in enumerate(cubic_cores()):
        for marks,color,n,weight in classified_orbits(core)[0]:
            paths=set(marks)
            if sector=='matching':
                selected=len(paths)==3 and len({vertex for e in paths for vertex in edges[e]})==6
            elif sector=='three':
                selected=len(paths)==3
            elif sector=='paired':
                selected=sorted(Counter(marks).values())==[2,2,2]
            else:
                assert sector=='two'
                selected=len(paths)<=2
                if selected:
                    assert len(paths)==2 and len({vertex for e in paths for vertex in edges[e]})==4
            if selected:
                colors[core,marks]+=1
                orders[core,marks]+=n

    if unit:
        from flint import fmpz_mat
        from mixed_unordered_hard import labeled_cubes,placement
        lengths=(1,)*9
        value=F(0)
        for core,marks in colors:
            probe=PrimitiveChamber(core)
            cubes,stabilizer=labeled_cubes(probe,marks,lengths,True)
            density=placement(probe,cubes,stabilizer,lengths,True)
            cycle=fmpz_mat(probe.cycles)
            determinant=int((cycle*cycle.transpose()).det())
            assert determinant==(45,36,56,75,81)[core]
            value+=density/determinant**4
        print('Unit density on all selected placements',value)
        return value

    shapes={}
    for row in read('mixed_bare_double_core.json')['shapes']:
        key=0,tuple(row['marks'])
        if key in colors:
            total=sum(F(value) for color,n,value in row['chambers'])
            assert total==F(row['subtotal']) and len(row['chambers'])==colors[key]
            shapes[key]=(total,F(0),F(0))
    single_counts=Counter()
    for line in (folder/'mixed_bare_single_core.jsonl').read_text().splitlines():
        marks,color,weight,rational,zeta3=json.loads(line)
        key=2,tuple(marks)
        if key in colors:
            single_counts[key]+=1
            a,b,_=shapes.get(key,(F(0),F(0),F(0)))
            shapes[key]=(a+F(rational),b+F(zeta3),F(0))
    for line in (folder/'mixed_bare_single_sum.jsonl').read_text().splitlines():
        row=json.loads(line)
        key=2,tuple(row['paths'])
        if key not in colors:continue
        assert row['nonzero']==colors[key] and row['labeled']==orders[key]
        total=tuple(map(F,row['period']))+(F(0),)
        if single_counts[key]==colors[key]:assert shapes[key]==total
        shapes[key]=total
        single_counts[key]=colors[key]
    assert all(n==colors[k] for k,n in single_counts.items())
    census=(len(colors),sum(colors.values()),sum(orders.values()))
    assert census=={'matching':(11,328,42720),'two':(5,98,14256),
                    'paired':(15,546,84320),'three':(63,2809,458000)}[sector]
    prism=PrimitiveChamber(3)
    def add(marks,coefficients):
        canonical=min(tuple(sorted(g[e][0] for e in marks)) for g in prism.actions)
        key=3,canonical
        if key not in colors:return
        columns=tuple(map(F,coefficients))
        total=columns+(F(0),)*(3-len(columns))
        assert shapes.setdefault(key,total)==total

    if sector=='two':
        for name in ('mixed_bare_prism_colors.json','mixed_prism_one_star.json'):
            row=read(name)
            add(row.get('paths',row.get('hard_paths')),row['sum'])
        shapes[4,MARKS]=PERIOD
        assert all(sorted(Counter(marks).values())==[3,3] for _,marks in colors)
    else:
        for word,row in read('mixed_prism_direct_periods.json').items():
            key=3,tuple(map(int,word))
            if key in colors:
                assert len(row['colors'])==colors[key] and row['nonzero_hard_orders']==orders[key]
                add(key[1],row['period'])
        for occup,row in read('mixed_prism_empty_stars.json').items():
            marks=tuple(e for c,e in zip(occup,(3,5,7)) for _ in range(int(c)))
            add(marks,row['sum'])
        for occup in ('312','222'):
            row=read('mixed_prism_empty_stars_'+occup+'.json')
            add(tuple(e for c,e in zip(occup,(3,5,7)) for _ in range(int(c))),row['partial_sum'])
        for name,prefix,field in (('opposite_face','00011','sum'),('rung_sum','00222','columns')):
            for last,row in read('mixed_prism_'+name+'.json')['placements'].items():
                add(tuple(map(int,prefix+last)),row[field])
        for word,row in read('mixed_bipartite_three_paths.json').items():
            key=4,tuple(map(int,word))
            if key in colors:shapes[key]=tuple(map(F,row['period']))
        for word,row in read('mixed_bipartite_matching_periods.json').items():
            key=4,tuple(map(int,word))
            if key in colors:
                assert len(row['colors'])==colors[key] and row['total']['nonzero_hard_orders']==orders[key]
                shapes[key]=tuple(map(F,row['total']['period']))
    missing=set(colors)-set(shapes)
    if missing:
        raise ValueError(f'Missing integrated placement periods: {sorted(missing)}')
    assert set(shapes)==set(colors)
    total=tuple(sum(col) for col in zip(*shapes.values(),strict=True))
    print('All-core totals, placements, nonzero colors, hard orders',len(shapes),sum(colors.values()),sum(orders.values()))
    print('Rational, zeta(3), zeta(5)',*total)
    print('Signed bounds in units of 10^-20',bounds(total),'Ward weight five',F(3,16)*total[2])
    return total


def build(program,tmp,control=False):
    probe=PrimitiveChamber(4,MARKS)
    batch=ColorMoments(MARKS,ARMS,4)
    edges=list(probe.core)
    edges[4],edges[6],edges[8]=(2,3),(1,3),(1,2)
    checks=[]
    for metric in (METRIC,(11,5,2,7,9,6,8,3,4)):
        right=list(map(F,metric))
        delta=sum(right[i]*right[j] for i,j in combinations(ARMS,2))
        for i in ARMS:right[i]=delta/right[i]
        jac=prod(right[i] for i in ARMS)/sum(right[i] for i in ARMS)**3
        checks.append((metric,pow(probe.kirchhoff(metric),-13,PRIME),
                       residue(jac*kirchhoff(edges,right)**13),
                       tuple(fmpq(str(x)) for x in right),[0,F()]))
    eligible=orders=nonzero_orders=0
    for color,count in chamber_rows(MARKS,4):
        if control and color!=probe.color:continue
        eligible+=1
        orders+=count
        before=len(batch.nonzero)
        batch.add(color,F(count,len(probe.actions)),checks)
        if len(batch.nonzero)>before:nonzero_orders+=count
    assert (eligible,len(batch.nonzero),orders,nonzero_orders)==((1,1,576,576) if control else (27,17,12960,7776))
    print('Color orbits, nonzero, hard orders, nonzero hard orders:',eligible,len(batch.nonzero),
          orders,nonzero_orders,flush=True)
    raw=tmp/'weights'
    count=0
    with raw.open('w') as stream:
        for poly in batch.projected():
            for _,_,_,right,state in checks:state[1]+=F(str(poly(*right)))
            for ex,c in poly.terms():
                powers=tuple(12-int(ex[i]) for i in PRISM)
                assert sum(powers)==52 and 0<min(powers)<=max(powers)<13
                factor=fmpq(prod(factorial(a-1) for a in powers),factorial(12))
                print(*powers,c*factor,0,file=stream)
                count+=1
    for _,_,_,_,state in checks:assert state[0]==residue(state[1])
    with raw.open('rb') as source:
        result=run([program,'13','--full-components'],stdin=source,capture_output=True,check=True)
    if control:assert result.stdout==gzip.decompress(TABLE.read_bytes())
    rows=[line.split() for line in result.stdout.decode().splitlines()]
    assert (count,len(rows))==((105218,37391) if control else (105222,37391))
    print('Projected monomials, nonzero reflected components:',count,len(rows),flush=True)
    return rows


def integral(tmp,complete,ibp,wheel):
    rows=build(complete,tmp)
    seeds=[tuple(map(int,r[:9])) for r in rows]
    values=reduced(seeds,13,ibp,wheel,tmp)
    total=tuple(sum(F(r[9])*values[i][c] for i,r in enumerate(rows)) for c in range(3))
    lo=hi=total[0]
    for v,n in zip(total[1:],(3,5),strict=True):
        a,b=sorted(v*z for z in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    assert total==PERIOD
    assert 663316<10**20*lo<=10**20*hi<663317
    print('Two disjoint triples, three columns:',*total,flush=True)
    print('Strict bounds, units of 10^-20:',(10**20*lo)//1,-((-10**20*hi)//1),flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    sectors=parser.add_mutually_exclusive_group()
    sectors.add_argument('--all-cores',dest='sector',action='store_const',const='two',help='sum all integrated two-path sectors')
    sectors.add_argument('--perfect-matchings',dest='sector',action='store_const',const='matching',help='sum all cores on three disjoint hard paths')
    sectors.add_argument('--paired',dest='sector',action='store_const',const='paired',help='sum all three-path placements of two hard marks on each path')
    sectors.add_argument('--three-paths',dest='sector',action='store_const',const='three',help='sum all placements with three occupied hard paths')
    parser.add_argument('--unit',action='store_true',help='compute the all-core sector at nine unit lengths')
    args=parser.parse_args()
    if args.sector:
        all_core_periods(args.sector,args.unit)
        raise SystemExit
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        complete,ibp,wheel=tmp/'complete',tmp/'ibp',tmp/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_prism_local_ibp.cpp',ibp)
        compile_program('quartic_wheel_bank.cpp',wheel)
        if args.check:
            from mixed_bare_symbolic import classified_orbits
            supported=[item for item in classified_orbits(4)[0] if len(set(item[0]))<=2]
            assert len(supported)==17 and {item[0] for item in supported}=={MARKS}
            build(complete,tmp,True)
        integral(tmp,complete,ibp,wheel)
