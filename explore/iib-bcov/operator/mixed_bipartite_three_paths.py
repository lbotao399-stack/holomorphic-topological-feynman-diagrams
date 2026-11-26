"""Integrate all color words of each mixed bipartite three-path occupancy."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations
from math import factorial,prod
from multiprocessing import get_context
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json
import os
import sqlite3

from flint import fmpq
from mixed_bipartite_local_ibp import compile_program,reduced,zeta_interval
from mixed_bipartite_one_star import ARMS,METRIC,PRIME,PRISM,PrimitiveChamber,kirchhoff,residue
from mixed_prism_one_star import chamber_rows
from mixed_prism_rung_sum import ColorMoments

COUNTS={'000113':54,'000115':54,'000117':54,'000155':54,'000177':54,'000556':54,
        '001133':38,'001177':38}
HERE=Path(__file__).resolve().parent
RESULT=HERE/'mixed_bipartite_three_paths.json'


@cache
def expected_colors(word):
    if word in COUNTS:return COUNTS[word]
    from mixed_bare_symbolic import classified_orbits
    marks=tuple(map(int,word))
    return sum(m==marks for m,_,_,_ in classified_orbits(4)[0])


def choose_chart(word):
    marks=tuple(map(int,word))
    probe=PrimitiveChamber(4,marks)
    positions={tuple(sorted(g[i][0] for i in marks)) for g in probe.actions}
    # Carry an empty star to the one used for the previously evaluated dual prism.
    marked=min(m for m in positions if set(m).isdisjoint(ARMS))
    return marked,len(positions)


def color_part(word,raw,chunk_size,part,total_parts,stream_roots=False):
    marks,layouts=choose_chart(word)
    probe=PrimitiveChamber(4,marks)
    batch=ColorMoments(marks,ARMS,4)
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
    accumulated={}
    def collect(chunk):
        for degree,poly in zip(chunk.coefficients,chunk.projected(),strict=True):
            if degree in accumulated:
                accumulated[degree]+=poly
            else:
                accumulated[degree]=poly
    nonzero_colors=nonzero_orders=orders=color_count=0
    for ordinal,(color,count) in enumerate(chamber_rows(marks,4)):
        if ordinal%total_parts!=part:continue
        color_count+=1
        orders+=count
        before=len(batch.nonzero)
        batch.add(color,F(count,len(probe.actions)),checks,collect if stream_roots else None)
        if len(batch.nonzero)>before:
            nonzero_orders+=count
            nonzero_colors+=1
        if chunk_size and len(batch.nonzero)>=chunk_size:
            collect(batch)
            batch=ColorMoments(marks,ARMS,4)
    print(word,'grounded color evaluations complete:',color_count,nonzero_colors,flush=True)
    if chunk_size or stream_roots:
        collect(batch)
        slices=accumulated.values()
    else:
        slices=batch.projected()
    count=0
    with raw.open('w') as stream:
        for poly in slices:
            for _,_,_,right,state in checks:state[1]+=F(str(poly(*right)))
            for ex,c in poly.terms():
                powers=tuple(12-int(ex[i]) for i in PRISM)
                assert sum(powers)==52 and 0<min(powers)<=max(powers)<13
                factor=fmpq(prod(factorial(a-1) for a in powers),factorial(12))
                print(*powers,c*factor,0,file=stream)
                count+=1
    for _,_,_,_,state in checks:assert state[0]==residue(state[1])
    return color_count,nonzero_colors,orders,nonzero_orders,count


def combine_weights(files,raw):
    combined={}
    for path in files:
        with path.open() as stream:
            for line in stream:
                *powers,coefficient,zero=line.split()
                assert zero=='0'
                key=bytes(map(int,powers))
                combined[key]=combined.get(key,0)+fmpq(coefficient)
    count=0
    with raw.open('w') as stream:
        for powers,coefficient in combined.items():
            if coefficient:
                print(*powers,coefficient,0,file=stream)
                count+=1
    return count


def finish_projection(word,program,raw,dest,color_count,nonzero_colors,orders,nonzero_orders,count):
    _,layouts=choose_chart(word)
    assert nonzero_colors==expected_colors(word) and orders==720*layouts
    with raw.open('rb') as stream:
        result=run([program,'13','--full-components'],stdin=stream,capture_output=True,check=True)
    stats={'color_orbits':color_count,'nonzero_colors':nonzero_colors,'placements':layouts,
           'hard_orders':orders,'nonzero_hard_orders':nonzero_orders,'monomials':count,
           'components':result.stdout.count(b'\n')}
    target=dest/(word+'.gz')
    pending=target.with_suffix('.partial')
    with gzip.open(pending,'wb') as stream:
        stream.write((json.dumps(stats)+'\n').encode())
        stream.write(result.stdout)
    pending.replace(target)
    print(word,stats,flush=True)


def project(word,program,tmp,dest,chunk_size=0,parts=1,stream_roots=False):
    raw=tmp/'weights'
    if parts==1:
        counts=color_part(word,raw,chunk_size,0,1,stream_roots)
    else:
        files=[tmp/f'weights-{i}' for i in range(parts)]
        with ProcessPoolExecutor(max_workers=parts,mp_context=get_context('spawn')) as pool:
            jobs=[pool.submit(color_part,word,path,chunk_size,i,parts,stream_roots) for i,path in enumerate(files)]
            partial=[job.result() for job in jobs]
        totals=list(map(sum,zip(*partial,strict=True)))
        totals[-1]=combine_weights(files,raw)
        for path in files:path.unlink()
        counts=totals
    finish_projection(word,program,raw,dest,*counts)


def bounds(columns):
    lo=hi=columns[0]
    for v,n in zip(columns[1:],(3,5),strict=True):
        a,b=sorted(v*z for z in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    return [(10**20*lo)//1,-((-10**20*hi)//1)]


def integrate(word,ibp,wheel,tmp,dest):
    with gzip.open(dest/(word+'.gz'),'rb') as stream:
        digest=sha256(stream.read()).hexdigest()
    with gzip.open(dest/(word+'.gz'),'rt') as stream:
        stats=json.loads(next(stream))
        rows=[line.split() for line in stream]
    seeds=[tuple(map(int,r[:9])) for r in rows]
    assert len(seeds)==len(set(seeds))==stats['components']
    with sqlite3.connect(dest/'angular.sqlite') as connection:
        connection.execute('CREATE TABLE IF NOT EXISTS wheels (degree INTEGER, key INTEGER, r TEXT, w3 TEXT, w5 TEXT, PRIMARY KEY(degree,key))')
        values=reduced(seeds,13,ibp,wheel,tmp,connection,prime_workers=3)
    cols=tuple(sum(F(r[9])*values[i][c] for i,r in enumerate(rows)) for c in range(3))
    stats.update(component_sha256=digest,period=list(map(str,cols)),bounds=bounds(cols),ward_zeta5=str(3*cols[2]/16))
    if RESULT.is_file():
        expected=json.loads(RESULT.read_text()).get(word)
        if expected is not None:assert stats==expected
    print(word,stats,flush=True)
    (dest/(word+'.json')).write_text(json.dumps(stats,indent=2)+'\n')


def aggregate(records):
    cols=tuple(sum(F(row['period'][c]) for row in records.values()) for c in range(3))
    return cols,bounds(cols),[F(3,16)*x for x in cols]


def check(records,scratch=None):
    from mixed_bare_symbolic import classified_orbits
    core=PrimitiveChamber(4).core
    colored=Counter()
    ordered=Counter()
    matching=Counter()
    for marks,color,n,w in classified_orbits(4)[0]:
        if len(set(marks))!=3:continue
        word=''.join(map(str,marks))
        if any(all(v not in core[i] for i in marks) for v in range(6)):
            colored[word]+=1
            ordered[word]+=n
        else:matching[word]+=1
    assert colored==COUNTS and matching=={'000558':28,'005588':8}
    for word,row in records.items():
        if scratch is not None:
            with gzip.open(scratch/(word+'.gz'),'rb') as stream:
                assert sha256(stream.read()).hexdigest()==row['component_sha256']
        cols=tuple(map(F,row['period']))
        assert row['nonzero_colors']==colored[word] and bounds(cols)==row['bounds']
        assert row['nonzero_hard_orders']==ordered[word]
        assert row['placements']==choose_chart(word)[1]
        assert F(row['ward_zeta5'])==F(3,16)*cols[2]


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('project','integrate','check'))
    parser.add_argument('words',nargs='*')
    parser.add_argument('--scratch',type=Path)
    parser.add_argument('--chunk',type=int,default=0,help='Project and collect this many nonzero colors at a time')
    parser.add_argument('--parts',type=int,default=1,help='Project disjoint color lists in separate workers')
    parser.add_argument('--stream-roots',action='store_true',help='Contract each squared current minor before processing the next')
    args=parser.parse_args()
    if args.mode=='check':
        records=json.loads(RESULT.read_text())
        check(records,args.scratch)
        total,interval,ward=aggregate(records)
        print('Verified placements and nonzero colors:',len(records),sum(r['nonzero_colors'] for r in records.values()))
        print('Bare columns and strict bounds in 10^-20 units:',*total,interval)
        print('Ward-normalized columns:',*ward)
    else:
        assert args.scratch is not None and args.words and all(len(w)==6 and all(c in '012345678' for c in w) and expected_colors(w) for w in args.words)
        args.scratch.mkdir(exist_ok=True)
        with TemporaryDirectory(dir=args.scratch) as directory:
            tmp=Path(directory)
            if args.mode=='project':
                if 'MIXED_PRISM_MOMENTS_LIB' not in os.environ:
                    native=tmp/'ordered-moments.so'
                    run(['c++','-std=c++17','-O3','-shared','-fPIC',
                         str(HERE/'ordered_path_moments.cpp'),'-o',str(native)],check=True)
                    os.environ['MIXED_PRISM_MOMENTS_LIB']=str(native.resolve())
                complete=tmp/'complete'
                compile_program('quartic_prism_completion.cpp',complete)
                for w in args.words:project(w,complete,tmp,args.scratch,args.chunk,args.parts,args.stream_roots)
            else:
                ibp,wheel=tmp/'ibp',tmp/'wheel'
                compile_program('quartic_prism_local_ibp.cpp',ibp)
                compile_program('quartic_wheel_bank.cpp',wheel)
                for w in args.words:integrate(w,ibp,wheel,tmp,args.scratch)
