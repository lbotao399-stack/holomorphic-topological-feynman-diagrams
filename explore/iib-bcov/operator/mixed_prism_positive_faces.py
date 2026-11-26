"""Recover positive prism periods from five-term identities at a triangular face."""
from argparse import ArgumentParser
from functools import lru_cache
from heapq import merge
from itertools import groupby
from math import factorial,lcm,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json
import gzip
import sqlite3

from flint import fmpq
from mixed_prism_scalar_functional import FACE_ORDER,HERE,batches,check_positive,compile_program,pack


def control(ibp):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    seeds={tuple(i+1 for i in ex):tuple(fmpq(str(c))*factorial(5)/prod(map(factorial,ex)) for c in values)
           for ex,values in zip(ref['exponents'],ref['coefficients'],strict=True)}
    def chart(p):
        g=[[0]*6 for _ in range(6)]
        for (i,j),a in zip(FACE_ORDER,p):g[i][j]=g[j][i]=a
        for i in range(5):g[i][5]=g[5][i]=12-sum(g[i])
        return g
    def label(i,j):return FACE_ORDER.index(tuple(sorted((i,j))))
    rows=[]
    candidates=set(seeds)
    for p in seeds:
        g=chart(p)
        for v in (0,1):
            for r in (2,3,4):
                terms=[(p,g[v][5]-g[v][r])]
                for j in {2,3,4}-{r}:
                    for decrease,sign in ((label(v,r),-1),(label(j,r),1)):
                        h=list(p);h[label(v,j)]+=1;h[decrease]-=1
                        terms.append((tuple(h),sign*g[v][j]))
                if any(min(h)<1 or max(h)>=6 for h,c in terms if c):continue
                rows.append([(h,c) for h,c in terms if c]);candidates.update(h for h,c in terms if c)
    ordered=sorted(candidates)
    result=run([ibp,'6','--canonical-powers'],input=''.join(str(pack(p))+'\n' for p in ordered),
               text=True,capture_output=True,check=True)
    aliases=dict(zip(ordered,map(int,result.stdout.splitlines()),strict=True))
    values={aliases[p]:c for p,c in seeds.items()}
    count=0
    for row in rows:
        if not all(aliases[p] in values for p,c in row):continue
        for i in range(3):
            if sum(c*values[aliases[p]][i] for p,c in row):
                raise ArithmeticError('Face identity disagrees with a measured propagator')
        count+=1
    if not count:raise ArithmeticError('No measured face identity')
    print('Independently integrated q=6 three-column face identities',count,flush=True)


def apply(plan,database,reference):
    db=sqlite3.connect(database,timeout=60.0)
    shared=sqlite3.connect(reference.resolve(strict=True).as_uri()+'?mode=ro',uri=True)
    sql='SELECT r,w3,w5 FROM prism WHERE key=?'
    @lru_cache(maxsize=65536)
    def value(key):
        row=db.execute(sql,(key,)).fetchone() or shared.execute(sql,(key,)).fetchone()
        if row is None:raise ValueError(f'Unknown angular anchor {key}')
        return tuple(map(fmpq,row))
    count=overlap=0
    try:
        opener=gzip.open if plan.suffix=='.gz' else open
        with opener(plan,'rt') as f:
            for lines in batches(f,512):
                cache={}
                for line in lines:
                    fields=list(map(int,line.split()))
                    key,pivot=fields[:2]
                    cols=[fmpq(0) for _ in range(3)]
                    for child,coefficient in zip(fields[2::2],fields[3::2],strict=True):
                        saved=cache.get(child) or value(child)
                        for i,c in enumerate(saved):cols[i]-=coefficient*c
                    cache[key]=tuple(c/pivot for c in cols)
                    check_positive(key,cache[key])
                old={}
                for group in batches(cache,500):
                    marks=','.join('?' for _ in group)
                    for key,*cols in db.execute(f'SELECT key,r,w3,w5 FROM prism WHERE key IN ({marks})',group):
                        old[key]=tuple(map(fmpq,cols))
                if any(cache[k]!=v for k,v in old.items()):raise ArithmeticError('Independently evaluated angular period differs')
                overlap+=len(old)
                db.executemany('INSERT OR IGNORE INTO prism VALUES (?,?,?,?)',
                               ((k,*map(str,v)) for k,v in cache.items()))
                db.commit()
                count+=len(lines)
                if count%65536==0:print('Exact face reconstructions',count,'independent overlaps',overlap,flush=True)
    finally:
        shared.close();db.close()
    print('Positive angular periods from five-term equations',count,'independent overlaps',overlap,flush=True)


def verify(archive,work):
    native,domain,angular,finite,scaled=(work/s for s in ('faces','domain','angular','finite','scaled'))
    compile_program('operator/mixed_prism_positive_ibp.cpp',native)
    compile_program('operator/mixed_prism_positive_domain.cpp',domain)
    for name,mode in ((angular,[]),(finite,['--finite-aliases'])):
        with name.open('wb') as output:run([domain,'13',*mode],stdout=output,check=True)
    common=[1]*3
    bits=[0]*3
    def ratio(s):
        n,separator,d=s.partition(b'/')
        return int(n),int(d) if separator else 1
    with gzip.open(archive,'rb') as stream,angular.open('rb') as labels:
        for line,label in zip(stream,labels,strict=True):
            k,*cols=line.split()
            if k!=label.strip() or len(cols)!=3:raise ValueError('Changed angular support')
            for j,c in enumerate(cols):
                n,d=ratio(c);bits[j]=max(bits[j],abs(n).bit_length());common[j]=lcm(common[j],d)
    # A row has at most seven terms, total absolute coefficient below 128.
    # Each scaled column entry has absolute value below 2^(bits + denominator bits).
    width=max((n+c.bit_length()+8+63)//64 for n,c in zip(bits,common))
    with gzip.open(archive,'rb') as stream,scaled.open('wb') as out:
        out.write(width.to_bytes(8,'little'))
        for line in stream:
            k,*cols=line.split()
            for c,d in zip(cols,common,strict=True):
                n,v=ratio(c)
                out.write((n*(d//v)).to_bytes(8*width,'little',signed=True))
    print('Exact scaled prism columns, common denominators and limb count:',common,width,flush=True)
    run([native,'13',angular,finite,'/dev/null',work,work/'unused','--all-finite',scaled],check=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('scratch',type=Path)
    p.add_argument('reference',type=Path,nargs='?',help='Optional prior angular prism database for reconstruction')
    p.add_argument('--plan',type=Path)
    p.add_argument('--check',action='store_true')
    p.add_argument('--verify',type=Path,help='Check every convergent five- and seven-term identity over the rationals')
    p.add_argument('--all-finite',action='store_true',help='Use every positive regrounding and the rectangular-face identities')
    args=p.parse_args()
    args.scratch.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(dir=args.scratch) as t:
        work=Path(t)
        if args.check:
            compile_program('quartic_prism_local_ibp.cpp',work/'local')
            control(work/'local')
        if args.verify is not None:
            verify(args.verify,work)
            raise SystemExit
        if args.reference is None and args.check and args.plan is None:
            raise SystemExit
        if args.reference is None:
            args.reference=args.scratch/'prism-13.sqlite'
        plan=args.plan
        if plan is None:
            compile_program('operator/mixed_prism_positive_domain.cpp',work/'domain')
            compile_program('operator/mixed_prism_positive_ibp.cpp',work/'face')
            for name,mode in (('angular',[]),('finite',['--finite-aliases'])):
                with (work/name).open('wb') as data:run([work/'domain','13',*mode],stdout=data,check=True)
            dbpath=args.scratch/'prism-13.sqlite'
            with sqlite3.connect(dbpath) as db,sqlite3.connect(args.reference.resolve(strict=True).as_uri()+'?mode=ro',uri=True) as shared:
                db.execute('CREATE TABLE IF NOT EXISTS prism (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
                query='SELECT key FROM prism ORDER BY key'
                with (work/'anchors').open('w') as f:
                    for (key,),_ in groupby(merge(db.execute(query),shared.execute(query))):print(key,file=f)
            plan=work/'face.plan'
            mode=['--all-finite'] if args.all_finite else []
            run([work/'face','13',work/'angular',work/'finite',work/'anchors',work,plan,*mode],check=True)
        apply(plan,args.scratch/'prism-13.sqlite',args.reference)
