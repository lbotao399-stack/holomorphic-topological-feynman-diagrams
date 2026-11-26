"""Contract positive triangular-prism powers against the completed angular periods."""
from argparse import ArgumentParser
from array import array
from collections import defaultdict
from contextlib import ExitStack
from functools import cache
from hashlib import sha256
from heapq import merge
from itertools import groupby, islice
from math import factorial, lcm, prod
from operator import itemgetter
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json
import fcntl
import os
import sqlite3

from flint import fmpq
from mixed_bipartite_local_ibp import CHAMBER, compile_program, reduced
from mixed_bipartite_fixed_spin import evaluate, verify_chart
from mixed_bipartite_star_bank import pack, unpack
from quartic_prism_physical_ibp import PROBES, FULL_PROBES, zeta_interval


HERE = Path(__file__).resolve().parent
FACE_ORDER = ((2,3),(0,2),(0,3),(1,3),(1,2),(0,4),(2,4),(3,4),(1,4))
ARCHIVE = HERE/'mixed_prism_positive_periods.gz'
METADATA = HERE/'mixed_prism_positive_periods.json'


def batches(items, size):
    items = iter(items)
    while batch := list(islice(items, size)):
        yield batch


def archive_reference(scratch):
    metadata=json.loads(METADATA.read_text())
    digest=metadata['uncompressed_sha256']
    dest=scratch/f'angular-{digest}.sqlite'
    fd=os.open(scratch,os.O_RDONLY)
    try:
        fcntl.flock(fd,fcntl.LOCK_EX)
        if dest.is_file():return dest
        with TemporaryDirectory(dir=scratch) as t:
            tmp=Path(t)/'complete.sqlite'
            db=sqlite3.connect(tmp)
            checksum=sha256()
            count=previous=0
            try:
                db.execute('CREATE TABLE prism (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
                db.execute('CREATE TABLE archive (digest TEXT PRIMARY KEY)')
                with gzip.open(ARCHIVE,'rb') as f:
                    for lines in batches(f,8192):
                        rows=[]
                        for line in lines:
                            checksum.update(line)
                            k,*cols=line.split();key=int(k)
                            if key<=previous or len(cols)!=3:raise ValueError('Malformed angular table')
                            rows.append((key,*(c.decode() for c in cols)))
                            previous=key
                        db.executemany('INSERT INTO prism VALUES (?,?,?,?)',rows)
                        count+=len(rows)
                if checksum.hexdigest()!=digest or count!=metadata['angular_classes']:
                    raise ArithmeticError('Positive prism archive checksum or count differs')
                db.execute('INSERT INTO archive VALUES (?)',(digest,))
                db.commit()
            finally:db.close()
            tmp.replace(dest)
    finally:os.close(fd)
    return dest


def read_bank(path, q, raw):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as source:
        for line in source:
            fields = line.split()
            if len(fields) == 2:
                key = int(fields[0])
                if not 0 <= key < (1 << 36):
                    raise ValueError(f'Expected a nine-nibble prism key, got {key}')
                powers = unpack(key, 9)
                weight = fmpq(fields[1])
            elif len(fields) == 10:
                powers = tuple(map(int, fields[:9]))
                if min(powers) < 1 or max(powers) >= q:
                    raise ValueError(f'Power out of range in {powers}')
                key = pack(powers)
                weight = fmpq(fields[9])
            else:
                raise ValueError('Expected a packed key and weight or nine face powers and a weight')
            if raw:
                weight *= fmpq(prod(factorial(p-1) for p in powers), factorial(q-1))
            yield key, weight


def known_periods(db, keys):
    found = {}
    for group in batches(keys, 500):
        marks = ','.join('?' for _ in group)
        for key, *values in db.execute(f'SELECT key,r,w3,w5 FROM prism WHERE key IN ({marks})', group):
            found[key] = tuple(map(fmpq, values))
    return found


@cache
def zeta_bounds():
    scale=10**60
    result=[]
    for n in (3,5):
        low,high=zeta_interval(n,192)
        result.append((int(scale*low//1),-int(-scale*high//1)))
    return scale,result


@cache
def sphere_bound(m):
    q=13
    if not 0<m<2*q:raise ValueError('Inadmissible rung sum')
    scale,_=zeta_bounds()
    return fmpq(scale*factorial(q-1)*(4*q)**(2*q)*(2*m)**m*(2*q-m)**(2*q-m),
                factorial(4*q-1)*3**(q+m))


def check_positive(key,cols):
    scale,enclosures=zeta_bounds()
    powers=unpack(key,9)
    m=powers[0]+powers[6]+powers[7]
    lower=upper=scale*cols[0]
    for coeff,(low,high) in zip(cols[1:],enclosures):
        lower+=coeff*(low if coeff>=0 else high)
        upper+=coeff*(high if coeff>=0 else low)
    if lower*prod(a**a for a in powers)<=sphere_bound(m):
        raise ArithmeticError(f'The scalar at {powers} missed its Euclidean lower bound')
    return lower,upper


def pooled_periods(keys,q,work,db,shared,anchor_keys,group):
    expansions,needed = {},defaultdict(set)
    for batch in batches(keys,group):
        local=reduced([unpack(k,9) for k in batch],q,work/'ibp',work/'wheel',work,
                      prism_keys=anchor_keys,return_leaves=True)
        for i,key in enumerate(batch):
            expansions[key]=local[i]
            for degree,terminal in local[i]:
                needed[degree].add(terminal)
    values={}
    for degree,terminals in sorted(needed.items()):
        if degree==-1:
            values[degree]=known_periods(db,terminals)
            if shared is not None:
                values[degree].update(known_periods(shared,terminals-values[degree].keys()))
            if values[degree].keys()!=terminals:raise ArithmeticError('Missing prism terminal')
            continue
        known={}
        for source in (db,shared):
            if source is None:continue
            for batch in batches(terminals-known.keys(),500):
                marks=','.join('?' for _ in batch)
                sql=f'SELECT key,r,w3,w5 FROM wheels WHERE degree=? AND key IN ({marks})'
                known.update((key,tuple(map(fmpq,cols))) for key,*cols in source.execute(sql,(degree,*batch)))
        for batch in batches(sorted(terminals-known.keys()),2048):
            for key in batch:verify_chart(key,q,degree)
            new=evaluate(q,degree,batch,work/'wheel',work,prime_workers=2)
            rows=[(degree,key,*map(str,c)) for key,c in new.items()]
            db.executemany('INSERT OR IGNORE INTO wheels VALUES (?,?,?,?,?)',rows)
            db.commit()
            known.update((key,tuple(fmpq(x.numerator,x.denominator) for x in cols))
                         for key,cols in new.items())
        values[degree]=known
    periods={}
    for key,expansion in expansions.items():
        result=[fmpq(0) for _ in range(3)]
        for (degree,terminal),weight in expansion.items():
            scale=fmpq(weight.numerator,weight.denominator)
            for i,c in enumerate(values[degree][terminal]):result[i]+=scale*c
        periods[key]=tuple(result)
    return periods


def contract(rows, q, work, scratch, reference=None, group=512, tally=False,bulk=False,consume=None):
    """Rows have angular-normalized weights and nine powers in FACE_ORDER."""
    ibp, wheel = work/'ibp', work/'wheel'
    total = [fmpq(0) for _ in range(3)]
    terms = missing = 0
    with ExitStack() as stack:
        db = sqlite3.connect(scratch/f'prism-{q}.sqlite',timeout=30.0)
        stack.callback(db.close)
        db.execute('CREATE TABLE IF NOT EXISTS prism (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
        db.execute('CREATE TABLE IF NOT EXISTS wheels (degree INTEGER,key INTEGER,r TEXT,w3 TEXT,w5 TEXT,PRIMARY KEY(degree,key))')
        shared = None
        anchor_keys = None
        complete=False
        if reference is not None:
            if q != 13: raise ValueError('The shared angular reference is at half-dimension 13')
            shared = sqlite3.connect(reference.resolve(strict=True).as_uri()+'?mode=ro', uri=True)
            stack.callback(shared.close)
            complete=shared.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='archive'").fetchone()[0]>0
        if q==13 and not complete:
            sources=(db,) if shared is None else (db,shared)
            anchor_keys = work/'known-prisms.bin'
            with anchor_keys.open('wb') as output:
                query='SELECT key FROM prism ORDER BY key'
                merged=merge(*(s.execute(query) for s in sources))
                for keys in batches((p[0] for p,_ in groupby(merged)),32768):
                    array('Q',keys).tofile(output)
            with Path(str(anchor_keys)+'.wheels').open('wb') as output:
                query='SELECT degree,key FROM wheels ORDER BY degree,key'
                merged=merge(*(s.execute(query) for s in sources))
                for keys in batches(((d<<32)|k for (d,k),_ in groupby(merged)),32768):
                    array('Q',keys).tofile(output)
        for batch in batches(rows, 32768 if bulk else 4096):
            query = ''.join(f'{row[0]}\n' for row in batch)
            run_result = run([ibp,str(q),'--canonical-convergent'],input=query,text=True,
                             capture_output=True)
            if run_result.returncode:
                raise ValueError(run_result.stderr.strip())
            aliases = list(map(int, run_result.stdout.splitlines()))
            weights = defaultdict(fmpq)
            for alias, row in zip(aliases, batch, strict=True):
                weights[alias] += row[1]
            weights = {k:v for k,v in weights.items() if v}
            periods = known_periods(shared if complete else db, weights)
            if complete and periods.keys()!=weights.keys():
                raise ArithmeticError('Missing a convergent angular class in the complete archive')
            if shared is not None and not complete:
                periods.update(known_periods(shared, weights.keys()-periods.keys()))
            for keys in batches(sorted(weights.keys()-periods.keys()), len(batch) if bulk else group):
                if bulk:
                    periods.update(pooled_periods(keys,q,work,db,shared,anchor_keys,group))
                else:
                    values = reduced([unpack(k,9) for k in keys],q,ibp,wheel,work,
                                     db if q==13 else None,prime_workers=2,wheel_reference=shared,
                                     prism_keys=anchor_keys)
                    for i,key in enumerate(keys):
                        periods[key] = tuple(fmpq(v.numerator,v.denominator) for v in values[i])
                db.executemany('INSERT OR IGNORE INTO prism VALUES (?,?,?,?)',
                               ((key,*map(str,periods[key])) for key in keys))
                db.commit()
                missing += len(keys)
            for key,weight in weights.items():
                for i,c in enumerate(periods[key]):
                    total[i] += weight*c
            if consume is not None:consume(batch,aliases,periods)
            terms += len(batch)
            if tally: print('Input rows',terms,'new angular periods',missing,flush=True)
    return tuple(total)


def colors(manifest,work,scratch,reference,fast=False):
    from mixed_bipartite_matching_periods import package

    record=json.loads(manifest.read_text())
    word=manifest.name.split('-',1)[0]
    names=[*record['colors'],'direct']
    total=[ [fmpq(0) for _ in range(3)] for _ in names ]
    union=0
    def source(i,name):
        previous=-1
        digest=sha256()
        count=0
        path=manifest.with_name(f'{word}-{name}.gz')
        with gzip.open(path,'rb') as f:
            for line in f:
                digest.update(line)
                key,w=line.split()
                key=int(key)
                if key<=previous:raise ValueError(f'Unordered raw prism bank {path}')
                previous=key
                count+=1
                yield key,i,fmpq(w.decode())
        row=record if name=='direct' else record['colors'][name]
        if count!=row['core_orbits'] or digest.hexdigest()!=row['weights_sha256']:
            raise ArithmeticError(f'Changed metric weights {path}')
    def rows():
        streams=[source(i,name) for i,name in enumerate(names)]
        for key,terms in groupby(merge(*streams,key=itemgetter(0)),key=itemgetter(0)):
            weights=[(i,w) for _,i,w in terms]
            combined=weights[-1][1] if weights[-1][0]==len(names)-1 else fmpq(0)
            if sum(w for i,w in weights if i<len(names)-1)!=combined:
                raise ArithmeticError(f'Individual colors differ from the combined numerator at {key}')
            yield key,1,weights
    def absorb(batch,aliases,periods):
        nonlocal union
        union+=len(batch)
        for key,cols in periods.items():check_positive(key,cols)
        for row,alias in zip(batch,aliases,strict=True):
            vals=periods[alias]
            for i,w in row[2]:
                for j,x in enumerate(vals):total[i][j]+=w*x
    if fast:
        from mixed_prism_fast_functional import integer_columns
        for batch in batches(rows(),32768):
            aliases,periods,denominators=integer_columns(scratch/'prism-index',[row[0] for row in batch])
            for key,c in periods.items():
                check_positive(key,tuple(fmpq(n,d) for n,d in zip(c,denominators,strict=True)))
            scale=lcm(*(set(int(w.q) for row in batch for _,w in row[2])))
            block=[[0,0,0] for _ in names]
            for row,key in zip(batch,aliases,strict=True):
                values=periods[key]
                for i,w in row[2]:
                    weight=int(w.p)*(scale//int(w.q))
                    for j,x in enumerate(values):block[i][j]+=weight*x
            for i,terms in enumerate(block):
                for j,(x,d) in enumerate(zip(terms,denominators,strict=True)):total[i][j]+=fmpq(x,scale*d)
            union+=len(batch)
            print('Verified raw prism powers',union,flush=True)
    else:
        contract(rows(),13,work,scratch,reference,tally=True,bulk=True,consume=absorb)
    if list(map(sum,zip(*total[:-1])))!=total[-1]:
        raise ArithmeticError('The sum of the integrated color periods differs from the combined prism')
    for tag,cols in zip(names,total,strict=True):
        entry=record if tag=='direct' else record['colors'][tag]
        entry.update(package(tuple(map(str,cols))))
    record['union_raw_powers']=union
    dest=manifest.with_name(word+'-period.json')
    pending=dest.with_suffix('.partial')
    pending.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    pending.replace(dest)
    print('All color periods, combined diagram and weight-five Ward coefficient',len(names)-1,
          *record['period'],record['ward_zeta5'],flush=True)


def complete_domain(work,scratch,reference,archive):
    """Evaluate and certify the full positive 26-dimensional scalar prism domain."""
    exe=work/'domain'
    compile_program('operator/mixed_prism_positive_domain.cpp',exe)
    keys=work/'angular-classes'
    with keys.open('wb') as out:
        census=run([exe,'13'],stdout=out,stderr=-1,text=True,check=True).stderr.strip()
    print(census,flush=True)
    counts=list(map(int,census.split(':')[-1].split()))
    if counts!=[67990131,5677268,2591953]:raise ArithmeticError('Changed positive prism census')
    archive.parent.mkdir(parents=True,exist_ok=True)
    count=0
    digest=sha256()
    best=None
    runner=None
    pending=archive.with_suffix('.partial')
    with keys.open() as inp,pending.open('wb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as output:
        def store(batch,aliases,periods):
            nonlocal count,best,runner
            for (key,_),alias in zip(batch,aliases,strict=True):
                if key!=alias:raise ArithmeticError('Noncanonical prism key')
                cols=periods[key]
                low,high=check_positive(key,cols)
                if best is None or low<best[0]:
                    runner=best
                    best=(low,high,key)
                elif runner is None or low<runner[0]:
                    runner=(low,high,key)
                line=(str(key)+' '+' '.join(map(str,cols))+'\n').encode()
                digest.update(line)
                output.write(line)
                count+=1
        contract(((int(line),1) for line in inp),13,work,scratch,reference,bulk=True,tally=True,consume=store)
    if count!=counts[-1]:raise ArithmeticError('Incomplete integral table')
    pending.replace(archive)
    print('All angular periods and uncompressed SHA-256:',count,digest.hexdigest(),flush=True)
    if best[1]<runner[0]:
        print('Unique minimum:',unpack(best[2],9),'strict bounds in units of 1e-60:',
              int(best[0].floor()),int(best[1].ceil()),flush=True)


def checks(work,scratch,reference):
    for powers, expected in zip(PROBES,zip(*FULL_PROBES),strict=True):
        result = contract([(pack(powers),fmpq(1))],12,work,scratch)
        assert tuple(map(str,result)) == tuple(map(str,expected))
        print('Independent prism propagator',powers,*result,flush=True)
    source = HERE/'mixed_bipartite_one_star.components.gz'
    with gzip.open(source,'rt') as f:
        rows = ((pack(tuple(map(int,s[:9]))),fmpq(s[9])) for line in f if (s:=line.split()))
        result = contract(rows,13,work,scratch,reference,tally=True,bulk=True)
    assert tuple(map(str,result)) == tuple(map(str,CHAMBER))
    print('Complete 37391-component physical control',*result,flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('scratch',type=Path,help='Writable scalar cache directory')
    parser.add_argument('banks',type=Path,nargs='*')
    parser.add_argument('--reference',type=Path,help='Existing project star13.sqlite, opened read-only')
    parser.add_argument('--raw',action='store_true',help='Apply product (a-1)! / (q-1)! to projective weights')
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--bulk',action='store_true',help='Pool equal terminal angular integrals across each 32768 input rows')
    parser.add_argument('--fast',action='store_true',help='Use the complete indexed archive for banks and color manifests')
    parser.add_argument('--colors',type=Path,help='Contract and verify all color banks in a metric weight manifest')
    parser.add_argument('--domain',type=Path,help='Write every positive q=13 completed scalar class to this gzip table')
    args = parser.parse_args()
    args.scratch.mkdir(parents=True,exist_ok=True)
    if (args.banks or args.colors) and not (args.check or args.domain) and (args.fast or
            (args.reference is None and ARCHIVE.exists() and METADATA.exists())):
        from mixed_prism_fast_functional import measure
        if args.colors:colors(args.colors,None,args.scratch,None,fast=True)
        for bank in args.banks:
            print('Integrated prism',bank,*measure(args.scratch/'prism-index',bank,args.raw),flush=True)
        raise SystemExit
    if args.reference is None and (args.banks or args.colors) and ARCHIVE.exists() and METADATA.exists() and not args.domain:
        args.reference=archive_reference(args.scratch)
    with TemporaryDirectory(dir=args.scratch) as tmp:
        work = Path(tmp)
        compile_program('quartic_prism_local_ibp.cpp',work/'ibp')
        compile_program('quartic_wheel_bank.cpp',work/'wheel')
        if args.check: checks(work,args.scratch,args.reference)
        if args.colors: colors(args.colors,work,args.scratch,args.reference)
        if args.domain: complete_domain(work,args.scratch,args.reference,args.domain)
        for bank in args.banks:
            result = contract(read_bank(bank,13,args.raw),13,work,args.scratch,args.reference,tally=True,bulk=args.bulk)
            print('Integrated prism',bank,*result,flush=True)
