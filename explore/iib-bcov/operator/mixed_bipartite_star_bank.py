"""Exact prism and wheel boundary data for the convergent bipartite identities."""
from argparse import ArgumentParser
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from contextlib import ExitStack
from itertools import islice
from pathlib import Path
from subprocess import PIPE,Popen,run
from tempfile import TemporaryDirectory
import gzip
import json
import sqlite3
import shutil

from flint import fmpq,fmpz
from mixed_bipartite_local_ibp import compile_program, reduced
from mixed_bipartite_fixed_spin import evaluate,verify_chart
from quartic_bipartite_replay import build
from quartic_bipartite_stars import preferred_star

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE/'mixed_bipartite_positive_periods.gz'
COMPLETION=HERE/'mixed_bipartite_positive_completion.gz'


def pack(powers):
    key=0
    for a in powers:key=(key<<4)|a
    return key


def unpack(key,n):
    return tuple((key>>(4*i))&15 for i in reversed(range(n)))


def expansion(key,q):
    a=unpack(key,9)
    if 0 in a:
        i,j=divmod(a.index(0),3)
        k,l=(i+1)%3,(i+2)%3
        m,n=(j+1)%3,(j+2)%3
        u=(a[3*k+m],a[3*l+m],a[3*l+n],a[3*k+n])
        b=(a[3*k+j],a[3*i+m],a[3*l+j],a[3*i+n])
        return 'wheel',[(pack(u+b),Fraction(1))]
    terms=preferred_star(key,q)
    return 'prism',[(pack(p),v) for p,v in terms]


def load_keys(conn,table,keys):
    # All table names are fixed by the caller. Bounded IN queries permit reuse
    # of the same angular database in a resumed or partial calculation.
    found={}
    for offset in range(0,len(keys),400):
        group=keys[offset:offset+400]
        marks=','.join('?' for _ in group)
        for key,*cols in conn.execute(f'SELECT key,r,w3,w5 FROM {table} WHERE key IN ({marks})',group):
            found[key]=tuple(map(Fraction,cols))
    return found


def fill(wanted,q,scratch,chunk):
    with TemporaryDirectory(dir=scratch) as directory,sqlite3.connect(scratch/f'star{q}.sqlite') as cache:
        tmp=Path(directory)
        ibp,wheel=tmp/'ibp',tmp/'wheel'
        compile_program('quartic_prism_local_ibp.cpp',ibp)
        compile_program('quartic_wheel_bank.cpp',wheel)
        cache.execute('CREATE TABLE IF NOT EXISTS wheels (degree INTEGER,key INTEGER,r TEXT,w3 TEXT,w5 TEXT,PRIMARY KEY(degree,key))')
        for kind in ('wheel','prism'):
            cache.execute(f'CREATE TABLE IF NOT EXISTS {kind} (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
            known=frozenset(row[0] for row in cache.execute(f'SELECT key FROM {kind}'))
            todo=sorted(wanted[kind]-known)
            for offset in range(0,len(todo),chunk):
                keys=todo[offset:offset+chunk]
                if kind=='prism':
                    cuts=[unpack(k,9) for k in keys]
                    periods=reduced(cuts,q,ibp,wheel,tmp,cache if q==13 else None,prime_workers=2)
                    values=[periods[i] for i in range(len(keys))]
                else:
                    for k in keys:verify_chart(k,q,0)
                    periods=evaluate(q,0,keys,wheel,tmp,prime_workers=2)
                    values=[periods[k] for k in keys]
                cache.executemany(f'INSERT INTO {kind} VALUES (?,?,?,?)',
                                  ((k,*map(str,v)) for k,v in zip(keys,values,strict=True)))
                cache.commit()
                print(kind,'verified',offset+len(keys),'of',len(todo),flush=True)


def canonical_cuts(executable,q,keys):
    labels={}
    iterator=iter(keys)
    while batch:=list(islice(iterator,50000)):
        query=''.join(f'{k}\n' for k in batch)
        response=run([executable,str(q),'--canonical-powers'],input=query,
                     capture_output=True,text=True,check=True)
        labels.update(zip(batch,map(int,response.stdout.splitlines()),strict=True))
    return labels


def merge_equal_prisms(wanted,q,scratch):
    # The original positive prism cut and every regrounding of its completed
    # six-point graph have equal periods in the projective normalization.
    with TemporaryDirectory(dir=scratch) as folder,sqlite3.connect(scratch/f'star{q}.sqlite') as db:
        exe=Path(folder)/'canonical'
        compile_program('quartic_prism_local_ibp.cpp',exe)
        aliases=canonical_cuts(exe,q,wanted)
        db.execute('CREATE TABLE IF NOT EXISTS prism (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
        db.execute('CREATE TEMP TABLE aliases (key INTEGER PRIMARY KEY,representative INTEGER)')
        stream=db.execute('SELECT key FROM prism')
        while rows:=stream.fetchmany(50000):
            batch=[k for k, in rows]
            extra=canonical_cuts(exe,q,[k for k in batch if k not in aliases])
            db.executemany('INSERT INTO aliases VALUES (?,?)',
                           ((k,aliases[k] if k in aliases else extra[k]) for k in batch))
        db.execute('CREATE TEMP TABLE common (key INTEGER PRIMARY KEY,r TEXT,w3 TEXT,w5 TEXT)')
        db.execute('INSERT OR IGNORE INTO common SELECT representative,r,w3,w5 FROM prism JOIN aliases USING(key)')
        mismatch=db.execute('SELECT 1 FROM prism AS p JOIN aliases AS a USING(key) JOIN common AS c '
                            'ON c.key=a.representative WHERE p.r!=c.r OR p.w3!=c.w3 OR p.w5!=c.w5 LIMIT 1')
        if mismatch.fetchone():raise ArithmeticError('Regrounded prism periods disagree')
        db.execute('INSERT OR IGNORE INTO prism SELECT * FROM common')
        db.commit()
        print('Exact cached charts and distinct regroundings:',db.execute('SELECT count(*) FROM aliases').fetchone()[0],
              db.execute('SELECT count(*) FROM common').fetchone()[0],flush=True)
    return aliases


def table(base,q,scratch,chunk,compact=False):
    scratch.mkdir(exist_ok=True)
    seeds=set(map(int,base.with_suffix('.wheel').read_text().split()))
    wanted={'prism':set(),'wheel':set()}
    for key in seeds:
        kind,terms=expansion(key,q)
        wanted[kind].update(k for k,c in terms)
    print('Boundary orbits, distinct prism cuts, distinct wheel cuts:',len(seeds),
          len(wanted['prism']),len(wanted['wheel']),flush=True)
    aliases={}
    if compact:
        aliases=merge_equal_prisms(wanted['prism'],q,scratch)
        wanted['prism']=set(aliases.values())
        print('Canonical prism cuts:',len(wanted['prism']),flush=True)
    fill(wanted,q,scratch,chunk)
    dest=base.with_suffix('.boundary')
    with sqlite3.connect(scratch/f'star{q}.sqlite') as cache:
        with base.with_suffix('.wheel').open() as source,dest.with_suffix('.partial').open('w') as output:
            while True:
                batch=list(map(int,islice(source,1000)))
                if not batch:break
                data=[expansion(k,q) for k in batch]
                if compact:
                    data=[(kind,[(aliases[v],c) for v,c in terms]) if kind=='prism' else (kind,terms)
                          for kind,terms in data]
                known={kind:load_keys(cache,kind,list({v for t,terms in data if t==kind for v,c in terms}))
                       for kind in ('prism','wheel')}
                for k,(kind,terms) in zip(batch,data,strict=True):
                    print(k,*(sum(c*known[kind][v][col] for v,c in terms) for col in range(3)),file=output)
        dest.with_suffix('.partial').replace(dest)


def measure(base,q,scratch,source,chunk):
    scratch.mkdir(exist_ok=True)
    total=[fmpq(0) for _ in range(3)]
    with gzip.open(source,'rt') as stream:
        weights={int(k):fmpq(v) for k,v in map(str.split,stream)}
    count=len(weights)
    def add(key,values):
        w=weights.pop(key,None)
        if w is not None:
            for i,v in enumerate(values):total[i]+=w*fmpq(v)
    if q==13 and ARCHIVE.is_file():
        verify_archive()
        metadata=json.loads(Path(__file__).with_suffix('.json').read_text())
        files=[ARCHIVE]
        if COMPLETION.is_file() and 'completion' in metadata:
            stats=metadata['completion']
            verify_member(COMPLETION,stats['extra_positive']+stats['extra_contracted'],stats['sha256'])
            files.append(COMPLETION)
        for file in files:
            with gzip.open(file,'rt') as period:
                for line in period:
                    if not weights:break
                    key,*v=line.split()
                    add(int(key),v)
    else:
        for i in range(3):
            with base.with_suffix(f'.channel{i}').open() as period:
                for line in period:
                    key,v=line.split();w=weights.get(int(key))
                    if w is not None:total[i]+=w*fmpq(v)
        with base.with_suffix('.nodes').open() as interior:
            for line in interior:weights.pop(int(line),None)
        with base.with_suffix('.boundary').open() as boundary:
            for line in boundary:
                key,*v=line.split()
                add(int(key),v)
    coefficients={'prism':{},'wheel':{}}
    for key,w in weights.items():
        kind,terms=expansion(key,q)
        coeff=coefficients[kind]
        for cut,c in terms:
            coeff[cut]=coeff.get(cut,0)+w*fmpq(c.numerator,c.denominator)
    wanted={kind:{k for k,c in cuts.items() if c} for kind,cuts in coefficients.items()}
    if any(wanted.values()):
        fill(wanted,q,scratch,chunk)
        with sqlite3.connect(scratch/f'star{q}.sqlite') as cache:
            for kind in coefficients:
                keys=list(wanted[kind])
                for offset in range(0,len(keys),1000):
                    batch=keys[offset:offset+1000]
                    values=load_keys(cache,kind,batch)
                    assert len(values)==len(batch)
                    for k,cols in values.items():
                        for i,v in enumerate(cols):
                            total[i]+=coefficients[kind][k]*fmpq(v.numerator,v.denominator)
    print('Core orbits, extra star expansions and physical period columns:',count,len(weights),*total,flush=True)
    return tuple(map(str,total))


def solve(base,scratch):
    base=base.resolve()
    boundary=base.with_suffix('.boundary')
    if not boundary.is_file():raise ValueError('The exact boundary values are required')
    with TemporaryDirectory(dir=scratch) as directory:
        tmp=Path(directory)
        rank=build(tmp,'quartic_bipartite_adjoint',('z',))
        cert=build(tmp,'quartic_bipartite_certificate',('gmpxx','gmp'))
        recipe=base.with_suffix('.recipe')
        if not recipe.exists():
            run([rank,'discover',base,'999999937',recipe],check=True)
        primes=[]
        candidate=999999938
        while True:
            for _ in range(9 if not primes else 4):
                candidate-=1
                while not fmpz(candidate).is_prime():candidate-=1
                prime=str(candidate)
                table=base.with_suffix(f'.{prime}.tsv')
                run([rank,'boundary',base,prime,recipe,boundary,table],check=True)
                primes.append((prime,table))
            manifest=tmp/'manifest'
            manifest.write_text(''.join(f'{p} {table}\n' for p,table in primes))
            certified=True
            for col in range(3):
                dest=base.with_suffix(f'.channel{col}')
                proc=run([cert,base,boundary,'-',manifest,str(col),dest])
                if proc.returncode==2:
                    certified=False
                    break
                proc.check_returncode()
            if certified:break


def prepare(base,q,scratch):
    scratch.mkdir(exist_ok=True)
    with TemporaryDirectory(dir=scratch) as folder:
        exe=Path(folder)/'orbit-ibp'
        compile_program('operator/mixed_bipartite_positive_ibp.cpp',exe)
        run([exe,str(q),base],input=b'',check=True)


def publish(base):
    expected=json.loads(Path(__file__).with_suffix('.json').read_text())
    for tag,wanted in expected['table_sha256'].items():
        with base.with_suffix('.'+tag).open('rb') as source:
            digest=sha256()
            for block in iter(lambda:source.read(1048576),b''):digest.update(block)
        assert digest.hexdigest()==wanted
    temp=ARCHIVE.with_suffix('.partial')
    with ExitStack() as files:
        raw=files.enter_context(temp.open('wb'))
        out=files.enter_context(gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6))
        channels=[files.enter_context(base.with_suffix(f'.channel{i}').open('rb')) for i in range(3)]
        for rows in zip(*channels,strict=True):
            entries=[r.split() for r in rows]
            assert len({r[0] for r in entries})==1
            out.write(b' '.join([entries[0][0],*(r[1] for r in entries)])+b'\n')
        bound=files.enter_context(base.with_suffix('.boundary').open('rb'))
        for block in iter(lambda:bound.read(1048576),b''):out.write(block)
    temp.replace(ARCHIVE)


@lru_cache(maxsize=4)
def checked_member(path,identity,rows,expected_digest):
    digest=sha256()
    count=0
    with gzip.open(path,'rb') as source:
        for block in iter(lambda:source.read(1048576),b''):
            digest.update(block);count+=block.count(b'\n')
    if count!=rows or digest.hexdigest()!=expected_digest:
        raise ArithmeticError('Period archive failed its row count or checksum')


def verify_member(path,rows,digest):
    stat=path.stat()
    checked_member(path,(stat.st_dev,stat.st_ino,stat.st_size,stat.st_mtime_ns,stat.st_ctime_ns),rows,digest)


def verify_archive():
    expected=json.loads(Path(__file__).with_suffix('.json').read_text())
    verify_member(ARCHIVE,expected['interior_orbits']+expected['boundary_orbits'],expected['archive_sha256'])


def combine_columns(columns,target):
    """Preserve the earlier archive and require exact agreement on its overlap."""
    old={}
    with gzip.open(ARCHIVE,'rb') as source:
        for line in source:
            key,*values=line.split()
            old[int(key)]=values
    digest=sha256()
    extra_positive=extra_scalar=0
    with ExitStack() as stack:
        inputs=[stack.enter_context(path.open('rb')) for path in columns]
        raw=stack.enter_context(target.open('wb'))
        out=stack.enter_context(gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=6))
        for records in zip(*inputs,strict=True):
            split=[r.split() for r in records]
            keys=[r[0] for r in split]
            if not keys[0]==keys[1]==keys[2]:raise ValueError('Unequal period-channel ordering')
            key=int(keys[0])
            values=[r[1] for r in split]
            earlier=old.pop(key,None)
            if earlier is not None:
                if earlier!=values:raise ArithmeticError(f'Changed established period at {key}')
                continue
            line=b' '.join([keys[0],*values])+b'\n'
            digest.update(line);out.write(line)
            # Every contracted orbit has a leading zero in its canonical nine nibbles.
            if key>>32:extra_positive+=1
            else:extra_scalar+=1
    if old:raise ValueError('Incomplete support of the original exact archive')
    return {'extra_positive':extra_positive,'extra_contracted':extra_scalar,'sha256':digest.hexdigest()}


def verify_combined(executable,supplement):
    # Independent enumeration bypasses the triangular recipe, chart grouping
    # and binary equation list used to reconstruct the three columns.
    with Popen([executable,'13','--verify-positive'],stdin=PIPE,stdout=PIPE) as child:
        for path in (ARCHIVE,supplement):
            with gzip.open(path,'rb') as source:shutil.copyfileobj(source,child.stdin)
        child.stdin.close()
        summary=child.stdout.read().decode().strip()
        if child.wait():raise ArithmeticError('Full positive-orbit verification failed')
    return summary


def complete_domain(scratch):
    """Reconstruct the full q=13 domain using the committed interior/boundary table."""
    scratch.mkdir(exist_ok=True,parents=True)
    verify_archive()
    with TemporaryDirectory(dir=scratch) as directory:
        tmp=Path(directory)
        full=tmp/'full';recipe=tmp/'recipe';labels=tmp/'scalar-labels'
        orbit=tmp/'orbits';closure=tmp/'closure';wheel=tmp/'wheel'
        compile_program('quartic_bipartite_orbit_ibp.cpp',orbit)
        compile_program('operator/mixed_bipartite_positive_ibp.cpp',closure)
        compile_program('quartic_wheel_bank.cpp',wheel)
        run([orbit,'13',full],check=True)
        with full.with_suffix('.wheel').open('rb') as source,labels.open('wb') as output:
            run([closure,'13','--contact-map'],stdin=source,stdout=output,check=True)
        keys=tmp/'anchors.keys'
        anchors=tmp/'anchors.tsv'
        with gzip.open(ARCHIVE,'rb') as src,keys.open('wb') as key,anchors.open('wb') as table_file:
            for line in src:
                key.write(line.split(maxsplit=1)[0]+b'\n');table_file.write(line)
        frontier=tmp/'extra'
        with keys.open('rb') as key:
            run([closure,'frontier',full,labels,frontier],stdin=key,check=True)
        table(frontier,13,scratch,5000,compact=True)
        with frontier.with_suffix('.boundary').open('rb') as src,keys.open('ab') as key,anchors.open('ab') as data:
            for row in src:
                key.write(row.split(maxsplit=1)[0]+b'\n');data.write(row)
        with keys.open('rb') as key:
            run([closure,'plan',full,labels,recipe],stdin=key,check=True)
        needed=set()
        for suffix in ('.angular','.unused'):
            needed.update(map(int,recipe.with_suffix(suffix).read_text().split()))
        for key in needed:verify_chart(key,13,0)
        angular=evaluate(13,0,needed,wheel,tmp,prime_workers=2)
        direct=tmp/'angular.tsv'
        with direct.open('w') as source:
            for k,v in sorted(angular.items()):print(k,*v,file=source)
        channels=[tmp/f'channel{c}' for c in range(3)]
        for c,path in enumerate(channels):
            with anchors.open('rb') as data:
                run([closure,'resolve',full,labels,recipe,str(c),direct,path],stdin=data,check=True)
        candidate=tmp/COMPLETION.name
        stats=combine_columns(channels,candidate)
        print(verify_combined(closure,candidate),flush=True)
        expected=json.loads(Path(__file__).with_suffix('.json').read_text())['completion']
        assert all(expected[k]==v for k,v in stats.items())
        candidate.replace(scratch/COMPLETION.name)
        print('Reproduced full-domain supplemental archive:',stats,flush=True)


def check_domain_signs():
    """Certify Euclidean positivity of both complete sets, including severe cancellation."""
    from math import factorial
    from quartic_prism_physical_ibp import zeta_interval
    verify_archive()
    expected=json.loads(Path(__file__).with_suffix('.json').read_text())['completion']
    verify_member(COMPLETION,expected['extra_positive']+expected['extra_contracted'],expected['sha256'])
    scale=10**60
    # On the whitened unit momentum sphere the nine edge quadrics sum to one.
    # The transitive equal-power minimum and AM-GM give 81^q Gamma(q)/Gamma(4q).
    spherical_bound=-((-scale*81**13*factorial(12))//factorial(51))
    intervals=[]
    for power in (3,5):
        lo,hi=zeta_interval(power,192)
        intervals.append(((scale*lo)//1,-((-scale*hi)//1)))
    counts=[0,0]
    smallest=[None,None]
    next_smallest=[None,None]
    least_keys=[[],[]]
    least_columns=[None,None]
    upper=[None,None]
    for path in (ARCHIVE,COMPLETION):
        with gzip.open(path,'rb') as stream:
            for row in stream:
                key,*columns=row.split()
                key=int(key)
                group=int(any(((key>>(4*j))&15)==0 for j in range(9)))
                nums=[]
                dens=[]
                for c in columns:
                    n,_,d=c.partition(b'/')
                    nums.append(int(n))
                    dens.append(int(d or 1))
                r=scale*nums[0]*dens[1]*dens[2]
                z3=nums[1]*dens[0]*dens[2]
                z5=nums[2]*dens[0]*dens[1]
                denominator=dens[0]*dens[1]*dens[2]
                lo=(r+z3*intervals[0][z3<0]+z5*intervals[1][z5<0])//denominator
                if lo<spherical_bound:raise ArithmeticError(f'Spherical lower bound failed at {key}')
                counts[group]+=1
                if smallest[group] is None or lo<smallest[group]:
                    next_smallest[group]=smallest[group]
                    smallest[group]=lo
                    least_keys[group]=[key]
                    least_columns[group]=columns
                    upper[group]=-((-(r+z3*intervals[0][z3>=0]+z5*intervals[1][z5>=0]))//denominator)
                elif columns==least_columns[group]:
                    least_keys[group].append(key)
                elif next_smallest[group] is None or lo<next_smallest[group]:
                    next_smallest[group]=lo
    assert counts==[1386282,972294]
    for group in range(2):
        if upper[group]>=next_smallest[group]:
            raise ArithmeticError('Extremal scalar periods have overlapping enclosures')
        print('Strictly positive scalar periods, minimizing powers, minimum enclosure in units of 10^-60',
              counts[group],[unpack(k,9) for k in sorted(least_keys[group])],smallest[group],upper[group],flush=True)


def control(scratch):
    from quartic_bipartite_physical import matrix_orderings
    from quartic_bipartite_position_ibp import orbits
    from quartic_harmonic_period import period
    base=scratch/'matching4'
    prepare(base,4,scratch)
    table(base,4,scratch,1000)
    solve(base,scratch)
    weights=scratch/'matching4.all.gz'
    expected=[Fraction(),Fraction(),Fraction()]
    all_rows=orbits(4)
    assert len(all_rows)==35
    with gzip.open(weights,'wt') as out:
        for i,edge_powers in enumerate(all_rows):
            key=min(pack(tuple(edge_powers[j] for j in o)) for o in matrix_orderings())
            coefficient=(-1)**i*(i*i+1)
            print(key,coefficient,file=out)
            direct=period(edge_powers,4,i%2)
            for j,w in enumerate(((),(3,),(5,))):
                expected[j]+=coefficient*Fraction(str(direct.get(w,0)))
    measured=measure(base,4,scratch,weights,1000)
    assert measured==tuple(map(str,expected))==('85487/320','-21965/48','6793/24')
    print('All 35 independent q=4 harmonic integrals passed the signed contraction.')


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('base',type=Path)
    parser.add_argument('q',type=int)
    parser.add_argument('scratch',type=Path)
    parser.add_argument('--chunk',default=10000,type=int)
    parser.add_argument('--prepare',action='store_true',help='generate the finite positive-orbit incidence system')
    parser.add_argument('--check',action='store_true',help='replay the independent q=4 harmonic integral controls')
    parser.add_argument('--archive',action='store_true',help='pack the verified q=13 exact period tables')
    parser.add_argument('--check-domain-signs',action='store_true',help='certify positivity of both full scalar domains')
    parser.add_argument('--complete-domain',action='store_true',help='reconstruct and check every q=13 positive period')
    parser.add_argument('--canonical-stars',action='store_true',help='Group prism boundary cuts by conformal regrounding')
    parser.add_argument('--interior',action='store_true',help='solve from the completed exact boundary table')
    parser.add_argument('--weights',type=Path,help='contract the three verified channels with a physical orbit bank')
    options=parser.parse_args()
    if options.check_domain_signs:
        if options.q!=13:parser.error('the positive archive uses half-dimension 13')
        check_domain_signs()
    elif options.complete_domain:
        if options.q!=13:parser.error('the complete-domain archive uses half-dimension 13')
        complete_domain(options.scratch)
    elif options.archive:
        publish(options.base)
        verify_archive()
    elif options.check:
        control(options.scratch)
    elif options.prepare:
        prepare(options.base,options.q,options.scratch)
    elif options.weights:
        measure(options.base,options.q,options.scratch,options.weights,options.chunk)
    elif options.interior:
        solve(options.base,options.scratch)
    else:
        table(options.base,options.q,options.scratch,options.chunk,options.canonical_stars)
