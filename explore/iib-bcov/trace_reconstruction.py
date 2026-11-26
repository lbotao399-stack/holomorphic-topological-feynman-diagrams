"""Recover the three-loop primary counterterm from seven remaining graph projections."""
from argparse import ArgumentParser
from fractions import Fraction as F
from itertools import chain,combinations
from math import comb
from pathlib import Path
from subprocess import run,Popen,PIPE
from tempfile import TemporaryDirectory
import json
import sys
from flint import fmpz,fmpq_mat,nmod_mat
from fourfold_bulk import CORES,symmetries

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'operator'))
from sparse_trace_rows import fixtures as sparse
from three_loop_soft_jet import exact_boundary
from bulk_three_hexagon import potentials

TESTS={'four':((2,3,0),(2,4,0),(2,5,0),(2,6,0)),
       'three':((2,3,0),(2,4,0)), 'two':((2,3,0),)}


def waves():
    p=[tuple(a*int(k==j) for k in range(5)) for j,m in enumerate((2,2,2,2,1))
       for a in ((1,j+2) if m==2 else (1,))]
    p.append(tuple(-sum(x) for x in zip(*p)))
    return dict(zip(TESTS,(p,*sparse())))


def system(soft=False):
    data=json.loads((ROOT/'operator'/'three_loop_constraints.json').read_text())
    rows=data['rows'][:]
    rhs=[list(map(F,r)) for r in data['reduced_periods']]
    with TemporaryDirectory() as td:
        exe=Path(td)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for name,p in waves().items():
            if soft and name=='two':
                from soft_trace_measurement import trace_row, hexagon_row
                row=trace_row(exe)
                hexagon=hexagon_row()
                stored=json.loads((ROOT/'operator'/'soft_trace_measurement.json').read_text())
                assert row==stored['row'] and str(hexagon)==stored['hexagon']
                rows.append(row)
                rhs.append([hexagon,F(0)])
                continue
            arr=exact_boundary(exe,p)
            hx=potentials(p)
            for i,j,k in TESTS[name]:
                assert not p[i][k] and not p[j][k]
                rows.append([a[i][k]-a[j][k] for a in arr])
                rhs.append([hx[i][k]-hx[j][k],F(0)])
    assert int(nmod_mat(rows,1009).det())==(515 if soft else 338)
    return rows,rhs


def load_rule(path,core,prime):
    lines=path.read_text().splitlines()
    if lines[0]!=f'# {core} U power 12 prime {prime} times and two weights':
        raise ValueError('Cubature header differs from the requested core or prime.')
    rule=[tuple(map(int,line.split())) for line in lines[1:]]
    if len(rule)!={'tetrahedron':2540,'double-square':841}[core] or any(len(r)!=8 for r in rule):
        raise ValueError('Incomplete core cubature.')
    return rule


def integrate(exe,core,name,prime,root,start=0,stop=None,*,tests=None):
    import fcntl
    rule=load_rule(root/f'core-rule-{core}-{prime}.tsv',core,prime)
    if tests is None:tests=TESTS[name]
    file=root/f'primary-{core}-{name}-{prime}.tsv'
    args=[str(exe),str(prime),core,name,*map(str,chain.from_iterable(tests))]
    # Store the node coordinates with every completed raw numerator, so a rule
    # can be regenerated or its nodes reordered without reusing a wrong value.
    with file.open('a+') as ledger:
        fcntl.flock(ledger,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ledger.seek(0)
        saved={}
        for line in ledger:
            vals=tuple(map(int,line.split()))
            if len(vals)!=6+len(tests) or vals[:6] in saved:
                raise ValueError('Bad or repeated numerator checkpoint.')
            saved[vals[:6]]=vals[6:]
        pending=[r[:6] for r in rule[start:stop] if r[:6] not in saved]
        if pending:
            with Popen(args,stdin=PIPE,stdout=PIPE,text=True) as worker:
                try:
                    for node in pending:
                        print(*node,file=worker.stdin,flush=True)
                        record=worker.stdout.readline()
                        vals=tuple(map(int,record.split()))
                        if len(vals)!=len(tests) or not all(0<=v<prime for v in vals):
                            raise ArithmeticError('Oracle did not return the requested transverse numerators.')
                        saved[node]=vals
                        print(*node,*vals,sep='\t',file=ledger,flush=True)
                finally:
                    worker.stdin.close()
                    if worker.wait():
                        raise RuntimeError('The core numerator process failed.')
        if len(saved)<len(rule):return None
        if set(saved)!={r[:6] for r in rule}:raise ValueError('Checkpoint nodes do not match the rule.')
    inverse=pow(len(symmetries(CORES[core])),-1,prime)
    return [[sum(saved[r[:6]][i]*r[6+k] for r in rule)*inverse%prime for k in (0,1)]
            for i in range(len(tests))]


def reconstruct(exe,core,name,root,*,tests=None,height_source=None,limits=None):
    from core_cubature import rational_table,rule
    from period_bounds import bounds,core_bounds
    graph=CORES[core]
    table_file=root/f'core-periods-{core}.tsv'
    moments=rational_table(graph,12,table_file)
    if tests is None:tests=TESTS[name]
    if limits is None:
        if height_source is None:height_source=waves()[name]
        if name=='three':
            from ordered_period_bounds import ordered_bounds
            limits=ordered_bounds(core,table_file,height_source,tests)
        else:
            limits=bounds(height_source,tests,graph,core_bounds(graph,table_file))
    nums=[[0,0] for q in limits]
    modulus=1
    candidate=999983
    while modulus<=2*max(bound for row in limits for den,bound in row):
        while not fmpz(candidate).is_prime():candidate-=2
        prime=candidate
        candidate-=2
        file=root/f'core-rule-{core}-{prime}.tsv'
        if not file.exists():
            values=rule(graph,12,prime,moments)
            # Unique per-process temporary name permits concurrent callers to
            # generate the same deterministic rule safely.
            from tempfile import NamedTemporaryFile
            with NamedTemporaryFile('w',dir=root,delete=False) as out:
                print('#',core,'U power',12,'prime',prime,'times and two weights',file=out)
                for row in values:print(*row,sep='\t',file=out)
                tmp=Path(out.name)
            tmp.replace(file)
        result=integrate(exe,core,name,prime,root,tests=tests)
        for i,row in enumerate(result):
            for k,v in enumerate(row):
                den,bound=limits[i][k]
                new=(den*v-nums[i][k])*pow(modulus,-1,prime)%prime
                nums[i][k]+=modulus*new
        modulus*=prime
        print(core,name,'prime',prime,'CRT bits',modulus.bit_length(),flush=True)
    values=[]
    for i,row in enumerate(nums):
        out=[]
        for k,v in enumerate(row):
            den,bound=limits[i][k]
            v=v-modulus if 2*v>modulus else v
            assert abs(v)<=bound
            out.append(F(v,den))
        values.append(out)
    # Use a separate small prime for the entire projection outside the CRT.
    check=integrate(exe,core,name,1009,root,tests=tests)
    assert check==[[v.numerator*pow(v.denominator,-1,1009)%1009 for v in row] for row in values]
    record=dict(core=core,fixture=name,probes=tests,verification_prime=1009,
                modulus=str(modulus),bounds=limits,values=[[str(v) for v in row] for row in values])
    file=root/f'primary-{core}-{name}-exact.json'
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile('w',dir=root,delete=False) as out:
        json.dump(record,out,indent=2)
        print(file=out)
        tmp=Path(out.name)
    tmp.replace(file)
    return values


def checked_periods(directory,core,name,tests):
    from core_cubature import residue
    records=[]
    for method in ('exact','polynomial'):
        path=directory/f'primary-{core}-{name}-{method}.json'
        if not path.exists():continue
        record=json.loads(path.read_text())
        assert record['core']==core and record['fixture']==name and record['verification_prime']==1009
        assert record['probes']==[list(t) for t in tests] and len(record['values'])==len(tests)
        assert all(len(row)==2 for row in record['values'])
        if method=='polynomial':
            from core_numerator_reconstruction import verify_periods
            verify_periods(directory,core,name,record)
        else:
            assert len(record['bounds'])==len(tests) and all(len(r)==2 for r in record['bounds'])
            assert int(record['modulus'])>2*max(b for r in record['bounds'] for d,b in r)
            for period,limits in zip(record['values'],record['bounds']):
                for x,(d,b) in zip(map(F,period),limits):
                    assert d%x.denominator==0 and abs(x)*d<=b
        records.append(record['values'])
    if not records:raise FileNotFoundError(f'No exact {name} measurement on {core}.')
    assert all(values==records[0] for values in records)
    check=integrate(None,core,name,1009,directory,0,0,tests=tests)
    assert check==[[residue(F(x),1009) for x in row] for row in records[0]]
    return records[0]


def measured_rows(directory, name, rows, rhs):
    family='two' if name=='soft-two' else name
    start=11
    for key,tests in TESTS.items():
        if key==family:break
        start+=len(tests)
    tests=TESTS[family]
    stop=start+len(tests)
    values=[row[:] for row in rhs[start:stop]]
    cores={}
    for core in ('tetrahedron','double-square'):
        cores[core]=checked_periods(directory,core,name,tests)
        for row,period in zip(values,cores[core]):
            for k,x in enumerate(period):row[k]+=F(x)
    return dict(fixture=name, probes=tests, rows=rows[start:stop], cores=cores,
                hexagon=[str(r[0]) for r in rhs[start:stop]],
                reduced_periods=[[str(v) for v in row] for row in values])


def solve_complete(directory, soft=False):
    rows,lower=system(soft)
    rhs=lower[:11]
    for name in TESTS:
        if soft and name=='two':name='soft-two'
        rhs.extend(measured_rows(directory,name,rows,lower)['reduced_periods'])
    matrix=fmpq_mat(rows)
    target=fmpq_mat([[str(x) for x in row] for row in rhs])
    answer=matrix.solve(target)
    assert matrix*answer==target
    return [[str(x) for x in row] for row in answer.tolist()]


def solve_reserved(directory):
    from core_cubature import residue
    prime=1009
    rows,lower=system(soft=True)
    rhs=[[residue(v,prime) for v in pair] for pair in lower]
    cursor=11
    for family,tests in TESTS.items():
        name='soft-two' if family=='two' else family
        for core in ('tetrahedron','double-square'):
            # An empty node range reads existing full cubatures and launches no worker.
            period=integrate(None,core,name,prime,directory,0,0,tests=tests)
            if period is None:raise ValueError(f'Incomplete {name} cubature on {core}.')
            for i,(r,s) in enumerate(period,cursor):
                rhs[i]=[(rhs[i][0]+r)%prime,(rhs[i][1]+s)%prime]
        cursor+=len(tests)
    matrix=nmod_mat(rows,prime)
    target=nmod_mat(rhs,prime)
    answer=matrix.solve(target)
    assert matrix*answer==target
    coefficients=[[int(x) for x in row] for row in answer.tolist()]
    assert all(a and b for a,b in coefficients)
    assert all((a*d-b*c)%prime for (a,b),(c,d) in combinations(coefficients,2))
    return dict(prime=prime, basis='three_loop_soft_jet.DOUBLE_SOFT',
                determinant=int(matrix.det()), rhs=rhs,
                rational_and_zeta=coefficients)


def coefficient_certificate(directory):
    from core_cubature import residue
    from three_loop_soft_jet import DOUBLE_SOFT
    rows,lower=system(soft=True)
    rhs=[[str(x) for x in row] for row in lower[:11]]
    measurements=[]
    for name in ('four','three','soft-two'):
        data=measured_rows(directory,name,rows,lower)
        rhs.extend(data['reduced_periods'])
        measurements.append(data)
    matrix=fmpq_mat(rows)
    target=fmpq_mat(rhs)
    answer=matrix.solve(target)
    assert matrix*answer==target
    result=[[str(x) for x in row] for row in answer.tolist()]
    reserved=json.loads((ROOT/'operator'/'three_loop_reserved.json').read_text())
    for key,actual in (('rhs',rhs),('rational_and_zeta',result)):
        assert [[residue(F(x),1009) for x in row] for row in actual]==reserved[key]
    assert int(nmod_mat(rows,1009).det())==reserved['determinant']==515
    # Euler's transform of eta(3). The kth finite difference equals the
    # integral of (1-x)^k log(x)^2/2 on [0,1], hence lies in (0,1].
    lower=F(0)
    k=0
    while True:
        difference=sum((-1)**j*F(comb(k,j),(j+1)**3) for j in range(k+1))
        lower+=F(4,3*2**(k+1))*difference
        k+=1
        upper=lower+F(4,3*2**k)
        pairs=[[F(a)+F(b)*z for z in (lower,upper)] for a,b in result]
        if all(x*y>0 for x,y in pairs):break
    signs=[1 if x>0 else -1 for x,y in pairs]
    factors={}
    for core in ('tetrahedron','double-square'):
        if all((directory/f'primary-{core}-{name}-polynomial.json').exists()
               for name in ('four','three','soft-two')):
            from core_numerator_reconstruction import common_core_factor
            factors[core]=common_core_factor(directory,core)
    return dict(basis=reserved['basis'],basis_words=DOUBLE_SOFT,
                sign='C in Phi=C; the bulk action contains -C', verification_prime=1009,
                determinant=str(matrix.det()),rows=rows,measurements=measurements,core_common_factors=factors,
                rhs=rhs,rational_and_zeta=result,coefficient_signs=signs,
                zeta3_euler_terms=k,zeta3_interval=[str(lower),str(upper)])


if __name__=='__main__':
    cli=ArgumentParser()
    sub=cli.add_subparsers(dest='mode',required=True)
    sub.add_parser('system')
    solve=sub.add_parser('solve')
    solve.add_argument('directory',type=Path)
    solve.add_argument('--soft',action='store_true')
    reserved=sub.add_parser('reserved')
    reserved.add_argument('directory',type=Path)
    certificate=sub.add_parser('certificate')
    certificate.add_argument('directory',type=Path)
    measurements=sub.add_parser('measurements')
    measurements.add_argument('fixture',choices=(*TESTS,'soft-two'))
    measurements.add_argument('directory',type=Path)
    crt=sub.add_parser('reconstruct')
    crt.add_argument('exe',type=Path)
    crt.add_argument('core',choices=('tetrahedron','double-square'))
    crt.add_argument('fixture',choices=TESTS)
    crt.add_argument('directory',type=Path)
    case=sub.add_parser('integrate')
    case.add_argument('exe',type=Path)
    case.add_argument('core',choices=('tetrahedron','double-square'))
    case.add_argument('fixture',choices=TESTS)
    case.add_argument('prime',type=int)
    case.add_argument('directory',type=Path)
    case.add_argument('--start',type=int,default=0)
    case.add_argument('--stop',type=int)
    args=cli.parse_args()
    if args.mode=='system':
        m,h=system()
        print('18 independent integer rows, determinant 338 modulo 1009; hexagons:',
              *(str(x[0]) for x in h[11:]))
    elif args.mode=='solve':
        print(json.dumps(solve_complete(args.directory,args.soft),indent=2))
    elif args.mode=='reserved':
        print(json.dumps(solve_reserved(args.directory),indent=2))
    elif args.mode=='certificate':
        print(json.dumps(coefficient_certificate(args.directory),indent=2))
    elif args.mode=='measurements':
        print(json.dumps(measured_rows(args.directory,args.fixture,*system(soft=args.fixture=='soft-two')),indent=2))
    elif args.mode=='reconstruct':
        print(reconstruct(args.exe,args.core,args.fixture,args.directory))
    else:
        print(args.core,args.fixture,args.prime,
              integrate(args.exe,args.core,args.fixture,args.prime,args.directory,args.start,args.stop))
