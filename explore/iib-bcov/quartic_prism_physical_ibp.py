"""Reconstruct the complete physical prism weight-five pair; all caches are temporary."""
from argparse import ArgumentParser
from fractions import Fraction
from gzip import open as gzip_open
from io import BytesIO
from math import comb, factorial, prod
import json
from pathlib import Path
from shutil import copyfileobj
from subprocess import PIPE, Popen, run
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
TOTAL = (Fraction(-13914257373753679,45754008300000),
         Fraction(-1633443038914372549,1045387581638400000))
CONTROL = (Fraction(-7741325,27216),Fraction(-150875,6048))
PROBES = ((5,5,5,5,5,5,6,6,6),(4,6,5,6,5,5,6,6,5))
VALUES = (Fraction(575,92681981263872),Fraction(227,27259406254080))

FULL_TOTAL = (
    (Fraction(43896298966379923306171160083,586747539003225600000000000),
     Fraction(476340584507847676272946020919,1239210802374812467200000000000)),
    (Fraction(768730384185255903297563,3841799362521120000000),
     Fraction(758352372011658482516989,737625477604055040000000)),
    TOTAL,
)
FULL_PROBES = (
    (Fraction(-7039683385921,4612448552370359500800000),
     Fraction(-80381826781559,39237843587872849920000000)),
    (Fraction(-427123877,104634249567660933120),Fraction(-31992971,5838964819623936000)),
    VALUES,
)


def compile_program(source,destination,extra_libraries=()):
    # Homebrew installs GMP outside Clang's default system include paths.
    flags = []
    for prefix in (Path('/opt/homebrew'),Path('/usr/local')):
        if (prefix/'include/gmpxx.h').is_file():
            flags = ['-I'+str(prefix/'include'),'-L'+str(prefix/'lib')]
            break
    run(['c++','-std=c++17','-O3',*flags,str(HERE/source),'-lgmpxx','-lgmp',*extra_libraries,
         '-o',str(destination)],check=True)


def stream(args,source,output=PIPE):
    with Popen(args,stdin=PIPE,stdout=output) as job:
        copyfileobj(source,job.stdin)
        job.stdin.close()
        text = job.stdout.read().decode() if job.stdout else ''
        if job.wait():
            raise RuntimeError('Prism replay failed')
    return text


def sums(text,expected,rows):
    line = next(s for s in text.splitlines() if s.startswith('weight sums '))
    _,_,count,nonzero,k,j = line.split()
    assert int(count) == rows and (Fraction(k),Fraction(j)) == expected
    print(line,flush=True)
    return int(nonzero)


def check_forcer(form,forcer):
    # The native reconstruction needs only C++17/GMP. The optional independent
    # propagator run uses this repository's existing Python-FLINT/FORCER setup.
    from quartic_master_weights import two_cuts
    import sympy as sp
    with TemporaryDirectory() as scratch:
        tmp = Path(scratch)
        for name in ('forcer','forcer.h'):
            (tmp/name).symlink_to((forcer/name).resolve(strict=True))
        for a,columns in zip(PROBES,zip(*FULL_PROBES)):
            result = two_cuts(tmp,form.resolve(strict=True),'prism',12,{a:Fraction(1)})
            z3,z5=sp.symbols('z3 z5')
            target=columns[0]+columns[1]*z3+columns[2]*z5
            assert sp.expand(result-target)==0
            print('Both cuts:',a,result,flush=True)


def check(form=None,forcer=None):
    with TemporaryDirectory() as scratch:
        tmp = Path(scratch)
        native,reflect = tmp/'ibp',tmp/'reflection'
        compile_program('quartic_prism_physical_ibp.cpp',native)
        compile_program('quartic_prism_completion.cpp',reflect)
        for q,extra in ((6,1200),(12,1800)):
            start = [native,str(q)]
            pivots,rows,table = (tmp/(s+str(q)) for s in ('triangular','rows','table'))
            run([*start,pivots],check=True)
            run([*start,pivots,rows,str(extra)],check=True)
            solve = [*start,pivots,rows,'--solve',table]
            if q == 6:
                bank = json.loads((HERE.parent/'nonSP/verify/prism_forcer/periods.json').read_text())
                inputs = '\n'.join(' '.join(map(str,(*[v+1 for v in ex],i-740,i*i%97-48)))
                                   for i,ex in enumerate(bank['exponents']))
                output = stream(solve,BytesIO(inputs.encode()))
                sums(output,CONTROL,1481)
                continue
            with gzip_open(HERE/'quartic_prism.weights.gz','rb') as source:
                output = stream(solve,source)
            assert 'exact pass 6 nonzero 267543 max numerator and denominator bits 62 83' in output
            assert sums(output,TOTAL,1598657) == 1282191
            run([*start,'--check-all',table],check=True)
            reduced = tmp/'reduced'
            with gzip_open(HERE/'quartic_prism.weights.gz','rb') as source, reduced.open('wb') as out:
                stream([reflect,'12','1925','--reflected-weights'],source,out)
            with reduced.open('rb') as source:
                assert sums(stream([*start,'--sum',table],source),TOTAL,166419) == 163316
            probes = '\n'.join(' '.join(map(str,(*a,int(i==0),int(i==1)))) for i,a in enumerate(PROBES))
            sums(stream([*start,'--sum',table],BytesIO(probes.encode())),VALUES,2)
    if form is not None:
        check_forcer(form,forcer)


def zeta_interval(power, terms=64):
    # Euler transformation of the alternating zeta series. The finite
    # differences are moments of (1-t)^j times a positive probability density.
    factor = 1-Fraction(1,2**(power-1))
    lower = sum(sum((-1)**k*Fraction(comb(j,k),(k+1)**power) for k in range(j+1))
                / (2**(j+1)*factor) for j in range(terms))
    return lower, lower+1/(2**terms*factor)


def check_signs():
    z3,z5 = zeta_interval(3),zeta_interval(5)
    intervals = ((Fraction(-2382,10**12),Fraction(-2381,10**12)),
                 (Fraction(-1043,10**14),Fraction(-1042,10**14)))
    for (a,b,c),(lo,hi) in zip(zip(*FULL_TOTAL),intervals):
        assert b>0 and c<0
        assert lo < a+b*z3[0]+c*z5[1] <= a+b*z3[1]+c*z5[0] < hi


def ordinary_zeta_interval(power, n=2*10**6, scale=10**40):
    # Convexity bounds the tail below by the integral from n+1 plus its
    # half endpoint, and above by integrating the cells centered at k>n.
    partial = sum(scale//k**power for k in range(1,n+1))
    lower = Fraction(1,(power-1)*(n+1)**(power-1))+Fraction(1,2*(n+1)**power)
    upper = Fraction(2**(power-1),(power-1)*(2*n+1)**(power-1))
    return Fraction(partial,scale)+lower,Fraction(partial+n,scale)+upper


def physical_response_bounds():
    """Bounds for the full elementary probe from 39 exact shapes and 11 signed shapes."""
    z3,z5 = ordinary_zeta_interval(3),ordinary_zeta_interval(5)
    a,b,c = (col[1] for col in FULL_TOTAL)
    a += -Fraction(1,6490767360000)+Fraction(197793221359,36248733480960000000)
    b -= Fraction(507937,111878807040000)
    assert b>0 and c<0
    j39_upper = a+b*z3[1]+c*z5[0]
    certificate = json.loads((HERE/'quartic_primitive_moments.json').read_text())
    k33_lower = Fraction(certificate['K33']['outward_minus_J_times_1e20'][0],10**20)
    minus_j_lower = -j39_upper+k33_lower
    # 141 P1 - 56 P2 = 8923136 J4 for the two integer measurements.
    factors = (Fraction(9072,11),Fraction(8923136,197),432,1728)
    lower = tuple(int(10**18*factor*minus_j_lower) for factor in factors)
    assert lower == (9435927830,518233567709,4942628863,19770515453)
    print('Full lower bounds, scaled by 10^18: complex absolute, integer maximum,'
          ' real difference, minus 1728 J4:',*lower,flush=True)
    return lower


def verify_physical_ranges(total):
    """Independent ordinary-series signs for the complete three-column reduction."""
    zetas=((Fraction(1),Fraction(1)),ordinary_zeta_interval(3),ordinary_zeta_interval(5))
    combinations=(
        ('K4',1,0,-310607099214), ('J4',0,1,-1383499096),
        ('B slope',Fraction(1,576),Fraction(-1,9),-385526315),
        ('B=31A/5',Fraction(31,2880),Fraction(1,5),-3620040124),
        ('B=4A/5',Fraction(1,720),Fraction(4,5),-1538198026),
        ('complex probe',0,Fraction(9072,11),-1141009435895),
        ('integer probe 1',896,92160,-405807237541622),
        ('integer probe 2',2256,72704,-801315534069517),
        ('real difference',0,-432,597671609277),
    )
    for name,pk,pj,floor in combinations:
        coeff=[pk*k+pj*j for k,j in total]
        lo=sum(min(c*x,c*y) for c,(x,y) in zip(coeff,zetas))
        hi=sum(max(c*x,c*y) for c,(x,y) in zip(coeff,zetas))
        assert floor<lo*10**20<=hi*10**20<floor+1, name
        print(name,'between',floor,'and',floor+1,'times 10^-20',flush=True)


def full_controls(q,channel):
    if q==12:
        rows = [(*powers,value) for powers,value in zip(PROBES,FULL_PROBES[channel])]
    else:
        bank = json.loads((HERE.parent/'nonSP/verify/prism_forcer/periods.json').read_text())
        rows = [(*[k+1 for k in ex], Fraction(c[channel])*factorial(5)/prod(map(factorial,ex)))
                for ex,c in zip(bank['exponents'],bank['coefficients'])]
    return BytesIO(('\n'.join(' '.join(map(str,row)) for row in rows)).encode())


def build_period_tables(native,q,tmp,reselect=False):
    if q not in (6,12):
        raise ValueError('The complete prism table requires q=6 or q=12.')
    start=[native,'--all-periods',str(q)]
    pivot,sparse,final,prefix=[tmp/(n+str(q)) for n in ('direct','sparse','final','table')]
    run([*start,pivot],check=True)
    with gzip_open(HERE/'quartic_prism_rank_rows.gz','rb') as cert:
        for dim in (6,12):
            for name,path in (('sparse',sparse),('final',final)):
                qtxt,section,size=cert.readline().decode().split()
                assert (int(qtxt),section)==(dim,name)
                data=cert.read(int(size))
                assert len(data)==int(size)
                if q==dim:
                    path.write_bytes(data)
        assert cert.read()==b''
    solve=[*start,pivot,'--full-solve',sparse,final,prefix]
    # Tied direct rows can change with the C++ hash-table implementation.
    if reselect or run(solve).returncode:
        run([*start,pivot,'--sparse',sparse,str(2000 if q==6 else 40000)],check=True)
        run([*start,pivot,'--final-rank',sparse,final,str(100 if q==6 else 3500)],check=True)
        run(solve,check=True)
    return prefix


def check_full(reselect=False):
    check_signs()
    physical_response_bounds()
    with TemporaryDirectory() as scratch:
        tmp = Path(scratch)
        native = tmp/'ibp'
        compile_program('quartic_prism_physical_ibp.cpp',native)
        for q in (6,12):
            start = [native,'--all-periods',str(q)]
            prefix = build_period_tables(native,q,tmp,reselect)
            run([*start,'--check-integers',prefix],check=True)
            reduced = tmp/'grouped'
            if q==12:
                with gzip_open(HERE/'quartic_prism.weights.gz','rb') as bank:
                    out = stream([*start,'--group',reduced],bank)
                assert 'physical rows 1598657 reflected groups 241543' in out
            for ch in range(3):
                table = str(prefix)+str(ch)
                count = 1481 if q==6 else 2
                out = stream([*start,'--compare',table],full_controls(q,ch))
                assert f'independent cuts {count}' in out
                print(out.strip(),flush=True)
                if q==12:
                    with gzip_open(HERE/'quartic_prism.weights.gz','rb') as bank:
                        sums(stream([*start,'--sum',table],bank),FULL_TOTAL[ch],1598657)
                    with reduced.open('rb') as groups:
                        out=stream([*start,'--group-sum',table],groups)
                    line=next(l for l in out.splitlines() if l.startswith('grouped sum '))
                    _,_,n,k,j=line.split()
                    assert int(n)==241543 and (Fraction(k),Fraction(j))==FULL_TOTAL[ch]


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--all-periods',action='store_true',help='solve and verify all three physical columns')
    parser.add_argument('--reselect',action='store_true',help='regenerate the full sparse and dense pivot rows')
    parser.add_argument('--form',type=Path)
    parser.add_argument('--forcer',type=Path)
    args = parser.parse_args()
    if (args.form is None) != (args.forcer is None):
        parser.error('Use --form and --forcer together')
    if args.reselect and not args.all_periods:
        parser.error('--reselect requires --all-periods')
    if args.all_periods:
        check_full(args.reselect)
        if args.form is not None:
            check_forcer(args.form,args.forcer)
    else:
        check(args.form,args.forcer)
