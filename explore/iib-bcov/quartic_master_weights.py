"""Assemble the two genus-four BCOV numerators on the primitive nine-edge cores."""
from argparse import ArgumentParser
from collections import defaultdict
from gzip import GzipFile, open as open_gzip
from hashlib import file_digest
from itertools import islice
import json
from math import comb, factorial, lcm, prod
from pathlib import Path
from shutil import copyfileobj
from subprocess import PIPE, Popen, run
import sys
from tempfile import TemporaryDirectory

from flint import fmpq
import numpy as np

from quartic_planar_cores import cubic_cores, graph_shapes
from quartic_planar_periods import core_data, polynomials, reduction
from four_loop import AUTOMORPHISMS, K33, PERMS, PRISM, two_cuts


def orbit_groups(poly,perms,keep,total=48):
    columns = defaultdict(fmpq)
    shifts = np.arange(32,-1,-4,dtype=np.uint64)
    ordering = np.asarray(keep)
    # python-flint's terms() first creates lists for every monomial and
    # coefficient; positional access keeps these arrays bounded by the batch.
    terms = ((poly.monomial(i),poly.coefficient(i)) for i in range(len(poly)))
    while batch:=list(islice(terms,8192)):
        powers = np.array([ex for ex,v in batch],dtype=np.uint64)[:,ordering]+1
        assert powers.min()>=1 and powers.max()<12 and np.all(powers.sum(axis=1)==total)
        orbit = np.full(len(batch),2**36,dtype=np.uint64)
        for permutation in perms:
            candidate = np.sum(powers[:,permutation]<<shifts,axis=1,dtype=np.uint64)
            np.minimum(orbit,candidate,out=orbit)
        for key,(_,value) in zip(orbit,batch):
            columns[int(key)] += value
    return {k:v for k,v in columns.items() if v}


def unpack(key):
    return tuple((key>>i)&15 for i in range(32,-1,-4))


def assemble(name,output):
    core = cubic_cores()[3 if name=='prism' else 4]
    shapes = list(graph_shapes(core))
    assert len(shapes)==(28 if name=='prism' else 11)
    moves,kind,keep = reduction(core)
    assert kind==name and not moves
    permutations = np.asarray(PERMS if name=='prism' else AUTOMORPHISMS,dtype=np.intp)
    base = defaultdict(fmpq)
    retained = defaultdict(fmpq)
    for i,edges in enumerate(shapes):
        reduced,polynomial,rest = polynomials(edges,components=True)
        assert reduced==core
        for key,value in orbit_groups(polynomial,permutations,keep,43).items():
            base[key] += value
        for key,value in orbit_groups(rest,permutations,keep).items():
            retained[key] += value
        print(name,i+1,'shapes',len(base),len(retained),flush=True)
    ctx,u,rho,pair = core_data(core)
    # The two resistance factors are automorphism invariant. Collect all
    # four-slot placements before multiplying. Restore the original edge
    # order on a representative; its orbit average has the same integral.
    inverse = [keep.index(i) for i in range(len(core))]
    representatives = {tuple(unpack(k)[i]-1 for i in inverse):v for k,v in base.items() if v}
    polynomial = ctx.from_dict(representatives)
    k = orbit_groups(polynomial*rho/8,permutations,keep)
    j = orbit_groups(polynomial*pair,permutations,keep)
    for key,value in k.items():
        retained[key] += value
    columns = [{k:24*factorial(4)**2*v for k,v in retained.items() if v},
               {k:48*v for k,v in j.items() if v}]
    with output.open('wb') as raw, GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as handle:
        for key in sorted(columns[0].keys()|columns[1].keys()):
            powers = unpack(key)
            norm = fmpq(prod(factorial(a-1) for a in powers),factorial(11))
            row = [*powers,*(norm*col.get(key,0) for col in columns)]
            handle.write((' '.join(map(str,row))+'\n').encode('ascii'))
    print(name,'scalar moment orbits',len(columns[0]),len(columns[1]),flush=True)


def control(name,form,forcer):
    import sympy as sp
    from fractions import Fraction as F
    prism = name=='prism'
    core = cubic_cores()[3 if prism else 4]
    reduced,terms = polynomials(next(graph_shapes(core)))
    moves,kind,keep = reduction(reduced)
    assert not moves and kind==name
    perms = np.asarray(PERMS if prism else AUTOMORPHISMS,dtype=np.intp)
    row = orbit_groups(terms[1][0],perms,keep)
    powers = (4,4,5,9,6,4,5,5,6) if prism else (4,4,5,4,5,6,5,9,6)
    key = sum(v<<(32-4*j) for j,v in enumerate(powers))
    assert row[key]==(fmpq(-1,14) if prism else fmpq(9,700))
    with TemporaryDirectory(prefix='bcov-genus-four-control-') as tmp:
        folder = Path(tmp)
        for file in ('forcer','forcer.h'):
            (folder/file).symlink_to((forcer/file).resolve(strict=True))
        norm = F(prod(factorial(a-1) for a in powers),factorial(11))
        result = two_cuts(folder,form.resolve(strict=True),name,12,{powers:norm})
    print(name,powers,result,flush=True)
    z3,z5 = sp.symbols('z3 z5')
    if prism:
        assert result == (-sp.Rational(531372573405073,2048959660011264000000)
                          -sp.Rational(80588077,116154175737600)*z3+sp.Rational(15,14226212)*z5)
    else:
        assert result == (-sp.Rational(186964375676833,9579551657195520000)
                          -sp.Rational(53055501437,1016349037704000)*z3
                          +sp.Rational(1693,21339318)*z5)


def audit(name,weights):
    core = cubic_cores()[3 if name=='prism' else 4]
    samples = [[1]*9]+[[1+int(i==j) for i in range(9)] for j in range(9)]
    predicted = [[fmpq(0),fmpq(0)] for _ in range(2)]
    for edges in graph_shapes(core):
        point = polynomials(edges,samples=samples)
        for i in range(2):
            predicted[0][i] += point[0][i]
            predicted[1][i] += sum(row[i] for row in point[1:])/9
    sums = [[fmpq(0),fmpq(0)] for _ in range(2)]
    with open_gzip(weights,'rt') as file:
        for line in file:
            *powers,k,j = line.split()
            exponent = [int(a)-1 for a in powers]
            undo_cut = fmpq(factorial(11),prod(factorial(a) for a in exponent))
            doubled = fmpq(sum(2**a for a in exponent),9)
            for i,entry in enumerate((k,j)):
                weight = fmpq(entry)*undo_cut
                sums[0][i] += weight
                sums[1][i] += weight*doubled
    assert predicted==sums
    print(name,'direct nine-unit-length numerators',*sums[0],flush=True)
    print(name,'direct mean single-doubled-length numerators',*sums[1],flush=True)


def hepp_zero(name,weights,binary):
    # The native integer recursion divides out the cut normalization. The
    # output gives the signed and absolute tropical moments in each channel.
    core = PRISM if name=='prism' else K33
    header = (' '.join(str(i) for edge in core for i in edge)+'\n').encode('ascii')
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'nonSP/verify/hepp_series'))
    from period import HeppSeries
    powers = (4,4,5,9,6,4,5,5,6) if name=='prism' else (4,4,5,4,5,6,5,9,6)
    ex = tuple(a-1 for a in powers)
    series = HeppSeries(core,12,{ex:1})
    expected = next(series.moments(0))[1]
    norm = fmpq(prod(factorial(k) for k in ex),factorial(11))
    sample = (' '.join(map(str,(*powers,-norm,2*norm)))+'\n').encode('ascii')
    checked = run([binary.resolve(strict=True)],input=header+sample,capture_output=True,check=True)
    assert list(map(fmpq,checked.stdout.decode('ascii').split()))==[1,-expected,expected,2*expected,2*expected]
    with Popen([binary.resolve(strict=True)],stdin=PIPE,stdout=PIPE) as child:
        child.stdin.write(header)
        with open_gzip(weights,'rb') as lines:
            copyfileobj(lines,child.stdin)
        child.stdin.close()
        result=child.stdout.read().decode('ascii').split()
        if child.wait() or len(result)!=5:
            raise RuntimeError('The tropical moment computation failed.')
    print(name,'distinct moment rows',result[0],flush=True)
    center,radius,normalization = series.expansion('gegenbauer')
    for label,signed,bound in zip(('K4','J4'),result[1::2],result[2::2]):
        print(label,'raw moment',signed,'absolute',bound,'zeroth Gegenbauer term',
              fmpq(signed)*normalization,'tail bound',fmpq(bound)*(1-normalization),flush=True)


def gegenbauer(name,weights,order):
    """Evaluate exact finite Gegenbauer partial sums and their absolute remainder bounds."""
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'nonSP/verify/hepp_series'))
    from period import HeppSeries
    if order<0:
        raise ValueError('The series order must be nonnegative.')
    core = PRISM if name=='prism' else K33
    partial = [fmpq(0),fmpq(0)]
    bounds = [fmpq(0),fmpq(0)]
    count = 0
    with open_gzip(weights,'rt') as file:
        while lines:=list(islice(file,64)):
            columns = [{},{}]
            for line in lines:
                *indices,k,j = line.split()
                ex = tuple(int(a)-1 for a in indices)
                norm = fmpq(prod(factorial(a) for a in ex),factorial(11))
                for col,value in zip(columns,(k,j)):
                    col[ex] = fmpq(value)/norm
            count += len(lines)
            for i,col in enumerate(columns):
                series = HeppSeries(core,12,col)
                for power,term,value,error in series.terms(order):
                    if power==order:
                        partial[i] += value
                        bounds[i] += error
    print(name,'distinct moment rows',count,flush=True)
    for label,p,b in zip(('K4','J4'),partial,bounds):
        print(label,'order',order,'partial sum',p,'absolute tail bound',b,flush=True)
    return partial,bounds


def higher_moment(name,weights,binary,order):
    """Expand the Kirchhoff polynomial orbitwise before evaluating an exact raw moment."""
    if not 1<=order<=3:
        raise ValueError('The native higher-moment evaluator supports orders one through three.')
    core = PRISM if name=='prism' else K33
    header = (' '.join(str(i) for edge in core for i in edge)+'\n').encode('ascii')
    den = 1
    with open_gzip(weights,'rt') as file:
        for line in file:
            *p,k,j = line.split()
            undo = fmpq(factorial(11),prod(factorial(int(a)-1) for a in p))
            den = lcm(den,int((undo*fmpq(k)).denominator),int((undo*fmpq(j)).denominator))
    with Popen([binary.resolve(strict=True),str(order),str(den)],stdin=PIPE,stdout=PIPE) as child:
        child.stdin.write(header)
        with open_gzip(weights,'rb') as file:
            copyfileobj(file,child.stdin)
        child.stdin.close()
        result=[]
        for degree,line in enumerate(child.stdout,1):
            row=line.decode('ascii').split()
            if len(row)!=4:
                raise ValueError('Unexpected native moment row.')
            print(name,'order',degree,'input and expanded rows',*row[:2],flush=True)
            print('K4 J4 raw moments',*row[2:],flush=True)
            result.append(row)
        if child.wait() or len(result)!=order:
            raise RuntimeError('The higher tropical moment computation failed.')
    return tuple(map(fmpq,result[-1][2:]))


def sign_bounds(name,certificate):
    """Enclose the integral against the positive measure of the second BCOV period."""
    from fractions import Fraction as F
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'nonSP/verify/hepp_series'))
    from moment_bounds import integer_enclosure
    data = json.loads(certificate.read_text())[name]
    with certificate.with_name(data['file']).open('rb') as weights:
        assert file_digest(weights,'sha256').hexdigest()==data['sha256']
    moments = [-F(row[1]) for row in data['raw_moments_K_J']]
    endpoints = integer_enclosure(moments,12,75 if name=='prism' else 81)
    coarse = integer_enclosure(moments[:2],12,75 if name=='prism' else 81)
    assert coarse[0]<endpoints[0]<endpoints[1]<coarse[1]
    lo,hi = endpoints
    scale = 10**20
    compact = [int(lo*scale),-(-hi*scale//1)]
    assert compact==data['outward_minus_J_times_1e20']
    print(name,'strict bounds on -10^20 J',*compact,flush=True)
    return endpoints


def physical_bounds(certificate):
    from fractions import Fraction as F
    lower = upper = F(0)
    for name in ('prism','K33'):
        lo,hi=sign_bounds(name,certificate)
        lower += 1728*lo
        upper += 1728*hi
    terms = 80
    zeta = sum(F(4,3*2**(k+1))*sum(F((-1)**j*comb(k,j),(j+1)**3) for j in range(k+1))
               for k in range(terms))
    a,b = F(507937,64744680000),F(173064182111,18355116780000000)
    lower += a*zeta-b
    upper += a*(zeta+F(4,3*2**terms))-b
    scale = 10**18
    compact = [int(scale*lower),-(-scale*upper//1)]
    assert compact==json.loads(certificate.read_text())['full_minus_1728_J_times_1e18']
    print('strict bounds on -10^18 * 1728 J4',*compact,flush=True)
    return lower,upper


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('core',choices=['prism','K33','both'])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('output',type=Path,nargs='?')
    group.add_argument('--control',nargs=2,type=Path,metavar=('FORM','FORCER'))
    group.add_argument('--audit',type=Path,metavar='WEIGHTS')
    group.add_argument('--hepp-zero',nargs=2,type=Path,metavar=('WEIGHTS','NATIVE'))
    group.add_argument('--gegenbauer',nargs=2,metavar=('WEIGHTS','ORDER'))
    group.add_argument('--native-moment',nargs=3,metavar=('WEIGHTS','NATIVE','ORDER'))
    group.add_argument('--bounds',type=Path,metavar='CERTIFICATE')
    args = parser.parse_args()
    if args.core=='both':
        if not args.bounds:
            parser.error('The combined projection uses --bounds CERTIFICATE.')
        physical_bounds(args.bounds)
    elif args.control:
        control(args.core,*args.control)
    elif args.audit:
        audit(args.core,args.audit)
    elif args.hepp_zero:
        hepp_zero(args.core,*args.hepp_zero)
    elif args.gegenbauer:
        gegenbauer(args.core,Path(args.gegenbauer[0]),int(args.gegenbauer[1]))
    elif args.native_moment:
        weights,binary,order=args.native_moment
        higher_moment(args.core,Path(weights),Path(binary),int(order))
    elif args.bounds:
        sign_bounds(args.core,args.bounds)
    else:
        assemble(args.core,args.output)
