"""Exact two-unit angular transitions on the residual mixed IIB bipartite graphs."""
from argparse import ArgumentParser
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from fractions import Fraction as F
from itertools import product,zip_longest
from math import factorial,lcm,prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import gzip
import json

from flint import fmpq
from mixed_bipartite_linear_bank import HERE,charts,completed,pack,primes_for,prism_charts,unpack,weights
from quartic_bipartite_wheel import geometry,jets_on_pole,rational_sum
from quartic_bipartite_sums import connection
from quartic_prism_physical_ibp import compile_program,zeta_interval


def height(q):
    lam=q-1
    last=2*lam+1
    mass=0
    poles=set()
    for b in range(1,q):
        s=0
        for k in range(q-b):
            coeff,roots=connection(q,b,k)
            poles.update(roots)
            # Simple partial fractions of each connection before multiplying the four rim edges.
            pieces=[coeff/prod(y-x for y in roots if y!=x) for x in roots] if roots else [coeff]
            assert all(v.denominator==1 for v in pieces)
            s+=sum(abs(int(v)) for v in pieces)
        mass=max(mass,s)
    # Relative angular factors. A positive power in this profile means a denominator.
    base=Counter({a:-1 for a in range(1,2*lam)})
    base[lam]+=3
    sectors=[(0,(),()),(0,(lam,lam+1),(1,2*lam)),(0,(lam,lam-1),(0,2*lam-1)),
             (1,(lam+1,lam+1),(lam,2*lam)),
             (2,(lam+1,lam+2,lam+2),(lam,2*lam,2*lam+1))]
    for jump,den,num in sectors:
        angular=base.copy()
        angular.update(den)
        angular.subtract(num)
        for a in range(1,last+1):
            connection_order=2*int(a in poles)+2*int(a-jump in poles)
            assert angular[a]+connection_order+3<=7
    # At most 13 inverse linear factors: four connections, three base angular,
    # three relative angular and three radial. The half-unit pole circles
    # bound every coefficient of all the sector partial fractions.
    lattice=lcm(*range(1,last+1))**13
    denominator=factorial(lam)**3*factorial(2*lam-1)*lattice
    limit=lattice*lam**3*20*16*mass**4*2**13*(last+2)**(2*lam+2)*7*last**2
    return denominator,limit


def certified(q,keys,program,directory,*,degree=2,proof=None,fixed_spin=False,prime_workers=1):
    ordered=sorted(keys)
    manifest=directory/'keys'
    manifest.write_text(''.join(f'{k}\n' for k in ordered))
    if proof is None:
        denominator,limit=height(q)
        order,last,mode=7,2*q-1,'--quadratic'
    else:
        order,last,denominator,limit=proof
        mode=f'--fixed-spin={degree}' if fixed_spin else f'--degree={degree}'
    primes=list(primes_for(limit))
    modulus=prod(primes)
    lifts=[(modulus//p)*((denominator*pow(modulus//p,-1,p))%p) for p in primes]
    def evaluate_field(p):
        kernel=directory/f'kernel-{p}'
        kernel.write_text(''.join(' '.join(map(str,row))+'\n' for row in weights(q,p,order,last)))
        dest=directory/f'field-{p}'
        run([program,str(q),str(p),kernel,manifest,dest,mode],check=True,capture_output=True)
        return dest
    if prime_workers==1:
        paths=list(map(evaluate_field,primes))
    else:
        with ThreadPoolExecutor(max_workers=prime_workers) as workers:
            paths=list(workers.map(evaluate_field,primes))
    values={}
    with ExitStack() as stack:
        streams=[stack.enter_context(path.open()) for path in paths]
        for key,*rows in zip_longest(ordered,*streams):
            assert key is not None and all(r is not None for r in rows)
            coeffs=[]
            for column in zip(*(r.split() for r in rows),strict=True):
                v=sum(int(r)*c for r,c in zip(column,lifts,strict=True))%modulus
                if 2*v>modulus:
                    v-=modulus
                assert abs(v)<=limit
                coeffs.append(v)
            assert len(coeffs)==order+1 and all(not v for j,v in enumerate(coeffs) if j not in (0,3,5))
            values[key]=tuple(F(coeffs[j],denominator) for j in (0,3,5))
    print('dimension, charts, angular evaluations and primes:',2*q,len(keys),order*last,len(primes),flush=True)
    return values


def reflected(program,q,inputs,kind,formal=False):
    data=''.join(' '.join(map(str,row))+'\n' for row in inputs)
    if formal and not isinstance(kind,int):
        raise ValueError('A formal chart census requires an explicit diagonal degree')
    label='--formal-numerator-wheel-charts' if formal else '--numerator-wheel-charts'
    flags=[label,str(kind)] if isinstance(kind,int) else ['--'+kind+'-charts']
    result=run([program,str(q),*flags],input=data,text=True,capture_output=True,check=True)
    return {tuple(map(int,r[:9])):(F(r[9]),int(r[11]),F(r[12]))
            for line in result.stdout.splitlines() if (r:=line.split())}


def control(complete,wheel,scratch):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    known={}
    forward={}
    inputs=[]
    for ex,columns in zip(ref['exponents'],ref['coefficients'],strict=True):
        seed=tuple(a+1 for a in ex)
        inputs.append((*seed,1,0))
        graph=completed(seed,6)
        expected=tuple(F(str(v))*factorial(5)/prod(map(factorial,ex)) for v in columns)
        key=min(prism_charts(graph,6))
        assert forward.setdefault(key,expected)==expected
        for chart in charts(graph,6,2):
            assert known.setdefault(pack(chart),expected)==expected
    direct=len(known)
    selected=reflected(complete,6,inputs,'quadratic-wheel')
    for seed,(_,k,scale) in selected.items():
        expected=tuple(x/scale for x in forward[seed])
        assert known.setdefault(k,expected)==expected
    assert direct==230 and len(selected)==113 and len(known)==295
    print('Independent direct charts and reflected seeds:',direct,len(selected),flush=True)
    assert certified(6,known,wheel,scratch)==known


def direct_period(spokes,ring,q):
    """Sum the individual rational functions through their Laurent coefficients."""
    lam=q-1
    contour,cuts=geometry()
    scalar=(0,0,0,0)
    sectors=[(scalar,(4,0,0,0),1,(),()),(scalar,(0,0,4,0),1,(),()),
             (scalar,(2,0,2,0),2,(),()),
             (scalar,(2,0,2,0),1,(lam,lam+1),(1,2*lam)),
             (scalar,(2,0,2,0),1,(lam,lam-1),(0,2*lam-1))]
    for high in ((0,1,1,0),(1,0,0,1)):
        for radius in ((1,0,3,0),(3,0,1,0)):
            sectors.append((high,radius,-2,(lam+1,lam+1),(lam,2*lam)))
        sectors.append((tuple(2*x for x in high),(2,0,2,0),1,
                        (lam+1,lam+2,lam+2),(lam,2*lam,2*lam+1)))
    total=defaultdict(fmpq)
    for high,radius,sign,den,num in sectors:
        angular=Counter({a:-1 for a in range(1,2*lam)})
        angular[lam]+=3
        angular.update(den)
        angular.subtract(num)
        radial={}
        for key,(inside,crossing,internal) in cuts.items():
            shift=sum(radius[i] for i in inside)+sum(high[i] for i in crossing)
            assert shift%2==0
            radial[key]=sum(q-spokes[i] for i in inside)-sum(ring[i] for i in internal)+shift//2
        assert min(radial.values())>0
        for ks in product(*(range(q-b) for b in ring)):
            profile=angular.copy()
            scale=fmpq(sign)
            for b,k,delta in zip(ring,ks,high,strict=True):
                coeff,roots=connection(q,b,k)
                scale*=fmpq(str(coeff))
                profile.update(a+delta for a in roots)
            rad={key:x+sum(ks[i] for i in cuts[key][1]) for key,x in radial.items()}
            for poles,n in Counter(tuple(sorted(rad[k] for k in row)) for row in contour).items():
                whole=profile.copy()
                whole.update(poles)
                partial=rational_sum(tuple(sorted((a,p) for a,p in whole.items() if p)))
                for k,v in partial.items():
                    total[k]+=n*scale*v
            jets_on_pole.cache_clear()
    assert all(v==0 or k in (0,3,5) for k,v in total.items())
    norm=fmpq(lam**3,factorial(lam)**3*factorial(2*lam-1))
    return tuple(F(str(norm*total[k])) for k in (0,3,5))


EXPECTED=(F(2946614573584234833492017791125139,58070236078413279589220352000000000),
          F(20683963022424458481150037,152405717430957838848000000),
          -F(810035414794129859,3927172681688256000))


def physical(complete,wheel,scratch,check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as f:
        data={tuple(map(int,r[:9])):r[:11] for line in f if (r:=line.split())[11]=='0'}
    assert len(data)==35563
    for kind,count in (('wheel',10940),('linear-wheel',10583)):
        removed=reflected(complete,13,data.values(),kind)
        assert len(removed)==count and removed.keys()<=data.keys()
        for seed in removed:
            del data[seed]
    assert len(data)==14040
    selected=reflected(complete,13,data.values(),'quadratic-wheel')
    assert len(selected)==5903 and selected.keys()<=data.keys()
    values=certified(13,{k for _,k,_ in selected.values()},wheel,scratch)
    if check:
        short=sorted(values,key=lambda k:(prod(13-b for b in unpack(k)[1]),k))[:6]
        for key in short:
            assert direct_period(*unpack(key),13)==values[key]
        print('Direct q=13 Laurent checks:',len(short),flush=True)
    total=tuple(sum(w*scale*values[k][i] for w,k,scale in selected.values()) for i in range(3))
    assert total==EXPECTED
    lo=hi=total[0]
    for v,n in zip(total[1:],(3,5)):
        a,b=sorted(v*z for z in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    assert -6146<10**20*lo<=10**20*hi<-6145
    print('5903 new components:',*total,flush=True)
    print('10^20 times the subtotal strictly between',(10**20*lo)//1,-((-10**20*hi)//1),flush=True)
    return total


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Test three exact period columns against independent propagators')
    parser.add_argument('--controls-only',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory() as work:
        directory=Path(work)
        complete,wheel=directory/'complete',directory/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_wheel_bank.cpp',wheel)
        if args.check or args.controls_only:
            control(complete,wheel,directory)
        if not args.controls_only:
            physical(complete,wheel,directory,args.check)
