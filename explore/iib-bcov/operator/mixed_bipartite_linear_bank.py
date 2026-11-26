"""Certify neighboring-spin wheel cuts of the mixed four-loop IIB Ward chamber."""
from argparse import ArgumentParser
from fractions import Fraction as F
import gzip
from itertools import zip_longest
import json
from math import comb, factorial, isqrt, lcm, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from mixed_bipartite_linear_wheel import HERE, charts, period
from quartic_prism_completion import completed, prism_charts
from quartic_prism_physical_ibp import compile_program, zeta_interval


def pack(wheel):
    n = 0
    for a in sum(wheel, ()):
        n = 16*n+a
    return n


def unpack(key):
    row = tuple((key>>(4*i)) & 15 for i in reversed(range(8)))
    return row[:4],row[4:]


def weights(q,p,order=6,last=None):
    """Partial fractions of rational cardinal functions on all admissible angular poles."""
    if last is None:
        last=2*(q-1)
    n=order*last
    harmonic = {(a,h):sum(pow(i,-h,p) for i in range(1,a))%p
                for a in range(1,last+1) for h in range(1,order+1)}
    regular = []
    for a in range(1,last+1):
        jet=[1]+[0]*(order-1)
        for b in range(1,last+1):
            if b==a:
                continue
            series=[(-1)**i*comb(order+i-1,i)*pow(b-a,-order-i,p)%p for i in range(order)]
            jet=[sum(jet[k-i]*series[i] for i in range(k+1))%p for k in range(order)]
        for j in range(n):
            jet=[((-j-a)*x+(jet[k-1] if k else 0))%p for k,x in enumerate(jet)]
        regular.append(jet)
    for j in range(n):
        scale=prod(pow(j+b,order,p) for b in range(1,last+1))%p
        scale=scale*pow((-1)**(n-1-j)*factorial(j)*factorial(n-1-j),-1,p)%p
        row=[0]*(order+1)
        for a,jet in enumerate(regular,1):
            inv=pow(-a-j,-1,p)
            leading=0
            for i,x in enumerate(jet):
                leading=(x-leading)*inv%p
                level=order-i
                c=scale*leading%p
                row[level]=(row[level]+c)%p
                row[0]=(row[0]-c*harmonic[a,level])%p
        yield row


def bounds(q):
    lam=q-1
    common=lcm(*range(1,2*lam))**9
    mass=1
    for b in range(1,lam):
        dots=lam-b
        mass=max(mass,comb(lam-1,b-1)*sum(abs(dots-k-j)*comb(dots,k)*comb(dots,j)
                                          for k in range(dots+1) for j in range(dots+1)))
    # Up to nine linear denominator factors survive in each partial fraction
    # of the four connections, angular prefactor and three radial integrals.
    denominator=factorial(lam)**3*factorial(2*lam-1)*common
    numerator=common*lam**3*80*mass**4*(2*lam+1)**(2*lam)*2**9*9*(2*lam+1)**2
    return denominator,numerator


def primes_for(limit):
    p=1000000007
    modulus=1
    while modulus<=2*limit:
        if all(p%i for i in range(3,isqrt(p)+1,2)):
            yield p
            modulus*=p
        p-=2


def certified(q,keys,exe,tmp):
    ordered=sorted(keys)
    keys_file=tmp/'keys'
    keys_file.write_text(''.join(f'{k}\n' for k in ordered))
    denominator,bound=bounds(q)
    moduli=list(primes_for(bound))
    modulus=prod(moduli)
    crt=[(modulus//p)*((denominator*pow(modulus//p,-1,p))%p) for p in moduli]
    files=[]
    for p in moduli:
        kernel,output=tmp/'kernel',tmp/f'field-{p}'
        kernel.write_text(''.join(' '.join(map(str,r))+'\n' for r in weights(q,p)))
        run([exe,str(q),str(p),kernel,keys_file,output,'--signed'],check=True,capture_output=True)
        files.append(output.read_text().splitlines())
    result={}
    for key,*rows in zip_longest(ordered,*files):
        scaled=[]
        for residues in zip(*(r.split() for r in rows),strict=True):
            v=sum(int(x)*a for x,a in zip(residues,crt,strict=True))%modulus
            if 2*v>modulus:
                v-=modulus
            assert abs(v)<=bound
            scaled.append(v)
        assert len(scaled)==7 and all(scaled[i]==0 for i in (1,2,4,6))
        result[key]=tuple(F(scaled[i],denominator) for i in (0,3,5))
    print(f'q={q}:',len(ordered),'integrals,',12*(q-1),'nodes,',len(moduli),'primes',flush=True)
    return result


def control(complete,exe,tmp):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    known={}
    forward={}
    inputs=[]
    for ex,columns in zip(ref['exponents'],ref['coefficients'],strict=True):
        seed=tuple(a+1 for a in ex)
        expected=tuple(F(str(x))*factorial(5)/prod(map(factorial,ex)) for x in columns)
        graph=completed(seed,6)
        key=min(prism_charts(graph,6))
        assert forward.setdefault(key,expected)==expected
        inputs.append(' '.join(map(str,(*seed,1,0))))
        for wheel in charts(graph,6):
            assert known.setdefault(pack(wheel),expected)==expected
    assert len(known)==613
    reflected=run([complete,'6','--linear-wheel-charts'],input='\n'.join(inputs),text=True,
                  capture_output=True,check=True)
    measurements=0
    for line in reflected.stdout.splitlines():
        row=line.split()
        key=int(row[11])
        expected=tuple(c/F(row[12]) for c in forward[tuple(map(int,row[:9]))])
        assert known.setdefault(key,expected)==expected
        measurements+=1
    assert certified(6,known,exe,tmp)==known
    print('q=6 direct charts and independent reflected seed measurements:',613,measurements,flush=True)


def physical(complete,exe,tmp,check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as f:
        data=[r for l in f if (r:=l.split())[11]=='0']
    inputs=''.join(' '.join(r[:11])+'\n' for r in data)
    wheel=run([complete,'13','--wheel-charts'],input=inputs,text=True,capture_output=True,check=True)
    scalar={tuple(map(int,l.split()[:9])) for l in wheel.stdout.splitlines()}
    direct={}
    for row in data:
        seed=tuple(map(int,row[:9]))
        if seed in scalar:
            continue
        candidates=list(charts(completed(seed,13),13))
        if candidates:
            direct[seed]=pack(min(candidates))
    inputs=''.join(' '.join(r[:11])+'\n' for r in data if tuple(map(int,r[:9])) not in scalar)
    reflected=run([complete,'13','--linear-wheel-charts'],input=inputs,text=True,capture_output=True,check=True)
    selected=[]
    for line in reflected.stdout.splitlines():
        row=line.split()
        seed=tuple(map(int,row[:9]))
        selected.append((seed,F(row[9]),int(row[11]),F(row[12])))
    keys={key for _,_,key,_ in selected}|set(direct.values())
    lookup=certified(13,keys,exe,tmp)
    if check:
        probes=sorted(lookup,key=lambda k:(prod(13-x for x in unpack(k)[1]),k))[:12]
        for key in probes:
            assert period(*unpack(key),13)==lookup[key]
    total=[F(),F(),F()]
    direct_sum=[F(),F(),F()]
    for seed,weight,key,scale in selected:
        contribution=tuple(scale*x for x in lookup[key])
        if seed in direct:
            assert contribution==lookup[direct[seed]]
            for i,c in enumerate(contribution):
                direct_sum[i]+=weight*c
        for i,c in enumerate(contribution):
            total[i]+=weight*c
    assert len(direct)==4484 and direct.keys()<={s for s,_,_,_ in selected}
    assert tuple(direct_sum)==(-F(148622573162634799711770478892069,19356745359471093196406784000000000),
                               -F(12519193368531968593421,609622869723831355392000),F(3404741402515643,109088130046896000))
    assert len(selected)==10583 and tuple(total)==(
        -F(21349904253811884693303442588528537,116140472156826559178440704000000000),
        -F(269760817092032201695082203,548660582751448219852800000),
         F(2934584497151629337,3927172681688256000))
    lo=hi=total[0]
    for c,h in zip(total[1:],(3,5)):
        a,b=sorted(c*x for x in zeta_interval(h,terms=100))
        lo+=a
        hi+=b
    print('components:',len(selected),'q=13 signed wheel columns:',*total,flush=True)
    assert -47677<10**20*lo<10**20*hi<-47676
    print('10^20 enclosure:',(10**20*lo)//1,-((-10**20*hi)//1))


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory() as work:
        tmp=Path(work)
        complete,exe=tmp/'complete',tmp/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_wheel_bank.cpp',exe)
        if args.check:
            control(complete,exe,tmp)
        physical(complete,exe,tmp,args.check)
