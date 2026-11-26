"""All remaining convergent polynomial-diagonal cuts in one mixed bipartite chamber."""
from argparse import ArgumentParser
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations,product
from math import comb,factorial,lcm,prod
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import json

from flint import fmpq
from mixed_bipartite_quadratic_bank import (
    HERE,certified,charts,compile_program,completed,connection,geometry,jets_on_pole,pack,
    prism_charts,rational_sum,reflected,unpack,zeta_interval,
)


@cache
def paths(q,degree):
    """Rational angular multipliers for the monomials of (r0²+r2²-2r0r2 x)^degree."""
    lam=q-1
    assert 0<=degree<=lam
    out=defaultdict(Counter)
    for p in range(degree+1):
        for steps in product((-1,1),repeat=p):
            j=sum(steps)
            if j<0:
                continue
            # Start with the scalar angular norm and the ratio for the final degree n+j.
            angular=Counter({a:-1 for a in range(1,2*lam)})
            angular[lam]+=3
            if j:
                angular[lam+j]+=2
                angular[lam]-=2
                angular.update(range(1,j+1))
                angular.subtract(range(2*lam,2*lam+j))
            spin=0
            for step in steps:
                angular[lam+spin]+=1
                angular[spin+(1 if step==1 else 2*lam-1)]-=1
                spin+=step
            assert spin==j
            profile=tuple(sorted((a,v) for a,v in angular.items() if v))
            for i in range(degree-p+1):
                r=2*i+p
                sign=(-1)**p*comb(degree,p)*comb(degree-p,i)
                out[r,j][profile]+=sign
    return {sector:{f:v for f,v in rows.items() if v} for sector,rows in out.items()}


def height(q,d,displacement=0):
    lam=q-1
    last=2*lam+d-1
    poles=set()
    mass=0
    for b in range(1,q):
        s=0
        for k in range(q-b):
            c,roots=connection(q,b,k)
            poles.update(roots)
            pieces=[c/prod(y-x for y in roots if y!=x) for x in roots] if roots else [c]
            assert all(x.denominator==1 for x in pieces)
            s+=sum(abs(int(x)) for x in pieces)
        mass=max(mass,s)
    order=inverse_factors=numerator_factors=angular_mass=0
    for (r,j),expansion in paths(q,d).items():
        angular_mass+=sum(abs(c) for c in expansion.values())*(1+(j>0))
        for term in expansion:
            profile=dict(term)
            assert all(1<=a<=last for a,v in term if v>0)
            assert all(-d<=a<=last+1 for a,v in term if v<0)
            for a in range(1,last+1):
                order=max(order,profile.get(a,0)+2*int(a in poles)+2*int(a-j in poles)+3)
            inverse_factors=max(inverse_factors,7+sum(max(0,v) for a,v in term))
            numerator_factors=max(numerator_factors,sum(max(0,-v) for a,v in term))
    # Every transfer uses a factor from the same half-unit contour about -a.
    # Four simple connections and three radial poles add seven inverse factors.
    last+=displacement
    lattice=lcm(*range(1,last+1))**inverse_factors
    denominator=factorial(lam)**3*factorial(2*lam-1)*lattice
    limit=(lattice*lam**3*20*angular_mass*mass**4*2**inverse_factors*
           (last+d+2)**numerator_factors*(order*last**2+displacement))
    return order,last,denominator,limit


def direct_period(spokes,ring,q,d,fixed_spin=False):
    lam=q-1
    contour,cuts=geometry()
    total=defaultdict(fmpq)
    for (r,j),expansion in paths(q,d).items():
        masks=(0,) if j==0 else (6,9)
        for mask in masks:
            shifted=[j*((mask>>i)&1) for i in range(4)]
            radial={}
            for key,(inside,crossing,internal) in cuts.items():
                extra=r*int(0 in inside)+(2*d-r)*int(2 in inside)+sum(shifted[i] for i in crossing)
                assert extra%2==0
                radial[key]=sum(q-spokes[i] for i in inside)-sum(ring[i] for i in internal)+extra//2
            assert min(radial.values())>0
            for ks in product(*(range(q-b) for b in ring)):
                connected=Counter()
                scale=fmpq(1)
                for b,k,shift in zip(ring,ks,shifted,strict=True):
                    c,roots=connection(q,b,k)
                    scale*=fmpq(str(c))
                    connected.update(a+shift for a in roots)
                rad={key:a+sum(ks[i] for i in cuts[key][1]) for key,a in radial.items()}
                groups=Counter(tuple(sorted(rad[i] for i in row)) for row in contour)
                for profile,number in expansion.items():
                    for roots,multiplicity in groups.items():
                        whole=connected.copy()
                        whole.update(dict(profile))
                        whole.update(roots)
                        offset=(d-shifted[0]) if fixed_spin else 0
                        profile_at_tail=tuple(sorted((a+offset,v) for a,v in whole.items() if v))
                        values=rational_sum(profile_at_tail)
                        if offset:
                            values[0]+=sum(prod(fmpq(n+a)**(-v) for a,v in whole.items() if v)
                                           for n in range(offset))
                        for h,v in values.items():
                            total[h]+=number*multiplicity*scale*v
                    jets_on_pole.cache_clear()
    assert all(h in (0,3,5) or not v for h,v in total.items())
    norm=fmpq(lam**3,factorial(lam)**3*factorial(2*lam-1))
    return tuple(F(str(norm*total[h])) for h in (0,3,5))


def low_dimension(complete,exe,tmp):
    ref=json.loads((HERE.parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    for d in range(2,6):
        direct={}
        forward={}
        inputs=[]
        for ex,vals in zip(ref['exponents'],ref['coefficients'],strict=True):
            seed=tuple(a+1 for a in ex)
            inputs.append((*seed,1,0))
            g=completed(seed,6)
            expected=tuple(F(str(x))*factorial(5)/prod(map(factorial,ex)) for x in vals)
            assert forward.setdefault(min(prism_charts(g,6)),expected)==expected
            for chart in charts(g,6,d):
                assert direct.setdefault(pack(chart),expected)==expected
        direct_count=len(direct)
        selected=reflected(complete,6,inputs,d)
        for seed,(_,k,scale) in selected.items():
            expected=tuple(v/scale for v in forward[seed])
            assert direct.setdefault(k,expected)==expected
        if direct:
            assert certified(6,direct,exe,tmp,degree=d,proof=height(6,d))==direct
        print('independent d, direct charts, reflected seed tests:',d,direct_count,len(selected),flush=True)


EXPECTED={
    3:(-F(364895270148783393589254058194439,10558224741529687198040064000000000),
       -F(84526443564565906262801573,914434304585747033088000000),
       F(110342188952376239,785434536337651200)),
    4:(-F(1791492694432349678304514107083,2277264159937775670165504000000000),
       -F(5772154012645353439573237,2743302913757241099264000000),
       F(12558440591753711,3927172681688256000)),
    5:(-F(65095820170004617900140044417,414787400560094854208716800000000),
       -F(28787513532770928424481,68582572843931027481600000),
       F(62632808252219,98179317042206400))}
COMPLETE_SUBTOTAL=(
    F(223034075807229941528734717303,1009917149189796166769049600000000),
    F(1753641788846304650723,2968942547356321536000000),
    -F(73446572862343,81816097535172000))


def subtotal():
    from mixed_bipartite_wheels import TOTAL
    from mixed_bipartite_quadratic_bank import EXPECTED as quadratic
    bubble=(F(25095295652441839438597,62374045197006744993792000000000),
            -F(383542229667277,1145907649856825856000000),F())
    linear=(-F(21349904253811884693303442588528537,116140472156826559178440704000000000),
            -F(269760817092032201695082203,548660582751448219852800000),
            F(2934584497151629337,3927172681688256000))
    assert tuple(map(sum,zip(bubble,TOTAL,linear,quadratic,*EXPECTED.values())))==COMPLETE_SUBTOTAL
    lo=hi=COMPLETE_SUBTOTAL[0]
    for c,n in zip(COMPLETE_SUBTOTAL[1:],(3,5),strict=True):
        a,b=sorted(c*x for x in zeta_interval(n,terms=96))
        lo+=a
        hi+=b
    assert -112241<10**20*lo<=10**20*hi<-112240
    return COMPLETE_SUBTOTAL


def physical(complete,exe,tmp,check):
    with gzip.open(HERE/'mixed_bipartite_one_star.components.gz','rt') as f:
        data={tuple(map(int,r[:9])):r[:11] for line in f if (r:=line.split())[11]=='0'}
    for kind,n in (('wheel',10940),(1,10583),(2,5903)):
        selected=reflected(complete,13,data.values(),kind)
        assert len(selected)==n
        for seed in selected:del data[seed]
    contributions=[]
    # Each other incident weight is at most q-1, while their sum with -d is 2q.
    for d in range(3,2*13-3):
        selected=reflected(complete,13,data.values(),d)
        for seed in selected:del data[seed]
        if d>5:
            assert not selected
            continue
        keys={k for _,k,_ in selected.values()}
        assert len(selected)==(2583,710,75)[d-3]
        values=certified(13,keys,exe,tmp,degree=d,proof=height(13,d))
        if check:
            key=min(keys,key=lambda k:(prod(13-b for b in unpack(k)[1]),k))
            assert direct_period(*unpack(key),13,d)==values[key]
        value=tuple(sum(w*f*values[k][i] for w,k,f in selected.values()) for i in range(3))
        assert value==EXPECTED[d]
        print('degree, components, columns:',d,len(selected),*value,flush=True)
        contributions.append(value)
    assert len(data)==4769
    if check:
        counts=[]
        for degree in range(1,2*13-3):
            charts=reflected(complete,13,data.values(),degree,formal=True)
            for seed in charts:del data[seed]
            counts.append(len(charts))
        assert counts==[87,232,271,218,109,12]+[0]*16 and len(data)==3840
        print('Remaining with divergent separated terms:',sum(counts),'without a polynomial diagonal:',len(data))
    return tuple(map(sum,zip(*contributions)))


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--controls-only',action='store_true')
    parser.add_argument('--subtotal',action='store_true',help='Check the complete evaluated subtotal and its sign')
    args=parser.parse_args()
    print('32622 disjoint components:',*subtotal(),flush=True)
    if args.subtotal:
        raise SystemExit
    with TemporaryDirectory() as directory:
        tmp=Path(directory)
        complete,exe=tmp/'complete',tmp/'wheel'
        compile_program('quartic_prism_completion.cpp',complete)
        compile_program('quartic_wheel_bank.cpp',exe)
        if args.check or args.controls_only:
            low_dimension(complete,exe,tmp)
        if not args.controls_only:
            print('remaining polynomial diagonals combined:',*physical(complete,exe,tmp,args.check))
