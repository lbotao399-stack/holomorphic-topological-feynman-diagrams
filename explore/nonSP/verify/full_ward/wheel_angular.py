"""Exact projective wheel monomials at any loop order, by angular harmonics."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial, prod

from flint import fmpq, fmpq_poly

X=fmpq_poly([0,1])


def rational(a):
    return fmpq(a.numerator,a.denominator)


def multiplicity(q):
    s=X
    return 2*s**2*prod((s**2-j*j for j in range(1,q-1)),start=fmpq_poly([1]))/factorial(2*q-2)


def radial_master(n,q):
    """Powers q-1 on the spokes and 1 on the rims; keys are zeta weights."""
    if not 2<=q<n:
        raise ValueError('The positive-index radial master converges for 2 <= q < n.')
    return {2*n-1-j:Fraction(str(c))*comb(2*n-2,n-1)
            for j,c in enumerate(multiplicity(q)) if c}


def sectors(spoke,rim,q):
    n=len(rim)
    drift=[2*q-2*spoke[i]-rim[i-1]-rim[i] for i in range(n)]
    assert sum(drift)==0
    cuts={}
    for bit in range(1,(1<<n)-1):
        inside=[(bit>>i)&1 for i in range(n)]
        cross=tuple(i for i in range(n) if inside[i]!=inside[(i+1)%n])
        a=sum(v*w for v,w in zip(inside,drift))
        if sum(rim[i] for i in cross)+a<=0:
            raise ValueError('Divergent radial sector.')
        cuts[bit]=(cross,a)
    schedules=[]
    for order in permutations(range(n)):
        bit=0
        chain=[]
        for i in order[:-1]:
            bit|=1<<i
            chain.append(bit)
        weight=Fraction(1,prod(len(cuts[bit][0]) for bit in chain))
        schedules.append((chain,weight))
    return cuts,schedules


def harmonic_sum(spoke,rim,q):
    """Return the reduced odd rational harmonic summand and its integer poles."""
    n=len(spoke)
    if len(rim)!=n or n<3 or sum(spoke)+sum(rim)!=n*q:
        raise ValueError('Expected n spoke powers and n rim powers, of total n*q.')
    if any(type(k) is not int or not 1<=k<q for k in (*spoke,*rim)):
        raise ValueError('Convergent projective powers are integers between one and q-1.')
    if sum(rim)<=q:
        raise ValueError('Divergent coincidence of all rim vertices away from the hub.')
    # Position-space coordinate exponents.
    spoke=tuple(q-j for j in spoke)
    rim=tuple(q-j for j in rim)
    d=[q-1-r for r in rim]
    cuts,schedules=sectors(spoke,rim,q)
    terms=Counter()
    for k in product(*(range(v+1) for v in d)):
        coefficient=prod((-1)**j*comb(v,j) for v,j in zip(d,k))
        base=[Fraction(j-v+z) for v,j in zip(d,k) for z in range(v+1)]
        shifts=[2*j-v for j,v in zip(k,d)]
        linear={bit:Fraction(a+sum(shifts[i] for i in cross),len(cross))
                for bit,(cross,a) in cuts.items()}
        for chain,weight in schedules:
            poles=tuple(sorted(base+[linear[bit] for bit in chain]))
            terms[poles]+=coefficient*weight
    terms={e:c for e,c in terms.items() if c}
    common=Counter()
    for poles in terms:
        common |= Counter(poles)
    factors={p:X+rational(p) for p in common}
    denominator=prod((factors[p]**v for p,v in common.items()),start=fmpq_poly([1]))
    numerator=fmpq_poly([])
    for poles,c in terms.items():
        term=prod((factors[p]**v for p,v in Counter(poles).items()),start=fmpq_poly([1]))
        numerator+=rational(c)*(denominator//term)
    numerator*=2**(n-1)*prod(factorial(v-1) for v in spoke)*multiplicity(q)
    cancel=numerator.gcd(denominator)
    numerator//=cancel
    denominator//=cancel
    assert numerator.degree()+2<=denominator.degree()
    assert numerator(-X)*denominator==-denominator(-X)*numerator
    poles={}
    for p in common:
        v=0
        while denominator((-rational(p)))==0:
            denominator//=factors[p]
            v+=1
        if v:
            assert p.denominator==1 and v<=2*n-1 and q-1+p>0
            poles[int(p)]=v
    assert denominator==1
    return numerator,poles


def integrate(spoke,rim,q):
    """Dictionary with key zero for the rational term; positive keys label zeta values."""
    n=len(spoke)
    num,poles=harmonic_sum(spoke,rim,q)
    zeta=Counter()
    for p,order in poles.items():
        # Taylor expansion at s=-p after removing its singular factor.
        other=prod(((X+h-p)**a for h,a in poles.items() if h!=p),start=fmpq_poly([1]))
        shifted=num(X-p)
        coefficients=[]
        for j in range(order):
            a=(shifted[j]-sum(other[k]*coefficients[j-k] for k in range(1,j+1)))/other[0]
            coefficients.append(a)
            weight=order-j
            c=Fraction(str(a))
            zeta[weight]+=c
            zeta[0]-=c*sum((Fraction(1,k**weight) for k in range(1,q-1+p)),Fraction())
    assert not zeta[1]
    assert all(not c or w==0 or (w%2==1 and w<=2*n-3) for w,c in zeta.items())
    return dict(sorted((w,c) for w,c in zeta.items() if c))


def cut_rows():
    import json
    from pathlib import Path
    bank=json.loads((Path(__file__).resolve().parents[1]/'prism_forcer'/'wheel.json').read_text())
    edges=list(map(tuple,bank['edges']))
    vertices=(0,1,4,3)
    sp=[edges.index(tuple(sorted((2,i)))) for i in vertices]
    rim=[edges.index(tuple(sorted((i,j)))) for i,j in zip(vertices,vertices[1:]+vertices[:1])]
    for a,values in bank['cut_integrals']:
        scale=Fraction(prod(factorial(j-1) for j in a),factorial(5))
        yield [a[j] for j in sp],[a[j] for j in rim],[scale*Fraction(v) for v in values]


def check_highest(rows):
    import os
    import subprocess
    from pathlib import Path
    from tempfile import TemporaryDirectory
    prime=1000003
    def check(binary,n,q,data):
        stream=[f'{n} {q} {prime}']
        for spoke,rim,value in data:
            stream.append(' '.join(map(str,[*spoke,*rim])))
        got=subprocess.check_output([binary],input='\n'.join(stream)+'\n',text=True).split()
        assert list(map(int,got))==[v.numerator*pow(v.denominator,-1,prime)%prime for _,_,v in data]
    with TemporaryDirectory() as temp:
        binary=Path(temp)/'wheel_highest'
        subprocess.run([os.environ.get('CXX','c++'),'-std=c++17','-O3',
                        Path(__file__).with_name('wheel_highest.cpp'),'-o',binary],check=True)
        check(binary,4,6,[(s,r,v[-1]) for s,r,v in rows])
        data=[]
        for spoke,rim in [([5,5,4,5,5],[3,2,2,2,2]),([6,4,5,5,4],[2,2,1,3,3]),
                          ([6,5,4,4,5],[1,3,2,2,3]),([3,5,5,5,3],[3,3,4,2,2]),
                          ([6,6,4,6,4],[2,2,1,2,2]),([4,4,4,4,4],[3,3,3,3,3])]:
            data.append((spoke,rim,integrate(spoke,rim,7).get(7,Fraction())))
        assert data[0][-1]==-Fraction(7,594)
        check(binary,5,7,data)
    print('Highest-weight residues match all 1421 four-loop cuts and six five-loop angular evaluations.')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true',help='check all four-loop cuts with complete angular summands')
    parser.add_argument('--highest',action='store_true',help='also compile and check the fast highest-weight evaluator')
    args=parser.parse_args()
    for n in range(3,8):
        for q in range(2,n):
            assert integrate([q-1]*n,[1]*n,q)==radial_master(n,q)
    print('Radial wheel periods for n=3 through 7 satisfy the closed formula.',flush=True)
    rows=list(cut_rows())
    selected=range(len(rows)) if args.all else [0,1,73,404,710,1213,1420]
    for i in selected:
        spoke,rim,expected=rows[i]
        exact=integrate(spoke,rim,6)
        assert [exact.get(k,0) for k in (0,3,5)]==expected
    print(len(selected),'wheel period comparisons with the independent four-loop cuts.')
    if args.highest:
        check_highest(rows)
