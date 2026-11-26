"""All weight channels of a four-wheel from a single Gegenbauer index."""
from collections import Counter,defaultdict
from functools import cache
from fractions import Fraction as F
from itertools import permutations,product
from math import comb,factorial,prod
from flint import fmpq as Q

from quartic_bipartite_sums import connection
from quartic_projective_ibp import K33,faces


def dual_powers(edges,powers):
    present=[(u,v,a) for (u,v),a in zip(edges,powers) if a]
    degree=Counter(u for u,v,a in present)+Counter(v for u,v,a in present)
    center=next(v for v,c in degree.items() if c==4)
    rim=sorted(v for v in degree if v!=center)
    weights={tuple(sorted((u,v))):a for u,v,a in present}
    order=next(p for p in permutations(rim) if all(tuple(sorted((p[i-1],p[i]))) in weights for i in range(4)))
    # Dual face i has original rim vertices order[i],order[i+1].
    spokes=tuple(weights[tuple(sorted((order[i],order[(i+1)%4])))] for i in range(4))
    ring=tuple(weights[tuple(sorted((center,order[i])))] for i in range(4))
    return spokes,ring


@cache
def geometry():
    contour=[]
    for i,j,k,l in permutations(range(4)):
        one=1<<i;two=(1<<i)|(1<<j);three=two|(1<<k)
        if (i-j)%2:
            contour.append((one,two,three))
        elif i<j:
            # No rim edge joins i and j. B(two)=B(i)+B(j).
            contour.append((one,1<<j,three))
    assert len(contour)==20
    cuts={}
    for mask in {v for row in contour for v in row}:
        vertices=tuple(i for i in range(4) if mask>>i&1)
        crossing=tuple(i for i in range(4) if ((mask>>(i-1 if i else 3))^(mask>>i))&1)
        internal=tuple(i for i in vertices if mask>>((i-1)%4)&1)
        assert len(crossing)==2
        cuts[mask]=vertices,crossing,internal
    return contour,cuts


def binom_integer(p,k):
    if p>=0:return comb(p,k) if k<=p else 0
    return (-1)**k*comb(k-p-1,k)


@cache
def jets_on_pole(profile,a,through):
    # Product of (n+b)^(-p) for b!=a, near n=-a.
    values=[Q(1)]+[Q(0)]*(through-1)
    for b,p in profile:
        if b==a:continue
        diff=Q(b-a)
        jets=[Q(binom_integer(-p,j))*diff**(-p-j) for j in range(through)]
        values=[sum((values[i]*jets[j-i] for i in range(j+1)),Q(0)) for j in range(through)]
    return values


@cache
def harmonic_finite(a,p):
    return sum((1/Q(n)**p for n in range(1,a)),Q(0))


def rational_sum(profile):
    """Hadamard constant with polynomial and simple-harmonic parts accumulated separately."""
    result=defaultdict(Q)
    for a,p in profile:
        if p<=0:continue
        for j,c in enumerate(jets_on_pole(profile,a,p)):
            k=p-j
            result[k]+=c
            result[0]-=c*harmonic_finite(a,k)
    # Constant finite sum of n^d, from n=0 to N-1, is zero for nonnegative d.
    return result


def period(spokes,ring,q):
    if sum(spokes)+sum(ring)!=4*q or max(spokes+ring)>=q or min(spokes+ring)<=0:
        raise ValueError('Positive convergent projective wheel powers are required.')
    lam=q-1
    contour,cuts=geometry()
    radial={key:sum(q-spokes[i] for i in vs)-sum(ring[i] for i in inside)
            for key,(vs,cross,inside) in cuts.items()}
    total=defaultdict(Q)
    for ks in product(*(range(q-b) for b in ring)):
        base=Counter({i:-1 for i in range(1,2*lam)})
        base[lam]+=3;weight=Q(1)
        for b,k in zip(ring,ks):
            c,roots=connection(q,b,k)
            base.update(roots);weight*=Q(str(c))
        cut={key:d+sum(ks[i] for i in cuts[key][1]) for key,d in radial.items()}
        assert min(cut.values())>0
        factors=Counter(tuple(sorted(cut[i] for i in flags)) for flags in contour)
        for roots,multiplicity in factors.items():
            coeffs=base.copy();coeffs.update(roots)
            profile=tuple(sorted((a,p) for a,p in coeffs.items() if p))
            summation=rational_sum(profile)
            for p,c in summation.items():total[p]+=multiplicity*weight*c
        # The pole coefficients are used once per polynomial kernel.
        jets_on_pole.cache_clear()
    assert not total[1]
    norm=Q(lam**3,factorial(q-1)**3*factorial(2*lam-1))
    values={k:c*norm for k,c in total.items() if c}
    assert all(k in (0,3,5) for k in values)
    return values


def terminal_controls():
    from pathlib import Path
    import json
    from four_loop import WHEEL
    rows=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())['wheel']
    for row in rows:
        p=row['powers'];u,b=dual_powers(WHEEL,p)
        values=period(u,b,3)
        weight=Q(prod(factorial(a-1) for a in p),2)
        assert [weight*values.get(k,0) for k in (0,3,5)]==list(map(Q,row['1,zeta3,zeta5'])),(p,u,b,values)
    print('All nine four-wheel boundary periods recovered from a single angular mode sum.',flush=True)


if __name__=='__main__':
    terminal_controls()
