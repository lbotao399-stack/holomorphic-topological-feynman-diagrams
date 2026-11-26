"""Finite sampling kernel for every positive four-wheel period in dimension 2q."""
from fractions import Fraction as F
from functools import cache
from math import comb,factorial,prod
from quartic_bipartite_wheel import binom_integer,harmonic_finite
from flint import fmpq as Q


def poles(q):
    lam=q-1
    return {a:(5 if a==lam else 6 if a<=2*lam-1 else 3) for a in range(1,3*q-4)}


def inverse_jets(a,through,roots):
    v=[Q(1)]+[Q(0)]*(through-1)
    for b,e in roots:
        step=Q(b-a)
        coeff=[Q(binom_integer(-e,i))*step**(-e-i) for i in range(through)]
        v=[sum((v[j]*coeff[i-j] for j in range(i+1)),Q(0)) for i in range(through)]
    return v


def kernel(q):
    shifts=poles(q)
    count=sum(shifts.values())-1
    common={a:inverse_jets(a,m,[(b,e) for b,e in shifts.items() if a!=b]+[(-v,-1) for v in range(count)])
            for a,m in shifts.items()}
    results=[]
    for node in range(count):
        scale=Q((-1)**(count-1-node),factorial(node)*factorial(count-1-node))
        scale*=prod(Q(node+a)**m for a,m in shifts.items())
        values=[Q(0)]*7
        for a,m in shifts.items():
            # Divide the common numerator by n-node to remove the missing cardinal factor.
            inv=Q(-1,a+node)
            leading=Q(0)
            for j,coefficient in enumerate(common[a]):
                leading=(coefficient-leading)*inv
                value=scale*leading
                power=m-j
                values[power]+=value
                values[0]-=value*harmonic_finite(a,power)
        assert not values[1]
        results.append(values)
    return results


def modular_kernel(q,p):
    # Rational interpolation is done once over the requested finite field.
    exponents=poles(q)
    count=sum(exponents.values())-1
    def inverse(k):return pow(k%p,-1,p)
    common={}
    for a,m in exponents.items():
        v=[1]+[0]*(m-1)
        for b,e in [(b,e) for b,e in exponents.items() if b!=a]+[(-j,-1) for j in range(count)]:
            coeff=[binom_integer(-e,k)*pow((b-a)%p,-e-k,p)%p for k in range(m)]
            v=[sum(v[j]*coeff[k-j] for j in range(k+1))%p for k in range(m)]
        common[a]=v
    harmonic={(a,k):sum(pow(i,-k,p) for i in range(1,a))%p for a,m in exponents.items() for k in range(1,m+1)}
    out=[]
    for j in range(count):
        scale=(-1)**(count-j-1)*inverse(factorial(j)*factorial(count-1-j))%p
        for a,m in exponents.items():scale=scale*pow(j+a,m,p)%p
        row=[0]*7
        for a,m in exponents.items():
            inv=inverse(-a-j);last=0
            for k,v in enumerate(common[a]):
                last=(v-last)*inv%p
                value=last*scale%p;order=m-k
                row[order]=(row[order]+value)%p
                row[0]=(row[0]-value*harmonic[a,order])%p
        assert not row[1]
        out.append(row)
    return out


if __name__=='__main__':
    from argparse import ArgumentParser
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('q',type=int)
    cli.add_argument('prime',type=int)
    args=cli.parse_args()
    for row in modular_kernel(args.q,args.prime):print(*row)
