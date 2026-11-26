"""Prism period by eliminating the common three-loop master at D=11."""
from collections import defaultdict
from functools import cache
from itertools import combinations
from pathlib import Path
import re

import sympy as sp
from flint import fmpq, fmpq_poly

HERE=Path(__file__).resolve().parent
MASTER='GschemeConstants(0,0)*BasicT1Integral'
G=[f'GschemeConstants({i},0)' for i in range(3)]
BUBBLES={f'{G[0]}^3':(0,0,0),f'{G[0]}^2*{G[1]}':(0,0,1),
         f'{G[0]}^2*{G[2]}':(0,0,2),f'{G[0]}*{G[1]}*{G[2]}':(0,1,2)}
HALF=sp.Rational(1,2)


class Rational:
    """A rational function of t=ep+7/2 with exact Taylor/Laurent extraction."""
    def __init__(self,numerator,denominator=1):
        n,d=fmpq_poly(numerator),fmpq_poly(denominator)
        g=n.gcd(d)
        self.numer,self.denom=n//g,d//g

    def __add__(self,other):
        return Rational(self.numer*other.denom+other.numer*self.denom,self.denom*other.denom)

    def __sub__(self,other):
        return Rational(self.numer*other.denom-other.numer*self.denom,self.denom*other.denom)

    def __mul__(self,other):
        return Rational(self.numer*other.numer,self.denom*other.denom)

    def __truediv__(self,other):
        return Rational(self.numer*other.denom,self.denom*other.numer)

    def valuation(self):
        if not self.numer:
            return sp.oo
        return next(i for i,c in enumerate(self.numer) if c)-next(i for i,c in enumerate(self.denom) if c)

    def series(self,last):
        if not self.numer:
            return {}
        a,b=self.numer.coeffs(),self.denom.coeffs()
        n=next(i for i,c in enumerate(a) if c)
        m=next(i for i,c in enumerate(b) if c)
        first=n-m
        a,b=a[n:],b[m:]
        result=[]
        for k in range(last-first+1):
            c=a[k] if k<len(a) else fmpq(0)
            for j in range(1,min(k+1,len(b))):
                c-=b[j]*result[k-j]
            result.append(c/b[0])
        return {first+i:sp.Rational(str(c)) for i,c in enumerate(result) if c}


def parse_coefficients(text):
    text=''.join(text.split()).replace('\\','')
    if text=='0':
        return {}
    names=[MASTER,*BUBBLES]
    for j,key in enumerate(names):
        text=text.replace(key,f'm{j}')
    pattern=r'm([0-4])\*rat\(([-+ep^*0-9]+),([-+ep^*0-9]+)\)'
    if not re.fullmatch(pattern+'(?:\\+'+pattern+')*',text):
        raise ValueError('unexpected FORCER master coefficient')
    result={}
    t=fmpq_poly([fmpq(-7,2),1])
    def polynomial(s):
        # Only integer monomials in ep enter from the validated FORM grammar.
        terms=re.findall(r'[+-]?[^+-]+',s)
        answer=fmpq_poly()
        for term in terms:
            if 'ep' not in term:
                answer+=int(term)
            else:
                c,_,exp=term.partition('ep')
                c=c.rstrip('*')
                v=-1 if c=='-' else (1 if c in ('+','') else int(c))
                answer+=v*t**(int(exp[1:]) if exp else 1)
        return answer
    for index,n,d in re.findall(pattern,text):
        name=names[int(index)]
        assert name not in result
        result[name]=Rational(polynomial(n),polynomial(d))
    return result


@cache
def coefficients(folder,cut,rank):
    return parse_coefficients((folder/f'c{cut}-n{rank}.out').read_text())


@cache
def bubble_product(js):
    """Return power, leading coefficient and first two regularized Taylor coefficients."""
    factor=sp.prod(j+1 for j in js)/sp.gamma(11*HALF)
    first=second=order=0
    for j in js:
        # Normalization (j+1) G(1+j*ep,1); normalized integration measure d^Dk/pi^(D/2).
        z=11*HALF
        a=1-7*HALF*j
        args=[(z-a,-1-j,1),(z-1,-1,1),(a+1-z,j+1,1),(a,j,-1),(2*z-a-1,-2-j,-1)]
        for argument,slope,power in args:
            base=sp.S.One if argument.is_integer else HALF
            factor*=sp.gamma(base)**power
            first+=power*slope*(-sp.EulerGamma if base==1 else -sp.EulerGamma-2*sp.log(2))
            second+=power*slope**2*(sp.pi**2/6 if base==1 else sp.pi**2/2)/2
            step=int(argument-base)
            for k in (range(step) if step>=0 else range(step,0)):
                x=base+k
                sign=power if step>=0 else -power
                if x==0:
                    assert slope
                    order+=sign
                    factor*=sp.Integer(slope)**sign
                else:
                    factor*=x**sign
                    first+=sign*sp.Rational(slope,x)
                    second-=sign*sp.Rational(slope**2,2*x*x)
    return order,sp.factor(factor),(1,sp.expand(first),sp.expand(second+first**2/2))


def cut_pair(rank,cut0,cut1,folder):
    p,q=(coefficients(folder,cut,rank) for cut in (cut0,cut1))
    a=p[MASTER]/(p[MASTER]-q[MASTER])
    b=Rational(1)-a
    assert a.valuation()>=0
    assert not (b*p[MASTER]+a*q[MASTER]).numer
    return bubble_limit({name:b*p.get(name,Rational(0))+a*q.get(name,Rational(0)) for name in BUBBLES})


def bubble_limit(bubble_weights):
    values=defaultdict(int)
    for name,js in BUBBLES.items():
        c=bubble_weights[name]
        order,leading,taylor=bubble_product(js)
        assert c.valuation()+order>=-2
        for i,v in c.series(-order).items():
            for j,t in enumerate(taylor):
                if i+j+order<=0:
                    values[i+j+order]+=v*leading*t
    assert all(sp.expand(v)==0 for k,v in values.items() if k<0)
    return sp.factor(sp.expand(values[0]))


@cache
def auxiliary(cut0,cut1,folder):
    first,second=(coefficients(folder,cut,4) for cut in (cut0,cut1))
    return {name:first.get(name,Rational(0))-second.get(name,Rational(0))
            for name in {MASTER,*BUBBLES}}


def period_from_vector(vector,auxiliary_cuts=(0,1),folder=HERE/'masters'):
    """Evaluate a new convergent critical prism combination with the stored cut difference."""
    zero=auxiliary(*auxiliary_cuts,folder)
    weight=vector.get(MASTER,Rational(0))/zero[MASTER]
    if weight.valuation()<0:
        raise ValueError('This cut needs an auxiliary difference with a higher pole order.')
    return bubble_limit({name:vector.get(name,Rational(0))-weight*zero[name] for name in BUBBLES})


def prism_coefficients(folder=HERE/'masters'):
    """Map each power of N to its integrated graph sum, before the common heat-kernel factor."""
    result={}
    for rank in (0,2,4):
        values={cut_pair(rank,i,j,folder) for i,j in combinations(range(3),2)}
        assert len(values)==1
        result[rank]=values.pop()
    return result


if __name__=='__main__':
    print(prism_coefficients())
