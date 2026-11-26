"""Six-residue projector onto the nontrivial planar lower master of a K3,3 cut."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import comb,factorial,prod
from pathlib import Path
from tempfile import TemporaryDirectory
import sympy as sp

from flint import fmpq_mpoly_ctx
from quartic_nonplanar_cut import AUTOMORPHISMS,endpoint_gram
from quartic_projective_ibp import K33
from four_loop import form,read_vectors,tagged
from shift import gram


@cache
def rising(x,n):
    if n>=0:return prod((x+j for j in range(n)),start=F(1))
    return F(1)/prod(x-j for j in range(1,1-n))


@cache
def disk_moment(i,j,k,tau,h,s):
    """Moment u^i v^j w^k B^(s+h) divided by the cut of B^s."""
    if i<0 and j>=0:i,j=j,i
    if j>=0 and i>j:i,j=j,i
    if i<0:raise ValueError('The one-contact moment requires one nonnegative power.')
    e=s+F(3,2)
    a=-3*s-F(5,2)
    slope=2*h+k+j
    base=F(-1,16)**h*(s+1)/(s+h+1)*rising(e+h,i)
    base/=rising(2*e+2*h,i+j)*rising(a+e,-slope)
    # Given v=w*x, the conditional mean of u is x+(w-1)(1-x).
    # Its transverse even moments and Chu-Vandermonde collapse to this sum.
    return base*sum(comb(i,l)*rising(e+h+i-l,j)*rising(a,-slope-h-tau-l)*rising(e,h+tau+l)
                    for l in range(i+1))


@cache
def cut_polynomial():
    gram16=endpoint_gram(0)
    ctx=fmpq_mpoly_ctx.get(('u','v','w'))
    pieces=defaultdict(dict)
    for ex,c in gram16.terms():
        ex=tuple(map(int,ex))
        retained=tuple(ex[i] for i in (0,1,2,3,5,6))
        pieces[retained][(ex[4],ex[7],ex[8])]=c/16
    pieces={ex:ctx.from_dict(coeffs) for ex,coeffs in pieces.items()}
    u,v,w=ctx.gens()
    base=pieces.pop((0,)*6)
    s=u+v
    d=u-v
    assert 16*base==w*(d*d+(w-1)*(s-w)**2-w*(w-1))
    return ctx,pieces,tuple((tuple(map(int,ex)),F(str(c))) for ex,c in base.terms())


@cache
def monomial_jet(target):
    """Coefficients of each power of B-B0 in the six cut variables."""
    ctx,all_parts,_=cut_polynomial()
    parts={ex:p for ex,p in all_parts.items() if all(x<=y for x,y in zip(ex,target))}
    @cache
    def build(ex,order):
        if order==0:return ctx.constant(int(not any(ex)))
        if order>sum(ex) or 4*order<sum(ex):return ctx.constant(0)
        return sum((p*build(tuple(y-x for y,x in zip(ex,beta)),order-1) for beta,p in parts.items()
                    if all(x<=y for x,y in zip(beta,ex))),ctx.constant(0))
    return [(order,[(tuple(map(int,r)),F(str(c))) for r,c in poly.terms()])
            for order in range(sum(target)+1) if (poly:=build(target,order))]


@cache
def hermite(k,t,s):
    p,q=max(0,-k),max(0,-t)
    coefficient=F(-1)**t*rising(2*s+1,k)*rising(2*s+1,t)/rising(4*s+2,k+t)
    length=max(k+p+t+q,p+q)
    residual=[0]*(length+1)
    for a,b,weight in ((k+p,t+q,1),(p,q,-coefficient)):
        for j in range(b+1):residual[a+j]+=weight*comb(b,j)*(-1)**(b-j)
    polynomial=[0]*length
    for j in reversed(range(length)):
        polynomial[j]=residual[j+1]/(j+2-p-q+4*s)
        residual[j]+=(j+1-p+2*s)*polynomial[j]
    assert not residual[0]
    return coefficient,p,q,polynomial


@cache
def moment(i,j,k,t,h,s):
    if i>=0 or j>=0:return 0,disk_moment(i,j,k,t,h,s)
    # The two square variables have the same Gram polynomial at this cut.
    if i>j:return moment(j,i,k,t,h,s)
    terms=cut_polynomial()[2]
    if i < -1:
        # Move the derivative off u^(i+1), leaving only higher u powers.
        factor=-(s+h)/(i+1)
        pair=[0,0]
        for (a,b,n),v in terms:
            if not a:continue
            lower=moment(i+a,j+b,k+n,t,h-1,s)
            for z in (0,1):pair[z]+=factor*a*v*lower[z]
        return tuple(pair)
    if h:
        if h>0:kp,tp,hp,factor=k,t,h-1,F(1)
        else:kp,tp,hp,factor=k-2,t-2,h,F(-16)
        increment=int(h>0)-int(h<0)
        top,lower=moment(i,j,k+2*increment,t+2*increment,h-increment,s)
        core=F(1,16)**increment
        lower*=core
        for (a,b,n),v in terms:
            if a or b:lower+=factor*v*disk_moment(a-1,b-1,kp+n,tp,hp,s)
        return top*core,lower
    top,p,q,polynomial=hermite(k,t,s)
    lower=0
    for degree,weight in enumerate(polynomial):
        if not weight:continue
        for (a,b,n),v in terms:
            if not(a or b):continue
            k0=degree-p+n
            remnant=(n-2)*disk_moment(a-1,b-1,k0,1-q,-1,s)
            remnant-=2*disk_moment(a-1,b-1,k0+1,-q,-1,s)
            lower-=s*weight*v*remnant
    return top,lower


def one_cut(powers,shifts,s):
    retained=tuple(powers[i]-1 for i in (0,1,2,3,5,6))
    if min(retained)<0:return F(0),F(0)
    a,b=powers[4],powers[7]
    total=[F(0),F(0)]
    for level,poly in monomial_jet(retained):
        h=shifts-level
        coefficient=prod(((s+shifts-j)/(j+1) for j in range(level)),start=F(1))
        for (i,j,k),weight in poly:
            value=moment(i-a,j-b,k,0,h,s)
            for z in (0,1):total[z]+=coefficient*weight*value[z]
    d=2*s+5
    scale=8**shifts/prod(((d+2*j-1)*(d+2*j-2)*(d+2*j-3) for j in range(shifts)),start=F(1))
    return tuple(scale*v for v in total)


def whole_cut(powers,shifts,s):
    """Transfer both remaining propagator derivatives before applying the moment."""
    retained=tuple(powers[i]-1 for i in (0,1,2,3,5,6))
    if min(retained)<0:return F(0),F(0)
    ctx,_,terms=cut_polynomial()
    from flint import fmpq
    base=ctx.from_dict({ex:fmpq(str(c)) for ex,c in terms})
    dots=sum(retained)
    h=shifts-dots
    combined=ctx.constant(0)
    for order,poly in monomial_jet(retained):
        coeff=prod(((s+shifts-j)/(j+1) for j in range(order)),start=F(1))
        combined+=fmpq(str(coeff))*ctx.from_dict({ex:fmpq(str(c)) for ex,c in poly})*base**(dots-order)
    exponents=[]
    for index,power in ((0,powers[4]),(1,powers[7])):
        for n in range(1,power):
            combined=(base*combined.derivative(index)+fmpq(str(s+h))*base.derivative(index)*combined)/n
            h-=1
        exponents.append(-min(power,1))
    top,lower=F(0),F(0)
    di,dj=exponents
    for ex,c in combined.terms():
        i,j,k=map(int,ex)
        u,v=moment(i+di,j+dj,k,0,h,s)
        w=F(str(c))
        top+=w*u
        lower+=w*v
    d=2*s+5
    prefactor=8**shifts/prod(((d+2*j-1)*(d+2*j-2)*(d+2*j-3) for j in range(shifts)),start=F(1))
    return top*prefactor,lower*prefactor


def subcut(powers,shifts,s,apply=one_cut):
    # The only six-line cuts with a nonsquare restricted Gram are the two
    # images of the missing pair {5,8}. The ordinary basis uses the first.
    exchange=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[5],g[8]}=={6,7})
    first=apply(powers,shifts,s)
    second=apply(tuple(powers[i-1] for i in exchange[1:]),shifts,s)
    assert first[0]==second[0]
    return first[0],first[1]+second[1]


def controls(folder,executable):
    examples=[(1,1,1,1,0,1,1,0),(1,1,1,1,0,1,1,1),(1,1,1,1,2,1,1,0),
              (2,1,1,1,-1,1,1,2),(1,2,1,2,1,0,2,1),(1,1,1,1,0,0,1,1),
              (0,1,1,1,1,1,1,1),(2,1,2,1,1,1,0,1),(1,)*8,(1,2,1,1,2,1,2,1),
              (2,1,1,1,2,2,1,1)]
    def expr(p):
        return '*'.join(f'1/(p{j}.p{j})^{v}' if v>0 else f'(p{j}.p{j})^{abs(v)}'
                        for j,v in enumerate(p,1) if v)
    master='GschemeConstants(0,0)*BasicT1Integral'
    units=read_vectors(form(folder,executable,'six_unit',tagged([(0,expr(examples[0])),(1,expr((1,)*8))]),edges=K33))
    def rational(pair):
        return sp.cancel(sp.sympify(pair[0].replace('^','**'))/sp.sympify(pair[1].replace('^','**')))
    unit=rational(units[0][master])
    lower_in_top=rational(units[1].get(master,('0','1')))
    top_unit=rational(units[1]['BasicNOIntegral'])
    ep=sp.Symbol('ep')
    for shift in (0,1):
        terms=[(i,expr(p)+(f'*{gram(K33)}' if shift else '')) for i,p in enumerate(examples)]
        rows=read_vectors(form(folder,executable,f'six_cut_{shift}',tagged(terms),edges=K33))
        for i,p in enumerate(examples):
            top=rational(rows[i].get('BasicNOIntegral',('0','1')))/top_unit
            supplied=(rational(rows[i].get(master,('0','1')))-lower_in_top*top)/unit
            for s in (F(-8,11),F(2,7),F(5,13)):
                d=2*s+5
                weight=F(8**shift,prod((d+2*j-1)*(d+2*j-2)*(d+2*j-3) for j in range(shift)))
                expected=tuple(F(v.subs(ep,F(2)-d/2))*weight for v in (top,supplied))
                actual=subcut(p,shift,s)
                assert actual==expected,(p,shift,s,actual,expected)
                assert subcut(p,shift,s,whole_cut)==expected
    print('66 independent simultaneous comparisons of two master coefficients',flush=True)


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('form',type=Path)
    parser.add_argument('forcer',type=Path)
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-planar-subcut-') as temp:
        folder=Path(temp)
        for item in ('forcer','forcer.h'):
            (folder/item).symlink_to((args.forcer/item).resolve(strict=True))
        controls(folder,args.form.resolve(strict=True))
