"""Double pole of the planar six-line master coefficient at physical dimension."""
from fractions import Fraction as F
from functools import cache
from math import comb
from flint import fmpq
from quartic_planar_subcut import AUTOMORPHISMS,cut_polynomial,monomial_jet


@cache
def pochhammer(center,slope,n):
    # Order, leading coefficient, logarithmic derivative of the regular unit.
    if n<0:
        p,a,b=pochhammer(center+n,slope,-n)
        return -p,1/a,-b
    power=0
    lead=F(1)
    log=F(0)
    for j in range(n):
        if center+j:
            lead*=center+j
            log+=F(slope,center+j)
        else:
            power+=1
            lead*=slope
    return power,lead,log


def multiply(factors):
    power,lead,log=0,F(1),F(0)
    for (p,a,b),sign in factors:
        power+=sign*p
        lead=lead*a if sign==1 else lead/a
        log+=sign*b
    return power,lead,log


@cache
def disk_poles(i,j,k,tau,h):
    if i<0 or 0<=j<i:return disk_poles(j,i,k,tau,h)
    slope=2*h+k+j
    c0=F(-1,16)**h*F(1,1+2*h)
    log0=F(-2)+F(2,1+2*h)
    common=[(pochhammer(1+h,-1,i),1),(pochhammer(2+2*h,-2,i+j),-1),
            (pochhammer(0,2,-slope),-1)]
    out=[F(0),F(0)]
    for l in range(i+1):
        power,lead,log=multiply(common+[(pochhammer(1+h+i-l,-1,j),1),
                   (pochhammer(-1,3,-slope-h-tau-l),1),(pochhammer(1,-1,h+tau+l),1)])
        assert power>=-2
        lead*=c0*comb(i,l)
        if power==-2:
            out[0]+=lead
            out[1]+=lead*(log+log0)
        elif power==-1:
            out[1]+=lead
    return tuple(out)


@cache
def hermite_poles(k,t):
    p,q=max(0,-k),max(0,-t)
    order,lead,log=multiply([(pochhammer(0,-2,k),1),(pochhammer(0,-2,t),1),
                             (pochhammer(0,-4,k+t),-1)])
    assert order>=0
    lead*=F(-1)**t
    c0=lead if order==0 else F(0)
    c1=lead*log if order==0 else lead if order==1 else F(0)
    length=max(k+p+t+q,p+q)
    residual=[[F(0) for _ in range(3)] for _ in range(length+1)]
    for a,b,x0,x1 in ((k+p,t+q,1,0),(p,q,-c0,-c1)):
        for j in range(b+1):
            integer=comb(b,j)*(-1)**(b-j)
            residual[a+j][1]+=integer*x0
            residual[a+j][2]+=integer*x1
    solution=[None]*length
    for j in reversed(range(length)):
        const=j-p-q
        neg,zero,one=residual[j+1]
        if const:
            m=neg/const;z=(zero+4*m)/const;u=(one+4*z)/const
        else:
            assert not neg
            m=zero/F(-4);z=one/F(-4);u=0
        solution[j]=(m,z)
        residual[j][0]+=(j-p)*m
        residual[j][1]+=(j-p)*z-2*m
        residual[j][2]+=(j-p)*u-2*z
    assert not any(residual[0][:2])
    return p,q,solution


@cache
def two_poles(k,t,h):
    # Return the epsilon^-3 and epsilon^-2 coefficients of the disk integral
    # after the unit eight-line master has been removed.
    terms=cut_polynomial()[2]
    if h:
        inc=1 if h>0 else -1
        scale=F(1,16)**inc
        kp,tp,hp,factor=(k,t,h-1,1) if h>0 else (k-2,t-2,h,-16)
        triple,double=two_poles(k+2*inc,t+2*inc,h-inc)
        double*=scale
        for (a,b,n),v in terms:
            if a or b:double+=factor*v*disk_poles(a-1,b-1,kp+n,tp,hp)[0]
        return scale*triple,double
    p,q,polynomial=hermite_poles(k,t)
    triple,double=F(0),F(0)
    for degree,(neg,zero) in enumerate(polynomial):
        if not neg and not zero:continue
        for (a,b,n),v in terms:
            if not(a or b):continue
            k0=degree-p+n
            left=disk_poles(a-1,b-1,k0,1-q,-1)
            right=disk_poles(a-1,b-1,k0+1,-q,-1)
            d2,d1=((n-2)*x-2*y for x,y in zip(left,right))
            triple+=F(1,2)*neg*v*d2
            double+=v*((F(1,2)*zero+neg)*d2+F(1,2)*neg*d1)
    return triple,double


def whole_poles(powers,shifts):
    retained=tuple(powers[i]-1 for i in (0,1,2,3,5,6))
    if min(retained)<0:return F(0),F(0)
    ctx,_,terms=cut_polynomial()
    base=ctx.from_dict({ex:fmpq(str(c)) for ex,c in terms})
    dots=sum(retained)
    h=shifts-dots
    combined=ctx.constant(0)
    first=ctx.constant(0)
    for order,poly in monomial_jet(retained):
        c0,c1=F(1),F(0)
        for j in range(order):
            scale=(F(2*shifts-1,2)-j)/(j+1)
            c0,c1=c0*scale,c1*scale-c0/(j+1)
        term=ctx.from_dict({ex:fmpq(str(c)) for ex,c in poly})*base**(dots-order)
        combined+=fmpq(str(c0))*term
        first+=fmpq(str(c1))*term
    shifts_uv=[]
    for axis,power in enumerate((powers[4],powers[7])):
        derivative=base.derivative(axis)
        for n in range(1,power):
            multiplier=fmpq(2*h-1,2)*derivative
            first=(base*first.derivative(axis)+multiplier*first-derivative*combined)/n
            combined=(base*combined.derivative(axis)+multiplier*combined)/n
            h-=1
        shifts_uv.append(-min(power,1))
    triple,double=F(0),F(0)
    for ex,c in combined.terms():
        i,j,k=map(int,ex)
        i+=shifts_uv[0];j+=shifts_uv[1]
        if i==j==-1:
            a,b=two_poles(k,0,h)
            triple+=F(str(c))*a
            double+=F(str(c))*b
        else:
            double+=F(str(c))*disk_poles(i,j,k,0,h)[0]
    # Only the doubly polar monomials can have a third-order pole.
    if shifts_uv==[-1,-1]:
        for ex,c in first.terms():
            i,j,k=map(int,ex)
            if not(i or j):double+=F(str(c))*two_poles(k,0,h)[0]
    norm,log=F(1),F(0)
    for r in range(shifts):
        norm*=F(8,(2*r+1)*(2*r+2)*(2*r+3))
        log+=sum(F(2,2*r+j) for j in (1,2,3))
    return norm*triple,norm*(double+triple*log)


def physical_poles(powers,shifts):
    exchange=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[5],g[8]}=={6,7})
    first=whole_poles(powers,shifts)
    second=whole_poles(tuple(powers[i-1] for i in exchange[1:]),shifts)
    return tuple(x+y for x,y in zip(first,second))


def terminal(folder,executable):
    from pathlib import Path
    import json
    import sympy as sp
    from quartic_projective_ibp import K33,faces,reduction
    from quartic_nonplanar_cut import projection
    from four_loop import form,read_vectors,tagged
    from shift import gram,k33,ordinary_coordinates,vector,FIELD
    source=json.loads(Path(__file__).with_name('quartic_terminal_periods.json').read_text())
    moves,name,keep=reduction(faces(K33,5)[1])
    assert name=='wheel' and not moves
    probes=[]
    for name in ('K33','wheel'):
        for row in source[name]:
            p=row['powers']
            if name=='wheel':
                shortened=[0]*8
                for j,k in enumerate(keep):shortened[k]=p[j]
                p=shortened[:5]+[0]+shortened[5:]
            images={tuple(p[i] for i in g) for g in AUTOMORPHISMS if p[g[0]]>0}
            for chosen in (min(images),max(images)):
                probes.append((name,row,chosen[1:]))
    def expr(p):
        return '*'.join(f'1/(p{j}.p{j})^{a}' for j,a in enumerate(p,1) if a)
    basic=[(i,expr([int(j not in missing) for j in range(1,9)])) for i,missing in enumerate(k33.MISSING)]
    units=read_vectors(form(folder,executable,'six_terminal_basis',tagged(basic),edges=K33))
    basis=[vector(units[i],0,K33) for i in range(5)]
    integral=[(i,expr(p)+'*'+gram(K33)) for i,(_,_,p) in enumerate(probes)]
    raw=read_vectors(form(folder,executable,'six_terminal_gram',tagged(integral),edges=K33))
    ep=sp.Symbol('ep');d=4-2*ep
    result=[]
    for i,(name,row,p) in enumerate(probes):
        coords=ordinary_coordinates(vector(raw[i],0,K33),basis)
        expressions=[sp.cancel(FIELD.to_sympy(q)*8/((d-1)*(d-2)*(d-3))) for q in coords]
        selected=[]
        for c,order in zip(expressions,(0,2,2,1,2)):
            v=sp.cancel(c*ep**order).subs(ep,0)
            assert v.is_Rational
            selected.append(F(v))
        nonplanar,planar,first,second,third=selected
        assert physical_poles(p,1)==(0,planar)
        assert projection(p,1)==nonplanar
        assert planar/3+first/2+third/6==0
        assert 14*planar/3-3*first+2*second-11*third/3==0
        assert 10*nonplanar+60*planar+20*first==F(row['1,zeta3,zeta5'][2])
        result.append({'core':name,'uncut':list(p),'top,planar_double,elementary_double':
                       list(map(str,selected[:3]))})
    return result


if __name__=='__main__':
    from argparse import ArgumentParser
    from tempfile import TemporaryDirectory
    from pathlib import Path
    import json
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('form',type=Path)
    parser.add_argument('forcer',type=Path)
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-planar-terminal-') as temp:
        folder=Path(temp)
        for item in ('forcer','forcer.h'):
            (folder/item).symlink_to((args.forcer/item).resolve(strict=True))
        table=terminal(folder,args.form.resolve(strict=True))
    path=Path(__file__).with_name('quartic_planar_terminal.json')
    content=json.dumps(table,indent=2)+'\n'
    if args.save:path.write_text(content)
    else:assert path.read_text()==content
    print('Recovered all zeta(5) contributions from three rational projections in',len(table),'cuts.',flush=True)
