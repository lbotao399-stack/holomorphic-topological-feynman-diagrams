"""Empty-star substitution in the remaining genus-five nine-edge quotient."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from math import lcm,prod

from flint import fmpz_mpoly_ctx
import sympy as sp

from genus_five_nine_series import CTX,G,beta,core_value,loads,numerator,tree_flows

EDGES=((0,1),(0,5),(1,2),(1,3),(2,4),(2,5),(3,4),(3,5),(4,5))
OUTSIDE=fmpz_mpoly_ctx.get(('l12','l14','l24','l25'))


def star_substitute(poly,lift):
    """Multiply by S^lift after replacing star lengths by XY/S, XZ/S, YZ/S."""
    if not poly:return poly
    ctx=poly.context()
    x=ctx.gens()
    groups=defaultdict(dict)
    for ex,c in poly.terms():
        v=list(ex)
        a,b,d=(v[j] for j in (3,6,7))
        v[3],v[6],v[7]=a+b,a+d,b+d
        groups[a+b+d][tuple(v)]=c
    degree=max(groups)
    S=x[3]+x[6]+x[7]
    mult=max(lift,degree)
    num=sum(ctx.from_dict(entries)*S**(mult-k) for k,entries in groups.items())
    result,rem=divmod(num,S**(mult-lift))
    assert not rem
    return result


def groups(poly):
    data=defaultdict(dict)
    for ex,v in poly.terms():
        k=tuple(map(int,ex))
        data[(k[0],k[1],*k[6:])][k[2:6]]=v
    return {k:OUTSIDE.from_dict(row) for k,row in data.items()}


@cache
def path_integral(k):
    a,b,c,d,e=(v+1 for v in k)
    return beta(a,b)*beta(14-a-b,14-c)*beta(14-d,14-e)


def integrate(marks,second):
    graph=G(6,EDGES)
    charge=loads(marks,6)
    trees=tree_flows(graph,charge)
    raw=numerator(graph,charge,next(iter(trees.values())))
    ctx=raw.context()
    x=ctx.gens()
    tree={m:prod(x[i] for i in m) for m in trees}
    current=[[sum((v[e][d]*tree[m] for m,v in trees.items()),ctx.constant(0))
              for d in range(2)] for e in range(9)]
    b=sp.Matrix([[int(v==i)-int(v==j) for i,j in EDGES] for v in range(5)])
    assert tuple(j for j in range(9) if j not in b.rref()[1])==(5,6,7,8)
    u=star_substitute(sum(tree.values()),1)
    h=star_substitute(raw,3)
    j=[[star_substitute(c,1) for c in row] for row in current]
    root=sum(c*c for c in j[3])
    insert=sum((j[7][d]+j[8][d])**2 for d in range(2))
    if second:insert+=2*sum((j[5][d]+j[7][d])**2 for d in range(2))
    # U_Y=U_Delta/S and the absolute star Jacobian is XYZ/S^3. The
    # distinguished resistance is XY/S, so all S factors cancel.
    prefactor=x[3]**2*x[6]**2*x[7]
    for lengths in ((2,3,7,5,11,4,6,3,2),(5,2,11,2,6,3,7,4,1)):
        L=list(map(sp.Integer,lengths))
        X,Y,Z=(L[e] for e in (3,6,7))
        den=X+Y+Z
        L[3],L[6],L[7]=X*Y/den,X*Z/den,Y*Z/den
        cond=sp.diag(*(1/t for t in L))
        tr=cond*b.T*(b*cond*b.T).inv()
        currents=tr*sp.Matrix(charge[:-1])
        H=-(cond-tr*b*cond).multiply_elementwise(currents*currents.T)
        U=int(u(*lengths))
        for i in range(9):
            assert all(currents[i,d]==F(int(j[i][d](*lengths)),U) for d in range(2))
        # Grounded Hessian is an affine minor. Divide by L_root for
        # the projective coefficient, then transform the three star lengths.
        projection=H.minor_submatrix(3,3).det()/L[3]*X*Y*Z/den**3
        numeric=F(int(prefactor(*lengths))*int(h(*lengths))**2,U**10)
        assert projection==numeric
    left=groups(prefactor*h**2)
    right=groups(root*insert)
    sums={tuple(a+b for a,b in zip(p,q)) for p in left for q in right}
    scale=lcm(*(path_integral(k).denominator for k in sums))
    weights={k:int(path_integral(k)*scale) for k in sums}
    reduced=defaultdict(lambda:OUTSIDE.constant(0))
    for p,A in left.items():
        for q,B in right.items():
            k=tuple(a+b for a,b in zip(p,q))
            powers=(sum(k[:3])-12,sum(k[3:])-13)
            reduced[powers]+=weights[k]*A*B
    terms={}
    for (l15,l45),poly in reduced.items():
        assert min(l15,l45)>=0 or not poly
        for ex,c in poly.terms():
            a,b,d,e=map(int,ex)
            terms[a,b,l15,d,e,l45]=c
    core=CTX.from_dict(terms)
    count,total=core_value(core,scale)
    print(marks,'convergent core monomials and orbits',len(core),count,*total,flush=True)
    expected={
        (2,4,0):(F(-750197707817,391877807915962368000000),F(25457,17278562959257600)),
        (4,2,0):(F(69230622026029,15675112316638494720000000),F(-341063,86392814796288000)),
        (0,2,4):(F(-314131519487,69667165851726643200000),F(1877,533288980224000)),
        (0,4,2):(F(-326684370491,489847259894952960000000),F(93371,345571259185152000))}
    assert total==expected[marks]
    return total


def summed_obstructions(channel):
    from genus_five_nine_series import integrate as series
    if channel==1:
        subtotal=(F(145775286798138013,26059874226411497472000000),-F(9817,863928147962880))
        labels=((2,4,0),(4,2,0),(3,4,0))
    elif channel==3:
        subtotal=(F(3963229030865953661,7088285789583927312384000000),-F(343723,34557125918515200))
        labels=((0,2,4),(0,4,2),(0,3,4))
    else:raise ValueError('Expected stable channel 1 or 3.')
    parts=(subtotal,integrate(labels[0],False),integrate(labels[1],True),series(labels[2]))
    total=tuple(map(sum,zip(*parts)))
    expected={1:F(352407233223796349,63175452670088478720000000),
              3:F(1172053333004275717,2147965390783008276480000000)}
    assert total==(expected[channel],0)
    return total[0]


if __name__=='__main__':
    for channel in (1,3):
        print("Rational sum after cancellation:",channel,summed_obstructions(channel),flush=True)
