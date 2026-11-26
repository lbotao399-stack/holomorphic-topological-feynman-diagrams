"""Four local double residues for the leading planar master pole."""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import comb,prod

from flint import fmpq, fmpq_mpoly_ctx
from quartic_nonplanar_cut import AUTOMORPHISMS,endpoint_gram,projection


@cache
def node_parts(axis):
    ctx=fmpq_mpoly_ctx.get(tuple(f'y{i}' for i in range(8))+('sigma','tau'))
    variables=list(ctx.gens())
    sigma,tau=variables[-2:]
    arguments=variables[:8]
    orientation = 1 if axis < 4 else -1
    arguments[axis] += orientation*(1+sigma+tau)
    node=endpoint_gram(0).compose(*arguments,-tau)
    ctx2=fmpq_mpoly_ctx.get(('sigma','tau'))
    blocks=defaultdict(dict)
    for exponent,coefficient in node.terms():
        key=tuple(map(int,exponent))
        blocks[key[:8]][key[8:]]=coefficient
    blocks={key:ctx2.from_dict(terms) for key,terms in blocks.items()}
    constant=blocks.pop((0,)*8)
    assert constant==prod(ctx2.gens())**2
    return constant,blocks


def node_residue(target,shifts,axis):
    n=sum(target)
    degree=2*(n-shifts)
    if degree<0:return Fraction(0)
    constant,blocks=node_parts(axis)
    ctx=constant.context()
    exponent=fmpq(2*shifts-1,2)
    parts=[(beta,sum(beta),poly*constant**(sum(beta)-1))
           for beta,poly in blocks.items() if all(a<=b for a,b in zip(beta,target))]
    @cache
    def jet(alpha):
        order=sum(alpha)
        if not order:return ctx.constant(1)
        value=ctx.constant(0)
        for beta,size,poly in parts:
            if all(a<=b for a,b in zip(beta,alpha)):
                previous=tuple(b-a for a,b in zip(beta,alpha))
                value+=((exponent+1)*size-order)*poly*jet(previous)
        return value/order
    result=fmpq(0)
    for (i,j),value in jet(target).terms():
        a,b=degree-int(i),degree-int(j)
        if min(a,b)>=0:result+=value*(-1)**(a+b)*comb(a+b,a)
    jet.cache_clear()
    return Fraction(str(result))


def vertex_projection(powers,shifts):
    if min(powers)<1:raise ValueError('Positive propagator powers are required.')
    exchange=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[5],g[8]}=={6,7})
    target=tuple(a-1 for a in powers)
    images=[target,tuple(target[i-1] for i in exchange[1:])]
    total=sum(node_residue(row,shifts,axis) for row in images for axis in (4,7))
    top=projection(powers,shifts)
    scale=Fraction(8**shifts,prod((2*j+1)*(2*j+2)*(2*j+3) for j in range(shifts)))
    return 2*top-Fraction(4)**(1-2*shifts)*scale*total/8


@cache
def infinity_polynomials(axis):
    """Homogenize in the cut momentum squared and diagonalize the leading conic."""
    ctx=fmpq_mpoly_ctx.get(tuple(f'y{i}' for i in range(9))+('w','x'))
    full=ctx.from_dict({(4-sum(e),)+tuple(map(int,e))+(0,):v
                        for e,v in endpoint_gram(0).terms()})
    variables=list(ctx.gens())
    w,x=variables[-2:]
    l,r=(2,4) if axis==4 else (1,3)
    arguments=variables[:]
    arguments[axis+1]-=x
    arguments[-2]=w-variables[l]-variables[r]+variables[6]+variables[7]
    translated=full.compose(*arguments).subs({9:1})
    c=fmpq_mpoly_ctx.get(tuple(f'y{i}' for i in range(9)))
    pieces=[c.from_dict({tuple(map(int,e[:9])):v for e,v in translated.terms() if e[-1]==k})
            for k in range(3)]
    leading=pieces[2]
    assert leading==1-4*c.gen(l)*c.gen(r)
    return (l,r),[[(tuple(map(int,e)),int(v)) for e,v in poly.terms()]
                  for poly in (pieces[1],leading*pieces[0])]


def infinity_residue(target,shifts,axis):
    """Normalized coefficient at x=infinity after the remaining w residue."""
    from math import factorial
    order=2*shifts-1
    if order<0 or min(target)<0:return Fraction(0)
    pair,polynomials=infinity_polynomials(axis)
    operators=[[(ex,c) for ex,c in p if all(k<=a for k,a in zip(ex,target))]
               for p in polynomials]
    @cache
    def coefficient(k,alpha):
        if k<0 or sum(alpha)>3*k:return 0
        if not k:return int(not any(alpha))
        total=0
        for i,factor in ((1,2*shifts+1-2*k),(2,4*(k-1)*(2*shifts+1-k))):
            if not factor:continue
            for ex,weight in operators[i-1]:
                if all(e<=a for e,a in zip(ex,alpha)):
                    total+=factor*weight*coefficient(k-i,tuple(a-e for a,e in zip(alpha,ex)))
        return total
    result=0
    binomial=1
    for j in range(min(target[k] for k in pair)+1):
        alpha=tuple(a-int(k in pair)*j for k,a in enumerate(target))
        result+=binomial*coefficient(order,alpha)
        binomial=binomial*2*(order+2*j)//(j+1)
    coefficient.cache_clear()
    return Fraction(result,2**order*factorial(order))


def infinity_projection(powers,shifts):
    target=(4*shifts-1-sum(a-1 for a in powers),)+tuple(a-1 for a in powers)
    if min(target)<0:return Fraction(0)
    exchange=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[5],g[8]}=={6,7})
    images=[target,tuple(target[i] for i in exchange)]
    total=sum(infinity_residue(row,shifts,axis) for row in images for axis in (4,7))
    scale=Fraction(8**shifts,prod((2*j+1)*(2*j+2)*(2*j+3) for j in range(shifts)))
    return -Fraction(4)**(1-2*shifts)*scale*total/8


@cache
def zero_endpoint_cuts():
    """Verify the complete affine Gram transformations in homogeneous variables."""
    ctx=fmpq_mpoly_ctx.get(tuple(f'y{i}' for i in range(9))+('w',))
    variables=list(ctx.gens())
    w=variables[-1]
    gram=ctx.from_dict({(4-sum(e),)+tuple(map(int,e)):v for e,v in endpoint_gram(0).terms()})
    cubic=ctx.from_dict({tuple(map(int,e[:9]))+(0,):v for e,v in gram.terms() if e[-1]==3})
    result=[]
    for axis in range(9):
        valid=[]
        for g in AUTOMORPHISMS:
            if g[0]!=axis:continue
            changed=[variables[i] for i in g]
            permuted_cubic=cubic.compose(*changed,w)
            for sign in (1,-1):
                translation=(sign*cubic-permuted_cubic)/4
                if translation.derivative(axis):continue
                if gram==gram.compose(*changed,sign*w+translation):valid.append((sign,g))
        assert len(valid)==4
        sign,g=valid[0]
        assert all(v==sign and {q[6],q[7]}=={g[6],g[7]} for v,q in valid)
        result.append((sign,g))
    assert [sign for sign,g in result]==[1,-1,-1,-1,-1,1,1,1,1]
    return tuple(result)


@cache
def recut_permutations():
    """Each negative disk pole is one endpoint on a different edge."""
    exchange=next(g for g in AUTOMORPHISMS if g[0]==0 and {g[5],g[8]}=={6,7})
    assignments=[]
    for axis,pair in ((5,{2,4}),(8,{1,3})):
        sign,g=zero_endpoint_cuts()[axis]
        assert sign==1 and {g[6],g[7]}==pair
        assignments.extend([g,tuple(exchange[i] for i in g)])
    return tuple(assignments)


def recut_projection(powers,shifts):
    from quartic_nonplanar_cut import homogeneous_projection
    cut_power=4*shifts-sum(a-1 for a in powers)
    if cut_power<1:return Fraction(0)
    full=(cut_power,)+tuple(powers)
    return sum(homogeneous_projection(tuple(full[i] for i in g[1:]),shifts,(0,))
               for g in recut_permutations())


def disjoint_projection(powers,shifts):
    from quartic_nonplanar_cut import homogeneous_projection
    cut_power=4*shifts-sum(a-1 for a in powers)
    if cut_power<1:return Fraction(0)
    full=(cut_power,)+tuple(powers)
    values=[]
    differences=[]
    for sign,g in zero_endpoint_cuts():
        p=tuple(full[i] for i in g[1:])
        e=[homogeneous_projection(p,shifts,(k,)) for k in (0,1)]
        differences.append(sign*(e[0]-e[1]))
        values.append(sum(e))
    assert len(set(differences))==1
    return sum(values[5:])/2


def controls():
    import json
    from pathlib import Path
    from quartic_planar_double_pole import physical_poles,two_poles,disk_poles
    rows=json.loads(Path(__file__).with_name('quartic_planar_terminal.json').read_text())
    for row in rows:
        p=row['uncut']
        if min(p)>0:
            expected=Fraction(row['top,planar_double,elementary_double'][1])
            actual=vertex_projection(p,1)
            assert actual==infinity_projection(p,1)==recut_projection(p,1)==disjoint_projection(p,1)==expected,(p,actual,expected)
    for level in range(5):
        for k in range(4*level+5):
            pure=Fraction(16**level*comb(k-2*level-1,2*level),4) if k>4*level else 0
            assert two_poles(k,0,-level)==(0,pure)
            for i in range(2*level+5):
                value=(Fraction(16**level,2)*(-1)**(i+k)*comb(i,2*level)*comb(i-2*level,2*level-k)
                       if 0<=2*level-k<=i-2*level else 0)
                assert disk_poles(i,-1,k,0,-level)[0]==value
    probes=((2,1,2,3,1,2,1,1),(3,1,2,2,3,2,1,1),(1,2,1,2,1,2,1,2),(1,1,1,1,1,1,1,1))
    for m in range(3):
        for p in probes:
            if sum(p)-8<m:continue
            actual=vertex_projection(p,m)
            assert infinity_projection(p,m)==recut_projection(p,m)==disjoint_projection(p,m)==actual,(p,m,actual,infinity_projection(p,m))
            assert physical_poles(p,m)==(0,actual),(p,m,actual,physical_poles(p,m))
    print('Four vertex residues reproduce every positive terminal control and the independent disk probes.')


if __name__=='__main__':
    controls()
