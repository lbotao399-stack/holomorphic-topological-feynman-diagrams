"""Endpoint residues of the nonplanar three-loop master at even dimension."""
from argparse import ArgumentParser
from functools import cache
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
from itertools import islice
from gzip import open as open_gzip
from math import comb,factorial,prod
from pathlib import Path
from tempfile import TemporaryDirectory

from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx
import numpy as np
import sympy as sp

from quartic_projective_ibp import AUTOMORPHISMS,K33
from quartic_prism_shift import G
from four_loop import form,read_vectors,tagged
from quartic_bipartite_shift import minor
from shift import gram


@cache
def baikov():
    cycle=sp.Matrix(G(6,K33).C)
    routing=cycle[:,[0,6,7,8]].inv()*cycle
    slots=list(combinations_with_replacement(range(4),2))
    metric=sp.Matrix([[routing[i,e]*routing[j,e]*(1+(i!=j)) for i,j in slots] for e in range(9)]
                     +[[int((i,j)==(1,2)) for i,j in slots]])
    coefficients=metric.inv()
    ctx=fmpq_mpoly_ctx.get(tuple(f'y{i}' for i in range(1,9))+('x',))
    variables=ctx.gens()
    sources=(1,*variables[:8],-variables[8]/2)
    entries=[sum(fmpq(str(c))*v for c,v in zip(coefficients.row(i),sources)) for i in range(10)]
    matrix=[[entries[slots.index(tuple(sorted((i,j))))] for j in range(4)] for i in range(4)]
    determinant=minor(matrix,range(4))
    coefficients={}
    for ex,c in determinant.terms():
        edge_ex=tuple(map(int,ex[:8]))
        coefficients.setdefault(edge_ex,fmpq_poly())[int(ex[8])]=c
    constant=coefficients.pop((0,)*8)
    assert constant==fmpq_poly([0,1,-1])**2/16
    return constant,coefficients


def projection(uncut,shifts):
    if min(uncut)<=0:return F(0)
    target=tuple(int(v)-1 for v in uncut)
    if sum(target)<shifts:return F(0)
    constant,parts=baikov()
    exponent=fmpq(2*shifts-1,2)
    # Clear precisely one power of the constant term for each derivative.
    entries=[(key,sum(key),coeff*constant**(sum(key)-1)) for key,coeff in parts.items()
             if all(k<=a for k,a in zip(key,target))]
    @cache
    def jet(alpha):
        n=sum(alpha)
        if not n:return fmpq_poly([1])
        value=fmpq_poly()
        for beta,size,poly in entries:
            if all(k<=a for k,a in zip(beta,alpha)):
                value+=((exponent+1)*size-n)*poly*jet(tuple(a-k for a,k in zip(alpha,beta)))
        return value/n
    p=jet(target)
    degree=2*sum(target)+1-2*shifts
    assert degree>0
    def residue(poly):
        return sum(poly[j]*comb(2*degree-2-j,degree-1) for j in range(min(poly.degree()+1,degree)))
    total=4**degree*(residue(p)+residue(p(fmpq_poly([1,-1]))))/8
    # The unit-index scalar master has twice the endpoint residue 4.
    prefactor=F(8**shifts,prod((2*j+3)*(2*j+2)*(2*j+1) for j in range(shifts)))
    return prefactor*F(str(total))


@cache
def endpoint_gram(endpoint):
    constant,parts=baikov()
    ctx=fmpq_mpoly_ctx.get(tuple(f'u{i}' for i in range(8))+('w',))
    variables=ctx.gens()
    gram=ctx.from_dict({e+(j,):16*v for e,p in {**parts,(0,)*8:constant}.items()
                        for j,v in enumerate(p) if v})
    if endpoint:
        # Translate the auxiliary coordinate to put the two-period tangent
        # cone into diagonal form. Formal residues are translation invariant.
        gram=gram.compose(*variables[:8],1-variables[-1]-sum(variables[:4])+sum(variables[4:8]))
    return gram


@cache
def homogeneous_endpoint(endpoint):
    gram=endpoint_gram(endpoint)
    ctx=gram.context()
    variables=ctx.gens()
    pieces=[ctx.from_dict({tuple(map(int,e)):v for e,v in gram.terms() if sum(e)==k})
            for k in (2,3,4)]
    assert sum(pieces)==gram
    indices=(4,7) if endpoint else (5,6)
    a=4*prod(variables[j] for j in indices)
    assert pieces[0]==variables[-1]**2-a
    polynomials=[pieces[1],pieces[2]*(1-a)]
    return indices,[[(tuple(map(int,e[:-1])),int(v)) for e,v in p.subs({8:1}).terms()] for p in polynomials]


def cut_tests():
    normals=[]
    for endpoint in (0,1):
        gram=endpoint_gram(endpoint)
        constraints=[(tuple(map(int,e[:-1])),2-int(e[-1])) for e,v in gram.terms() if e[-1]<2]
        weights=sorted(product(range(3),repeat=8),key=sum)
        accepted=[]
        for w in weights:
            if any(all(a<=b for a,b in zip(prev,w)) for prev in accepted):continue
            if all(sum(a*b for a,b in zip(w,e))>=k for e,k in constraints):accepted.append(w)
        normals.append(np.array(accepted,dtype=np.int16))
    assert list(map(len,normals))==[4,4]
    # Choose from all cuts separately for each orbit, retaining an ordinary
    # finite period in every case. Integer scores refer to actual exponents.
    tests=[]
    for group in AUTOMORPHISMS:
        extended=[]
        for matrix in normals:
            rows=[]
            for w in matrix:
                row=[0]*9
                for i,k in enumerate(w,1):row[group[i]]=k
                rows.append(tuple(row))
            extended.append(tuple(sorted(rows)))
        tests.append(tuple(sorted(extended)))
    tests=sorted(set(tests))
    assert len(tests)==9
    return tests


def cut_selection(emit=None):
    tests=cut_tests()
    path=Path(__file__).with_name('quartic_k33.weights.gz')
    counts=[0,0,0]
    transformations=np.array(AUTOMORPHISMS)
    def batch(exponents,raw):
        exponents=np.array(exponents,dtype=np.int16)
        good=np.zeros(len(exponents),np.bool_)
        for left,right in tests:
            good|=((exponents@np.array(left,np.int16).T<20).any(axis=1)
                   &(exponents@np.array(right,np.int16).T<20).any(axis=1))
        counts[0]+=len(exponents)
        counts[1]+=int(good.sum())
        counts[2]+=int((~good).sum())
        if emit:
            for idx in np.nonzero(~good)[0]:
                images=exponents[idx][transformations]
                highest=images[:,0].max()
                y=min(map(tuple,images[images[:,0]==highest,1:]))
                emit(tuple(map(int,y)),raw[idx].split()[-2:])
    with open_gzip(path,'rt') as stream:
        while lines:=list(islice(stream,4096)):
            batch([tuple(int(v)-1 for v in line.split()[:9]) for line in lines],lines)
    assert counts==[465564,395810,69754]
    print('physical orbits with certified zero-master choices and remaining orbits:',counts,flush=True)
    return counts


def homogeneous_projection(uncut,shifts,selected=(0,1)):
    target=tuple(i-1 for i in uncut)
    if min(target)<0:return F(0)
    r=sum(target)-2*shifts
    if r<0:return F(0)
    shape=tuple(v+1 for v in target)
    endpoints=[]
    bound=0
    for endpoint in selected:
        pair,polynomials=homogeneous_endpoint(endpoint)
        def slices(ex):
            dest=tuple(slice(k,None) for k in ex)
            source=tuple(slice(None,-k or None) for k in ex)
            return dest,source
        operators=[[(c,*slices(ex)) for ex,c in rows if all(k<=a for k,a in zip(ex,target))]
                   for rows in polynomials]
        weight=F(1)
        residues=[]
        for j in range(min(target[k] for k in pair)+1):
            key=tuple(v-int(k in pair)*j for k,v in enumerate(target))
            assert weight.denominator==1
            residues.append((key,int(weight)))
            weight*=F(-2*(2*shifts-1-2*r-2*j),j+1)
        norm0,norm1=0,1
        l1=[sum(abs(c) for c,_,_ in rows) for rows in operators]
        for k in range(1,r+1):
            norm0,norm1=norm1,abs(2*shifts+1-2*k)*l1[0]*norm1+abs(4*(k-1)*(2*shifts+1-k))*l1[1]*norm0
        bound+=norm1*sum(abs(c) for _,c in residues)
        endpoints.append((operators,residues))
    def one_prime(prime):
        value=0
        for operators,residues in endpoints:
            previous=np.zeros(shape,dtype=np.int64)
            current=np.zeros(shape,dtype=np.int64)
            current[(0,)*8]=1
            for k in range(1,r+1):
                result=np.zeros(shape,dtype=np.int64)
                for scale,rows,source in ((2*shifts+1-2*k,operators[0],current),
                                          (4*(k-1)*(2*shifts+1-k),operators[1],previous)):
                    if scale:
                        for coeff,dest,src in rows:
                            result[dest]+=(scale*coeff)*source[src]
                previous,current=current,result%prime
            value+=sum(c*int(current[key]) for key,c in residues)
        return value%prime
    modulus,value=1,0
    prime=2**31
    # The largest signed accumulation at one cell is bounded before every
    # modular operation, independently of the requested exponents.
    max_multiplier=max([1]+[abs(2*shifts+1-2*k)*sum(abs(c) for c,_,_ in rows[0])+
                   abs(4*(k-1)*(2*shifts+1-k))*sum(abs(c) for c,_,_ in rows[1])
                   for rows,_ in endpoints for k in range(1,r+1)])
    prime=min(prime,(2**63-1)//max_multiplier)
    while modulus<=2*bound:
        prime=int(sp.prevprime(prime))
        residue=one_prime(prime)
        value+=modulus*((residue-value)*pow(modulus,-1,prime)%prime)
        modulus*=prime
    if 2*value>modulus:value-=modulus
    assert abs(value)<=bound
    prefactor=F(8**shifts,prod((2*j+3)*(2*j+2)*(2*j+1) for j in range(shifts)))
    return prefactor*F(value,2**r*factorial(r)*8)*F(4)**(1-2*shifts)


def physical_control(folder,executable):
    powers=(9,6,5,6,4,5,4,5,4)
    shifts=10
    ep=sp.Symbol('ep')
    basis='*'.join(f'1/(p{j}.p{j})' for j in range(1,9))
    term='*'.join(f'1/(p{j}.p{j})^{a}' for j,a in enumerate(powers[1:],1))
    rows=read_vectors(form(folder,executable,'maxcut_D24',tagged([(0,basis),(1,term)]),edges=K33))
    def rational(pair):
        return sp.cancel(sp.sympify(pair[0].replace('^','**'))/sp.sympify(pair[1].replace('^','**')))
    ratio=sp.cancel((rational(rows[1]['BasicNOIntegral'])/rational(rows[0]['BasicNOIntegral'])).subs(ep,ep-shifts))
    # The maximal contour of the raised unit graph is the beta function
    # B(d-4+2m,d-4+2m); divide by the original B(d-4,d-4).
    d=4-2*ep
    x=d-4
    shift=sp.rf(x,2*shifts)**2/sp.rf(2*x,4*shifts)/16**shifts
    shift*=8**shifts/prod((d+2*j-1)*(d+2*j-2)*(d+2*j-3) for j in range(shifts))
    top=sp.cancel(ratio*shift).subs(ep,0)
    expected=homogeneous_projection(powers[1:],shifts)
    assert top==sp.Rational(expected)==0
    normalization=F(prod(factorial(a-1) for a in powers),factorial(11))
    print('physical D24 top master',expected,'its projective zeta(5) term',20*normalization*expected,flush=True)
    assert homogeneous_projection((1,3,3,3,9,9,9,10),10)==F(-1004423629,95274484236288000)
    return expected


def compare(folder,executable):
    examples=[(1,)*8,(2,2,1,1,1,1,1,1),(1,2,1,2,1,2,1,1),(2,1,3,1,2,1,1,2),
              (1,1,1,1,1,2,2,1),(1,2,1,1,3,1,2,2),
              (0,1,1,1,1,1,1,1)]
    ep=sp.Symbol('ep')
    def expr(p):
        return '*'.join(f'1/(p{j}.p{j})^{v}' for j,v in enumerate(p,1) if v)
    unit=read_vectors(form(folder,executable,'maxcut_unit',tagged([(0,expr(examples[0]))]),edges=K33))
    def rational(pair):
        return sp.cancel(sp.sympify(pair[0].replace('^','**'))/sp.sympify(pair[1].replace('^','**')))
    unit=rational(unit[0]['BasicNOIntegral'])
    counts=[]
    for shifts in (0,1,2):
        # Reduce the Gram insertion at exactly the same four-dimensional
        # master normalization; lower sectors vanish on this maximal cut.
        terms=[(i,expr(p)+(f'*{gram(K33)}^{shifts}' if shifts else '')) for i,p in enumerate(examples)]
        rows=read_vectors(form(folder,executable,f'maxcut_{shifts}',tagged(terms),edges=K33))
        expected=[]
        for i,p in enumerate(examples):
            row=rows.get(i,{})
            coefficient=rational(row.get('BasicNOIntegral',('0','1')))/unit
            d=4-2*ep
            normalized=sp.cancel(coefficient*8**shifts/prod((d+2*j-1)*(d+2*j-2)*(d+2*j-3)
                                                             for j in range(shifts)))
            assert normalized.as_numer_denom()[1].subs(ep,0)!=0
            value=projection(p,shifts)
            assert sp.Rational(value)==normalized.subs(ep,0)
            assert homogeneous_projection(p,shifts)==value
            expected.append(value)
        counts.append(tuple(expected))
    assert counts[1][4]==F(-1,12) and counts[2][5]==F(1,240)
    print('maximal cut against FORCER at pole orders 2, 3, 4:',counts,flush=True)
    return counts


if __name__=='__main__':
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('form',type=Path)
    parser.add_argument('forcer',type=Path)
    parser.add_argument('--physical',action='store_true')
    parser.add_argument('--cut-selection',action='store_true')
    args=parser.parse_args()
    with TemporaryDirectory(prefix='bcov-nonplanar-cut-') as temporary:
        folder=Path(temporary)
        for item in ('forcer','forcer.h'):
            (folder/item).symlink_to((args.forcer/item).resolve(strict=True))
        if args.cut_selection:cut_selection()
        elif args.physical:physical_control(folder,args.form.resolve(strict=True))
        else:compare(folder,args.form.resolve(strict=True))
