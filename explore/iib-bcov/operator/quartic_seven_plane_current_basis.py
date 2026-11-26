"""Twenty exact normal-source samples for seven-planar four-loop current minors."""
from fractions import Fraction as F
from functools import cache, reduce
from itertools import combinations, combinations_with_replacement, permutations
from math import gcd, lcm

from flint import fmpq_mat, fmpq_mpoly_ctx

from quartic_seven_plane_double import NORMAL, polynomials
from quartic_seven_plane_gale import canonical_sources, cubic_nodes, gale, monomials


@cache
def label_interpolation(normal,slopes):
    """Coordinates of every ordered cofactor in a basis of signed cofactor samples.

    The 3 x 4 normal matrix P uses four independent conserved source slots.
    Adding rows of P to the lift leaves its cofactor invariant. Consequently
    det(lift;P) divides the cofactor: its four linear coefficients are
    proportional to the signed 3-minors, which have gcd one. The quotient
    transforms by det(g)^3 under GL(3), hence is a homogeneous cubic in the
    four minors. There are no Pluecker relations for Gr(3,4).
    """
    orders=list(permutations(range(5)))
    rows=[cubic_row(tuple(normal[i] for i in order),tuple(slopes[i] for i in order)) for order in orders]
    matrix=fmpq_mat(rows).transpose()
    echelon,rank=matrix.rref()
    indices=[next(j for j in range(len(orders)) if echelon[i,j]) for i in range(rank)]
    basis=tuple(orders[j] for j in indices)
    coords={o:tuple(echelon[i,j] for i in range(rank)) for j,o in enumerate(orders)}
    if rank:
        selected=fmpq_mat([[matrix[i,j] for j in indices] for i in range(20)])
        assert selected*fmpq_mat([list(coords[o]) for o in orders]).transpose()==matrix
    else:
        assert not any(map(any,rows))
    return basis,coords


def cubic_row(normal,slopes):
    minors=gale(normal)
    det4=sum(x*y for x,y in zip(slopes,minors))
    return [det4*v for v in monomials(minors)]


@cache
def universal_samples():
    """Pure cubes, both polarizations of two coordinates, and three-coordinate sums."""
    samples=[]
    for node in cubic_nodes():
        source,_=canonical_sources(node)
        normal=tuple(tuple(row[1:]) for row in source[6:])
        slopes=tuple(row[0] for row in source[6:])
        samples.append((normal,slopes))
    matrix=fmpq_mat([cubic_row(normal,slopes) for normal,slopes in samples])
    assert matrix.nrows()==matrix.ncols()==20 and matrix.det()
    return samples,matrix.inv()


def expand(basis,poly):
    """Append a pivot when necessary; return the exact coordinates of this polynomial."""
    coords=[]
    for j,(pivot,p) in enumerate(basis):
        c=poly[pivot]
        if c:
            poly-=c*p
            coords.append((j,c))
    if poly:
        pivot,lead=next(iter(poly.terms()))
        coords.append((len(basis),lead))
        basis.append((pivot,poly/lead))
    return coords


def integer_poly(poly,ctx):
    """Primitive integral factor and its rational content for a nonzero polynomial."""
    den=lcm(*(int(c.denominator) for _,c in poly.terms()))
    terms=[(tuple(m),int(c*den)) for m,c in poly.terms()]
    content=reduce(gcd,(v for _,v in terms))
    return ctx.from_dict({m:v//content for m,v in terms}),F(content,den)


def polarized_currents(data,orders,**kinematics):
    normal=tuple(kinematics.get('normal',NORMAL))
    slopes=tuple(kinematics.get('slopes',(1,-1,0,0,0)))
    pair=kinematics.get('pair',(0,1))
    polarization=kinematics.get('polarization',(0,0,1))
    if (len(normal)!=5 or any(len(n)!=3 for n in normal) or len(slopes)!=5 or len(polarization)!=3
            or any(map(sum,zip(*normal))) or sum(slopes)):
        raise ValueError('Expected five conserved three-component normal momenta and five conserved slopes.')
    if len(pair)!=2 or pair[0]==pair[1] or not all(0<=i<5 for i in pair):
        raise ValueError('Expected two distinct observed normal slots.')
    charges=[sum(x*y for x,y in zip(n,polarization)) for n in normal]
    orders=list(orders)
    if not orders:
        return {},{}
    if any(sorted(order)!=list(range(5)) for order in orders):
        raise ValueError('Each normal order must be a permutation of the five labels.')
    samples,inverse=universal_samples()
    evaluation=fmpq_mat([cubic_row(tuple(normal[j] for j in order),tuple(slopes[j] for j in order))
                         for order in orders])
    coefficients=evaluation*inverse
    ctx=data[1]
    qq=fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
    bases=[[] for _ in range(4)]
    images=[]
    for i,(source,lift) in enumerate(samples):
        if not any(coefficients[j,i] for j in range(len(orders))):
            images.append([[] for _ in range(4)])
            continue
        ds,_=polynomials(data,tuple(range(5)),normal=source,slopes=lift)
        images.append([expand(basis,qq.from_dict(d.to_dict())) for basis,d in zip(bases,ds)])

    labels=[(r,j) for r,bank in enumerate(bases) for j in range(len(bank))]
    if not labels:
        return {},{}
    input_rows=[]
    for index,order in enumerate(orders):
        rows=[{} for _ in range(4)]
        for k,img in enumerate(images):
            a=coefficients[index,k]
            if a:
                for row,terms in zip(rows,img):
                    for j,c in terms:
                        row[j]=row.get(j,0)+a*c
        if order in (tuple(reversed(range(5))),(3,2,1,4,0)):
            direct,_=polynomials(data,order,original=True,**kinematics)
            for r,row in enumerate(rows):
                assert qq.from_dict(direct[r].to_dict())==sum(c*bases[r][j][1] for j,c in row.items())
        input_rows.append([rows[r].get(j,0) for r,j in labels])

    # V is linear in the ten pair resistances. Project the external orders into
    # their symmetric coefficient matrices before multiplying forest polynomials.
    matrix=fmpq_mat(input_rows)
    transposed=matrix.transpose()
    voltages={}
    for left,right in combinations(range(5),2):
        weighted=[]
        for order,row in zip(orders,input_rows):
            a,b=order[left],order[right]
            sign_a=int(a==pair[1])-int(a==pair[0])
            sign_b=int(b==pair[1])-int(b==pair[0])
            w=sign_a*charges[b]+sign_b*charges[a]
            weighted.append([w*c for c in row])
        gram=transposed*fmpq_mat(weighted)
        resistance=qq.from_dict(data[4][6+left][6+right].to_dict())
        for i,j in combinations_with_replacement(range(len(labels)),2):
            c=gram[i,j]
            if c:
                key=(*labels[i],*labels[j])
                factor=c if i!=j else c/2
                voltages[key]=voltages.get(key,0)+resistance*factor
    basis={(r,j):integer_poly(p,ctx) for r,bank in enumerate(bases) for j,(_,p) in enumerate(bank)}
    return basis,voltages
