"""Quadratic lift on the two prism orbits with path occupancies (2,2,1).

The source-independent current factors are divided out before the label sum.
The angular push reads each fourteen-variable FLINT term separately.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, lcm, prod
from operator import itemgetter

from flint import fmpq, fmpq_mat, fmpq_mpoly_ctx
from quartic_seven_plane_current_basis import cubic_row, expand, integer_poly, universal_samples
from quartic_seven_plane_prism import CTX, NORMAL, STARS, T, cleared_stars, graph_data, polynomials

MARKS=((3,3,5,5,7),(3,3,5,7,7))
PAIRS=tuple(combinations(range(5),2))
QCTX=fmpq_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
ANGULAR_SCALE=factorial(25)**3


def star_numerator(poly,stars,pole):
    """Substitute each star and multiply by the prescribed common S^pole."""
    ctx=poly.context()
    variables=ctx.gens()
    for indices in stars:
        if not poly:
            break
        s=sum(variables[i] for i in indices)
        slices=defaultdict(dict)
        for ex,c in poly.terms():
            m=list(map(int,ex))
            degree=sum(m[i] for i in indices)
            for i in indices:
                m[i]=degree-m[i]
            slices[degree][tuple(m)]=int(c)
        top=max(pole,max(slices))
        poly=sum(ctx.from_dict(row)*s**(top-d) for d,row in slices.items())
        if top>pole:
            poly,remainder=divmod(poly,s**(top-pole))
            assert not remainder
    return poly


def current_basis(currents,ctx):
    common=next(p for p in currents if p)
    for p in currents:
        if p:
            common=common.gcd(p)
    echelon=[]
    coordinates=[]
    for p in currents:
        quotient,rem=divmod(p,common)
        assert not rem
        coordinates.append(dict(expand(echelon,QCTX.from_dict(quotient.to_dict()))))
    factors=[integer_poly(row,ctx) for _,row in echelon]
    coords=fmpq_mat([[coords.get(i,0)*fmpq(scale.numerator,scale.denominator)
                      for i,(_,scale) in enumerate(factors)] for coords in coordinates])
    return common,[p for p,_ in factors],coords


@cache
def geometry(marks):
    if marks not in MARKS:
        raise ValueError('Expected two marked pairs and a singleton on the surviving prism paths.')
    data=graph_data(marks)
    ctx,adj,paths=data[1],data[3],data[5]
    stars=tuple(tuple(paths[i][0] for i in edges) for edges in STARS)
    samples,inverse=universal_samples()
    groups=[]
    for normal,slopes in samples:
        currents,_=polynomials(data,tuple(range(5)),normal=normal,slopes=slopes)
        groups.append([star_numerator(d,stars,3) for d in currents])
    families=[]
    assignment=[]
    for cycle,row in enumerate(zip(*groups)):
        if not any(row):
            continue
        for i,family in enumerate(families):
            if row==family:
                assignment.append((i,cycle))
                break
        else:
            assignment.append((len(families),cycle))
            families.append(row)
    pieces=[current_basis(row,ctx) for row in families]

    # Sum the adjugate entries for coincident current components.
    primaries=[path[0] for path in paths]
    other=set(range(14))-set(primaries)
    reduced={}
    for i,j in combinations_with_replacement(range(len(pieces)),2):
        a=sum(adj[r][s] for h,r in assignment for k,s in assignment if (h,k)==(i,j))
        poly=CTX.from_dict({tuple(m[x] for x in primaries):int(c) for m,c in a.terms()
                            if not any(m[x] for x in other)})
        reduced[i,j]=cleared_stars({0:poly},1).get(0,CTX.constant(0))
    return data,stars,pieces,inverse,reduced


def external_coefficients(data,stars,pieces,inverse,orders,normal,slopes,table):
    values=fmpq_mat([cubic_row(tuple(normal[j] for j in order),
                               tuple(slopes[j] for j in order)) for order in orders])
    evaluation=values*inverse
    currents=[evaluation*coords for _,_,coords in pieces]
    # The third component of the normal source tests the voltage; the first
    # two external labels select its endpoints. Each distance is a two-forest
    # numerator. Reduce the ten pair coefficients once per external frame.
    if table is None:
        table=fmpq_mat([[-fmpq((int(o[i]==0)-int(o[i]==1))*normal[o[j]][2]
                               +(int(o[j]==0)-int(o[j]==1))*normal[o[i]][2],2)
                        for i,j in PAIRS] for o in orders])
    assert table.nrows()==len(orders) and table.ncols()==len(PAIRS)
    rows,rank=table.rref()
    voltages=[]
    for k in range(rank):
        pivot=next(i for i in range(len(PAIRS)) if rows[k,i])
        poly=sum(rows[k,l]*data[4][6+i][6+j] for l,(i,j) in enumerate(PAIRS))
        if not poly:
            continue
        p,content=integer_poly(poly,data[1])
        scale=fmpq(content.numerator,content.denominator)
        voltages.append((star_numerator(p,stars,1),[table[i,pivot]*scale for i in range(len(orders))]))
    return currents,voltages


@cache
def simplex_weights(exponents,second_end):
    # The first marked path always has three segments. The remaining two
    # have lengths (3,2) or (2,3) in the two allowed prism geometries.
    powers=(sum(exponents[:3])+2,sum(exponents[3:second_end])+second_end-4,
            sum(exponents[second_end:])+len(exponents)-second_end-1)
    value,rem=divmod(ANGULAR_SCALE*prod(factorial(e) for e in exponents),
                     prod(factorial(e) for e in powers))
    assert not rem
    return powers,value


def position_integral(raw,paths):
    # Keep the star indices and integrate the ordered subsegments at fixed
    # total time on each of the three marked paths.
    rows=defaultdict(int)
    nontrivial=(paths[3],paths[5],paths[7])
    marked=itemgetter(*(i for path in nontrivial for i in path))
    second_end=len(paths[3])+len(paths[5])
    first=tuple(p[0] for p in paths)
    for n in range(len(raw)):
        m=raw.monomial(n)
        (a,b,c),weight=simplex_weights(marked(m),second_end)
        key=(m[first[0]],m[first[1]],m[first[2]],a,m[first[4]],b,m[first[6]],c,m[first[8]])
        rows[key]+=int(raw.coefficient(n))*weight
    return CTX.from_dict(dict(rows))/ANGULAR_SCALE


def degree_slices(poly):
    # Two subsegment exponents define an additive grading for multiplication.
    ctx=poly.context()
    buckets=defaultdict(dict)
    for m,c in poly.terms():
        buckets[int(m[4])+2*int(m[11])][m]=int(c)
    return {degree:ctx.from_dict(terms) for degree,terms in buckets.items()}


def integrated_square(p,q,kernel,blocks,denominator,paths):
    # Accumulate one homogeneous slice at a time. Within that slice multiply
    # by a voltage component only after summing all its current products.
    ctx=kernel.context()
    left=[degree_slices(poly) for poly in p]
    by_degree=defaultdict(lambda:defaultdict(list))
    voltage_slices=[]
    for channel,(lines,voltage) in enumerate(blocks):
        sliced_voltage=degree_slices(voltage)
        voltage_slices.append(sliced_voltage)
        for L,line in zip(left,lines):
            if not line:
                continue
            weighted=sum((q[j]*int(denominator*c) for j,c in line),ctx.constant(0))
            for b,right in degree_slices(weighted).items():
                for v in sliced_voltage:
                    for a,factor in L.items():
                        by_degree[a+b+v][channel,v].append((factor,right))
    total=CTX.constant(0)
    for groups in by_degree.values():
        raw=ctx.constant(0)
        for (channel,v),terms in groups.items():
            current=sum((a*b for a,b in terms),ctx.constant(0))
            raw+=current*voltage_slices[channel][v]
        if raw:
            total+=position_integral(raw*kernel,paths)
    return total/denominator


def project(marks,orders,normal=NORMAL,*,slopes=(1,-1,0,0,0),table=None):
    data,stars,pieces,inverse,adjugate=geometry(marks)
    orders=list(orders)
    if not orders:
        return CTX.constant(0)
    coefficients,voltages=external_coefficients(data,stars,pieces,inverse,orders,normal,slopes,table)
    result=CTX.constant(0)
    for (r,s),adj in adjugate.items():
        if not adj:
            continue
        blocks=[]
        denominator=1
        for voltage,weight in voltages:
            left,right=coefficients[r],coefficients[s]
            matrix=left.transpose()*fmpq_mat([[weight[i]*right[i,j] for j in range(right.ncols())]
                                              for i in range(len(orders))])
            lines=[]
            for i in range(matrix.nrows()):
                line=[]
                for j in range(i if r==s else 0,matrix.ncols()):
                    if matrix[i,j]:
                        entry=matrix[i,j]*(1+((r,i)!=(s,j)))
                        denominator=lcm(denominator,int(entry.denominator))
                        line.append((j,entry))
                lines.append(line)
            blocks.append((lines,voltage))
        curr,p,_=pieces[r]
        other,q,_=pieces[s]
        result+=integrated_square(p,q,curr*other,blocks,denominator,data[5])*adj
    jac=prod(prod(T[i] for i in edges)*sum(T[i] for i in edges)**2 for edges in STARS)
    return jac*result


if __name__=='__main__':
    from argparse import ArgumentParser
    from quartic_seven_plane_prism import cut_weights,cuts,frame
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('orbit',type=int,choices=(32,42))
    parser.add_argument('--frame',type=int,choices=(2,3))
    args=parser.parse_args()
    marks=MARKS[args.orbit==42]
    normal=NORMAL if args.frame is None else frame(args.frame)
    out=cut_weights(project(marks,permutations(range(5)),normal))
    print(marks,len(out),*cuts(out))
