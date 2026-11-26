"""Direct four-loop quadratic-lift prism periods from two independent star eliminations."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, permutations
from math import factorial, prod

from flint import fmpq, fmpq_mpoly_ctx, fmpz_mpoly_ctx
from quartic_seven_plane_double import (
    G, NORMAL, Q, fixture, forest_data, polynomial_det, polynomials,
)
from quartic_seven_plane_single import cuts, orbit
from quartic_seven_plane_double_tensor import frame
from quartic_scalar_response import interval
from quartic_planar_cores import cubic_cores
from fourfold_bulk import symmetries


CORE=((0,1),(0,2),(0,3),(1,2),(1,4),(2,5),(3,4),(3,5),(4,5))
STARS=((0,1,2),(4,6,8))
MARKS=((3,3,3,7,7),(3,5,5,5,7),(3,3,3,5,7))
CTX=fmpq_mpoly_ctx.get(tuple(f'w{i}' for i in range(9)))
T=CTX.gens()


@cache
def graph_data(marks,core=CORE):
    edges=fixture(core,marks)
    paths=[]
    for i in range(9):
        start=sum(map(len,paths))
        paths.append(tuple(range(start,start+marks.count(i)+1)))
    graph=G(11,edges)
    graph.C=[[row[i] for i,path in enumerate(paths) for _ in path] for row in G(6,core).C]
    ctx=fmpz_mpoly_ctx.get(tuple(f't{i}' for i in range(14)))
    times=ctx.gens()
    gram=[[sum(x*y*t for x,y,t in zip(row,col,times)) for col in graph.C] for row in graph.C]
    adj=[[(-1)**(i+j)*polynomial_det([row[:i]+row[i+1:] for k,row in enumerate(gram) if k!=j],ctx)
          for j in range(4)] for i in range(4)]
    _,u,dist=forest_data(edges,11)
    u=ctx.from_dict({k:int(v) for k,v in u.terms()})
    assert u==polynomial_det(gram,ctx)
    return graph,ctx,u,adj,dist,paths,None


def star(poly,arms):
    """Substitute t_i=T_j T_k/S, returning the numerator and its minimal S pole."""
    s=sum(T[i] for i in arms)
    slices=defaultdict(dict)
    for power,coef in poly.terms():
        m=list(map(int,power))
        degree=sum(m[i] for i in arms)
        for i in arms:
            m[i]=degree-m[i]
        slices[degree][tuple(m)]=coef
    pole=max(slices)
    result=sum(CTX.from_dict(row)*s**(pole-k) for k,row in slices.items())
    while True:
        quotient,remainder=divmod(result,s)
        if remainder:
            return result,pole
        result=quotient
        pole-=1


def cleared_stars(rows,pole):
    result={}
    for key,poly in rows.items():
        if not poly:
            continue
        for arms in STARS:
            poly,order=star(poly,arms)
            assert order<=pole
            poly*=sum(T[i] for i in arms)**(pole-order)
        result[key]=poly
    return result


def path_coefficients(poly,paths,shorts,long,integrate_long=False):
    # Put t_segment = T_path x_segment, sum(x_segment)=1 on each marked path.
    # The long-path factors have been combined. Short-path powers stay as keys.
    rows=defaultdict(lambda:defaultdict(F))
    for m,c in poly.terms():
        m=tuple(map(int,m))
        power=tuple(sum(m[i] for i in path)+(len(path)-1 if integrate_long and j==long else 0)
                    for j,path in enumerate(paths))
        key=tuple(m[i] for h in shorts for i in paths[h])
        value=F(int(c))
        if integrate_long:
            value*=F(prod(factorial(m[i]) for i in paths[long]),factorial(power[long]))
        rows[key][power]+=value
    return {k:CTX.from_dict({m:fmpq(x.numerator,x.denominator) for m,x in row.items() if x})
            for k,row in rows.items()}


def project(marks,orders,normal=NORMAL):
    """Two-star numerator after all five ordered fractional positions are integrated."""
    data=graph_data(marks)
    return project_currents(marks,data,(polynomials(data,order,normal=normal) for order in orders))


def project_currents(marks,data,current_terms):
    if marks not in MARKS:
        raise ValueError('This projection uses the triple insertion on one prism path.')
    paths=data[5]
    long=next(j for j,path in enumerate(paths) if len(path)==4)
    shorts=[j for j,path in enumerate(paths) if len(path)>1 and j!=long]
    index=0 if long==3 else 2
    outside=set(range(14))-set(paths[long])
    common=None
    total=data[1].constant(0)
    for ds,v in current_terms:
        assert ds==([ds[index],0,0,0] if long==3 else [0,0,ds[index],ds[index]])
        if not ds[index]:
            continue
        if common is None:
            common=data[1].constant(1)
            for f,k in ds[index].factor()[1]:
                if any(any(ex[i] for i in outside) for ex in f.monoms()):
                    common*=f**k
        local,remainder=divmod(ds[index],common)
        assert not remainder and all(not any(ex[i] for i in outside) for ex in local.monoms())
        total+=v*local**2
    if not total:
        return CTX.constant(0)
    adj=data[3][0][0] if long==3 else data[3][2][2]+2*data[3][2][3]+data[3][3][3]
    # The cycle Gram matrix depends on segmented times only through each path total.
    first=[path[0] for path in paths]
    other=set(range(14))-set(first)
    adj=CTX.from_dict({tuple(m[i] for i in first):int(c) for m,c in adj.terms()
                       if not any(m[i] for i in other)})
    adj=cleared_stars({():adj},1)[()]
    currents=cleared_stars(path_coefficients(common,paths,shorts,long),3)
    voltages=cleared_stars(path_coefficients(total,paths,shorts,long,True),1)
    squares={}
    for (j,p),(k,q) in combinations_with_replacement(currents.items(),2):
        key=tuple(a+b for a,b in zip(j,k))
        squares[key]=squares.get(key,0)+(1+(j!=k))*p*q
    result=CTX.constant(0)
    for j,square in squares.items():
        row=CTX.constant(0)
        for k,voltage in voltages.items():
            m=tuple(a+b for a,b in zip(j,k))
            pos=0
            den=1
            for h in shorts:
                size=len(paths[h])
                den*=factorial(sum(m[pos:pos+size])+size-1)
                pos+=size
            row+=voltage*fmpq(prod(map(factorial,m)),den)
        result+=square*row
    # At each star, U_old^{-13} gives S^{13}, the Jacobian gives product(T)/S^3,
    # both currents together cost S^{-6}, the voltage and the adjugate each S^{-1}.
    jac=prod(prod(T[i] for i in arms)*sum(T[i] for i in arms)**2 for arms in STARS)
    return jac*adj*result*prod(T[h]**(len(paths[h])-1) for h in shorts)


@cache
def parallel(p,q):
    assert 0<p<Q and 0<q<Q and p+q>Q
    return F(factorial(Q-p-1)*factorial(Q-q-1),factorial(2*Q-p-q-1))


@cache
def tetra_weight(p):
    return orbit(p),F(prod(factorial(x-1) for x in p),factorial(Q-1))


@cache
def integer_parallel(p,q):
    # All positive beta denominators divide (Q-1)!.
    value = parallel(p,q)*factorial(Q-1)
    assert value.denominator==1
    return value.numerator


@cache
def integer_tetra_weight(p):
    key,weight = tetra_weight(p)
    value = weight*factorial(Q-1)
    assert value.denominator==1
    return key,value.numerator


def cut_weights(poly):
    # The three parallel pairs and the tetrahedron factorial share one global
    # denominator. Sum coefficients with the same rational denominator first.
    bins = defaultdict(int)
    for monomial,c in poly.terms():
        p = tuple(int(x)+1 for x in monomial)
        pp = (p[2]+p[3]-Q,p[1]+p[8]-Q,p[6],p[0],p[5],p[4]+p[7]-Q)
        assert sum(pp)==3*Q
        key,scale = integer_tetra_weight(pp)
        factor = scale*integer_parallel(p[2],p[3])*integer_parallel(p[1],p[8])*integer_parallel(p[4],p[7])
        bins[key,int(c.denominator)] += int(c.numerator)*factor
    denominator = factorial(Q-1)**4
    result = defaultdict(F)
    for (key,d),j in bins.items():
        result[key] += F(j,d*denominator)
    return {p:c for p,c in result.items() if c}


def expected():
    from pathlib import Path
    import json
    records=json.loads(Path(__file__).with_name('quartic_seven_plane_prism_totals.json').read_text())
    assert CORE==cubic_cores()[3] and tuple(map(tuple,records['paths']))==MARKS
    actions=symmetries(CORE)
    assert len(actions)==records['automorphisms']==12
    sizes=[len({tuple(sorted(g[i] for i in m)) for g in actions}) for m in MARKS]
    assert sizes==records['placements']
    totals={}
    for key,values in records['frames'].items():
        rows=[tuple(map(F,pair)) for pair in values]
        totals[key]=tuple(sum(size*row[j]/12 for size,row in zip(sizes,rows)) for j in range(2))
    return records,totals


def five_orbit_totals():
    from quartic_seven_plane_single_projection import quadratic_tensor
    from quartic_seven_plane_prism_paired import MARKS as paired_marks
    records,triple=expected()
    paired=records['paired']
    assert tuple(map(tuple,paired['paths']))==paired_marks
    shapes=(*MARKS,*paired_marks)
    actions=symmetries(CORE)
    covered={tuple(sorted(g[h] for h in m)) for m in shapes for g in actions}
    sizes=[len({tuple(sorted(g[h] for h in m)) for g in actions}) for m in paired_marks]
    assert sizes==paired['placements'] and len(covered)==48
    assert records['normal_orders']*len(covered)==5760
    totals={}
    for name,rows in paired['frames'].items():
        rows=[tuple(map(F,row)) for row in rows]
        assert interval(rows[0],24)[1]<0<interval(rows[1],24)[0]
        totals[name]=tuple(a+sum(size*row[i]/12 for size,row in zip(sizes,rows))
                           for i,a in enumerate(triple[name]))
        assert totals[name]==tuple(map(F,records['five_orbit_totals'][name]))
    assert quadratic_tensor(totals)==tuple(tuple(map(F,row)) for row in records['five_orbit_tensor'])
    # The sparse restrictions contain c(c-1)(2c-1) times a quadratic even
    # about c=1/2. Its leading and central coefficients have the same strict sign.
    rows={key:[*records['frames'][key],*paired['frames'][key],totals[key]] for key in ('2','3')}
    for sign,second,third in zip((-1,-1,-1,-1,1,-1),rows['2'],rows['3'],strict=True):
        a=tuple(sign*(F(y)-5*F(x)) for x,y in zip(second,third))
        b=tuple(sign*(125*F(x)-9*F(y)) for x,y in zip(second,third))
        assert interval(a,24)[0]>0 and interval(b,24)[0]>0
    return totals


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--frame',type=int,choices=(2,3))
    args=parser.parse_args()
    records,totals=expected()
    key=str(args.frame) if args.frame else 'default'
    normal=frame(args.frame) if args.frame else NORMAL
    for i in (range(len(MARKS)) if args.all else [0]):
        w=cut_weights(project(MARKS[i],permutations(range(5)),normal))
        val=cuts(w)
        assert val==tuple(map(F,records['frames'][key][i]))
        print(MARKS[i],len(w),*val,flush=True)
    print('three-orbit sum',*totals[key],flush=True)
    print('five-orbit sum',*five_orbit_totals()[key],flush=True)
