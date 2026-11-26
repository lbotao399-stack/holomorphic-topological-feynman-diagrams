"""Evaluate the fourfold Ward prism orbit with its five hard waves on two edges."""
from collections import defaultdict
from argparse import ArgumentParser
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement
from math import factorial, lcm, prod
from pathlib import Path
import json

from flint import fmpq, fmpq_mat, fmpz_mpoly_ctx
from fourfold_quartic_orbits import grouped_orbits, separated_squares, squared_cofactor_polynomial
from quartic_seven_plane_double import primitive
from quartic_seven_plane_prism import CORE, CTX, MARKS, STARS, T, cut_weights, graph_data, star
from quartic_seven_plane_single import cuts


def split_paths(poly, marked, paths, stars=STARS):
    """Separate the segment variables on the marked paths from the free stars."""
    inside = [i for h in marked for i in paths[h]]
    outside = [j for j in range(9) if j not in marked]
    path_lengths = tuple(len(paths[h]) for h in marked)
    slices = defaultdict(dict)
    for ex,c in poly.terms():
        x = tuple(int(ex[i]) for i in inside)
        y = tuple(int(ex[paths[h][0]]) for h in outside)
        start = 0
        degrees = []
        for length in path_lengths:
            degrees.append(sum(x[start:start+length]))
            start += length
        slices[tuple(degrees)][x,y] = int(c)
    local = fmpz_mpoly_ctx.get(tuple('x'+str(i) for i in range(len(inside))))
    for degree,sparse in sorted(slices.items()):
        left = sorted({x for x,y in sparse})
        right = sorted({y for x,y in sparse})
        matrix = fmpq_mat([[sparse.get((x,y),0) for x in left] for y in right])
        rows,rank = matrix.rref()
        for i in range(rank):
            pivot = next(j for j in range(len(left)) if rows[i,j])
            l,coef = primitive([(x,F(str(rows[i,j]))) for j,x in enumerate(left) if rows[i,j]])
            r = {}
            for j,y in enumerate(right):
                if matrix[j,pivot]:
                    whole = [0]*9
                    for h,ex in zip(outside,y):
                        whole[h] = ex
                    r[tuple(whole)] = matrix[j,pivot]
            polynomial = CTX.from_dict(r)
            poles = []
            for arms in stars:
                polynomial,order = star(polynomial,arms)
                poles.append(order)
            yield degree,local.from_dict(dict(l)),polynomial,coef,poles


def project(data,main,voltage,denominator):
    paths = data[5]
    marked = [i for i,path in enumerate(paths) if len(path)>1]
    if len(marked)!=2:
        raise ValueError('This reduction requires two marked prism paths')
    m,n = (len(paths[k]) for k in marked)
    first,second = (list(split_paths(poly,marked,paths)) for poly in (main,voltage))
    poles = [max(p[4][i] for p in first)+max(p[4][i] for p in second) for i in range(2)]
    groups = defaultdict(lambda:CTX.constant(0))
    for a,x,r,c,s in first:
        for b,y,h,v,w in second:
            degree = tuple(i+j for i,j in zip(a,b))
            # Dirichlet moments for the 4-part and 3-part path splits.
            num = sum(int(coef)*prod(factorial(int(i)) for i in ex) for ex,coef in (x*y).terms())
            if num:
                norm = F(num,factorial(degree[0]+m-1)*factorial(degree[1]+n-1))
                orders = tuple(i+j for i,j in zip(s,w))
                groups[degree,orders] += fmpq(str(norm*c*v/denominator))*r*h
    result = CTX.constant(0)
    for (degree,orders),poly in groups.items():
        result += (poly*T[marked[0]]**(degree[0]+m-1)*T[marked[1]]**(degree[1]+n-1)
                   *prod(sum(T[i] for i in arms)**(p-q) for arms,p,q in zip(STARS,poles,orders)))
    for arms,pole in zip(STARS,poles):
        s = sum(T[i] for i in arms)
        if pole>10:
            result,rem = divmod(result,s**(pole-10))
            assert not rem
        elif pole<10:
            result *= s**(10-pole)
        result *= prod(T[i] for i in arms)
    return result


def project_squares(data,squares,voltage,denominator,stars=STARS):
    paths = data[5]
    marked = [i for i,path in enumerate(paths) if len(path)>1]
    dimensions = tuple(len(paths[k])-1 for k in marked)
    star_sums = [sum(T[i] for i in arms) for arms in stars]

    @cache
    def factors(degrees,*powers):
        return prod(T[k]**d for k,d in zip(marked,degrees))*prod(s**a for s,a in zip(star_sums,powers))

    result = CTX.constant(0)
    for p,weight in squares:
        rows = list(split_paths(p,marked,paths,stars))
        weights = list(split_paths(weight*voltage,marked,paths,stars))
        for i,j in combinations_with_replacement(range(len(rows)),2):
            a,x,r,c,s = rows[i]
            b,y,h,v,t = rows[j]
            local = x*y
            degree = tuple(i+j for i,j in zip(a,b))
            poles = tuple(i+j for i,j in zip(s,t))
            coefficient = CTX.constant(0)
            for exponent,ell,outside,scalar,extra in weights:
                moment = sum(int(f)*prod(factorial(int(k)) for k in ex) for ex,f in (local*ell).terms())
                if not moment:
                    continue
                norm = F(moment,prod(factorial(i+j+k) for i,j,k in zip(degree,exponent,dimensions)))
                powers = tuple(11-p-e for p,e in zip(poles,extra))
                assert all(a>=0 for a in powers)
                coefficient += fmpq(str(norm*scalar))*outside*factors(exponent,*powers)
            scale = c*v*F(1+(i!=j),denominator)
            result += r*h*(coefficient*fmpq(str(scale)))*factors(tuple(i+j for i,j in zip(degree,dimensions)))
    for s in star_sums:
        result,rem = divmod(result,s)
        assert not rem
    return result*prod(T[i] for arms in stars for i in arms)


def project_slices(data,squares,voltage,denominator,stars=STARS):
    """Yield disjoint marked-degree slices after the empty-star divisions."""
    paths = data[5]
    marked = tuple(i for i,p in enumerate(paths) if len(p)>1)
    dimensions = tuple(len(paths[i])-1 for i in marked)
    sums = [sum(T[i] for i in star) for star in stars]
    jacobian = prod(T[i] for star in stars for i in star)
    @cache
    def star_powers(powers):
        return prod(s**e for s,e in zip(sums,powers,strict=True))

    from ordered_path_moments import PositionMoments
    by_degree = defaultdict(list)
    for root,weight in squares:
        pieces = list(split_paths(root,marked,paths,stars))
        weights = list(split_paths(weight*voltage,marked,paths,stars))
        moments = PositionMoments([w[1] for w in weights])
        for i,j in combinations_with_replacement(range(len(pieces)),2):
            a,x,r,c,s = pieces[i]
            b,y,h,f,t = pieces[j]
            local = x*y
            degree = tuple(k+l+n for k,l,n in zip(a,b,dimensions))
            coefficients = defaultdict(list)
            for (u,_,outside,scalar,extra),moment in zip(weights,moments(local),strict=True):
                final = tuple(k+l for k,l in zip(degree,u))
                if moment:
                    powers = tuple(11-p-q-e for p,q,e in zip(s,t,extra,strict=True))
                    assert all(p>=0 for p in powers)
                    value = c*f*scalar*F((1+(i!=j))*moment,denominator*prod(map(factorial,final)))
                    coefficients[final].append((value,powers,outside))
            for final,terms in coefficients.items():
                by_degree[final].append((r,h,terms))
    for degree,blocks in sorted(by_degree.items()):
        poly = CTX.constant(0)
        for left,right,terms in blocks:
            c = sum(fmpq(str(value))*star_powers(powers)*outside for value,powers,outside in terms)
            poly += left*right*c
        # The divisors use exclusively the empty-star variables. Exact
        # division therefore commutes with extraction of the marked degrees.
        for star_sum in sums:
            poly,rem = divmod(poly,star_sum)
            assert not rem
        yield poly*jacobian*prod(T[i]**k for i,k in zip(marked,degree))


def prism_group(marks,order,weights,check=False):
    data = graph_data(marks)
    graph,ctx,_,_,dist,_,_ = data
    v = -sum(fmpq(str(weights[i][j]))*dist[6+i][6+j]
             for i in range(5) for j in range(5) if weights[i][j])/2
    if not v:
        return F(),F()
    den = lcm(*(int(c.denominator) for _,c in v.terms()))
    voltage = ctx.from_dict({m:int(c*den) for m,c in v.terms()})
    squares = separated_squares(graph,order)
    if not squares:
        return F(),F()
    result = project_squares(data,squares,voltage,den)
    if check:
        main = squared_cofactor_polynomial(graph,order)
        assert result == project(data,main,voltage,den)
    return cuts(cut_weights(result))


def recorded_orbit():
    path = Path(__file__).with_name('fourfold_quartic_periods.json')
    shape, = json.loads(path.read_text())['prism']['shapes']
    marks = tuple(shape['hard_paths'])
    assert marks == MARKS[0]
    groups = grouped_orbits(CORE,marks)
    recorded = [(tuple(key),tuple(map(F,row))) for key,row in shape['gale_groups']]
    assert [key for key,_ in groups] == [key for key,_ in recorded]
    assert len(recorded) == 30 and sum(any(row) for _,row in recorded) == 24
    channels = tuple(sum(col) for col in zip(*(row for _,row in recorded)))
    assert channels == tuple(map(F,shape['total'])) == (
        -F(233682040128782429,80879986579392000000000), F(85714400153,35661369744000000))
    # The upper integral-test bound on the zeta(3) tail fixes the sign.
    upper = sum(F(1,k**3) for k in range(1,65)) + F(1,2*64**2)
    assert channels[1] > 0 and channels[0] + upper*channels[1] < 0
    return marks,groups,recorded,channels


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--group',type=int,choices=range(30))
    mode.add_argument('--table-only',action='store_true')
    args = parser.parse_args()
    marks,groups,recorded,channels = recorded_orbit()
    indices = [] if args.table_only else range(30) if args.group is None else [args.group]
    for i in indices:
        key,(order,weights) = groups[i]
        row = prism_group(marks,order,weights,check=i==0)
        assert row == recorded[i][1]
        print(key,*row,flush=True)
    print('Recorded two-edge prism orbit',*channels)
