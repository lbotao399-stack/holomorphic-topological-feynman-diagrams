"""Contract the fourfold four-loop Ward integrand with the 3872 tetrahedron periods."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import factorial, lcm

from flint import fmpq, fmpq_mpoly_ctx
from fourfold_quartic_orbits import grouped_orbits, separated_squares, squared_cofactor_polynomial
from quartic_seven_plane_single import CORE, cuts, graph_data
from quartic_seven_plane_single_projection import Bubble, split_variables, orbital, factorials, Q, SCALE, shape_orbits


def project(data,main,voltage):
    """Take the local two-terminal bubble moment before multiplying exterior polynomials."""

    paths, split = data[5:]
    nleft = sum(map(len,paths[:4]))
    bubble = Bubble(nleft,split)
    unit = bubble.register((((0,)*nleft,1),))
    right = fmpq_mpoly_ctx.get(tuple('x'+str(i) for i in range(sum(map(len,paths))-nleft)))

    def parts(poly):
        return [(bubble.register(l),right.from_dict(dict(r)),degree,coef)
                for l,r,degree,coef in split_variables(poly,nleft)]

    main_parts, voltage_parts = parts(main), parts(voltage)
    good, bad = {}, {}
    for i,r,degree,coef in main_parts:
        for j,s,other,weight in voltage_parts:
            moment, lost = bubble.project(i,unit,j)
            product = r*s*fmpq(str(coef*weight))
            if moment:
                good[degree+other] = good.get(degree+other,0)+product*(fmpq(moment)/SCALE)
            for pos,x in lost:
                bad[pos] = bad.get(pos,0)+product*x
    assert not any(bad.values())
    return tetrahedron_projection(paths,nleft,good)


def tetrahedron_projection(paths,nleft,good):
    weights = defaultdict(F)
    for degree,poly in good.items():
        power = degree+nleft-Q
        for ex,x in poly.terms():
            k = tuple(map(int,ex))
            rest = tuple(sum(k[i-nleft]+1 for i in path) for path in paths[4:])
            assert power+sum(rest) == 3*Q
            weights[orbital((power,*rest))] += F(int(x.numerator),int(x.denominator))*F(
                factorial(power-1)*factorials(k),factorial(Q-1))
    return cuts({k:v for k,v in weights.items() if v})


def project_squares(data,squares,voltage):
    paths,split = data[5:]
    nleft = sum(map(len,paths[:4]))
    bubble = Bubble(nleft,split)
    right = fmpq_mpoly_ctx.get(tuple('x'+str(i) for i in range(sum(map(len,paths))-nleft)))
    good,bad = {},{}

    def parts(poly):
        return [(bubble.register(l),right.from_dict(dict(r)),degree,coef)
                for l,r,degree,coef in split_variables(poly,nleft)]

    for polynomial,weight in squares:
        rows = parts(polynomial)
        factors = parts(weight*voltage)
        for i,j in combinations_with_replacement(range(len(rows)),2):
            a,r,degree,c = rows[i]
            b,s,other,d = rows[j]
            finite,singular = {},{}
            for idx,t,shift,k in factors:
                value,lost = bubble.project(a,b,idx)
                scalar = fmpq(str(c*d*k*(1+(i!=j))))
                if value:
                    finite[shift] = finite.get(shift,0)+t*(scalar*value/SCALE)
                for pos,x in lost:
                    singular[pos] = singular.get(pos,0)+t*(scalar*x)
            exterior = r*s
            for shift,h in finite.items():
                total = degree+other+shift
                good[total] = good.get(total,0)+exterior*h
            for pos,h in singular.items():
                bad[pos] = bad.get(pos,0)+exterior*h
    assert not any(bad.values())
    return tetrahedron_projection(paths,nleft,good)


def integrate(marks,group=None,separate=False):
    data = graph_data(marks)
    graph,ctx,_,_,dist,_,_ = data
    total = [F(),F()]
    for index,(_,(order,weights)) in enumerate(grouped_orbits(CORE,marks)):
        if group is not None and index!=group:
            continue
        v = -sum(fmpq(str(weights[i][j]))*dist[6+i][6+j]
                 for i in range(5) for j in range(5) if weights[i][j])/2
        if not v:
            continue
        den = lcm(*(int(c.denominator) for _,c in v.terms()))
        voltage = ctx.from_dict({m:int(c*den) for m,c in v.terms()})
        if separate:
            result = project_squares(data,separated_squares(graph,order),voltage)
        else:
            main = squared_cofactor_polynomial(graph,order)
            result = project(data,main,voltage) if main else (F(),F())
        for i in range(2):
            total[i] += result[i]/den
    return tuple(total)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('shape',help='3+2, 2+3 or a zero-based single-double shape index')
    parser.add_argument('--group',type=int)
    parser.add_argument('--separate',action='store_true')
    args = parser.parse_args()
    if args.shape in ('3+2','2+3'):
        marks = (0,0,0,8,8) if args.shape=='3+2' else (0,0,8,8,8)
    else:
        all_marks = [m for m,_ in shape_orbits()]
        if not args.shape.isdecimal() or int(args.shape)>=len(all_marks):
            parser.error('Expected one of the 36 admissible single-double shapes')
        marks = all_marks[int(args.shape)]
    if args.group is not None and args.group not in range(len(grouped_orbits(CORE,marks))):
        parser.error('Group number outside the chosen shape')
    result = integrate(marks,args.group,args.separate)
    expect = {'3+2':(-F(10068993766933,17300531888640000000),F(184663,381405024000)),
              '2+3':(-F(618691001569,469890989568000000),F(5106053,4661616960000))}
    if args.shape in expect and args.group is None:
        assert result == expect[args.shape]
    print('Single-double four-loop hard orbit',marks,'rational and zeta(3) channels',*result)
