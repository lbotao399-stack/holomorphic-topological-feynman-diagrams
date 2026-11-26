"""Planar primary binary at a double zero after the massive shear z -> z + uxy."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from planar_cubic_primary import CORNERS, ROOT, labeled_vertex


@cache
def packets(terms, n):
    out = []
    for slots in combinations_with_replacement(range(len(terms)), n):
        weight = prod(terms[i][1] for i in slots)
        for i in set(slots):
            weight /= factorial(slots.count(i))
        out.append((tuple(terms[i][0] for i in slots), weight))
    return out


def stencils():
    # Each background mark has at least two derivatives. The three stencils
    # start at degrees 2,3,4 and their selected degrees already exhaust 6g.
    quadratic = (((0,0,1), F(1)),)
    cubic = tuple((p,F(prod(p),8)) for p in CORNERS)
    quartic = tuple(((x,y,0), F(1,2) if x*y else F(-1))
                    for x,y in product((-1,0,1), repeat=2) if x or y)
    # The omitted origin in 2(cosh x - 1)(cosh y - 1) has zero marked vertex.
    return quadratic, cubic, quartic


def grouped(genus, count):
    quadratic, cubic, quartic = stencils()
    letters = set(p for terms in (quadratic,cubic,quartic) for p,w in terms)
    actions = [{p: tuple(s*t for s,t in zip(signs,(p[i] for i in order))) for p in letters}
               for order in permutations(range(3)) for signs in CORNERS]
    weights = defaultdict(F)
    for (a,wa),(b,wb),(d,wd) in product(packets(quadratic,count), packets(cubic,2*(genus-count)),
                                         packets(quartic,count)):
        key = min(tuple(sorted(action[p] for p in a+b+d)) for action in actions)
        weights[key] += wa*wb*wd
    return {key:w for key,w in weights.items() if w}


def polynomial(genus, exe):
    out = []
    for count in range(genus+1):
        rows = grouped(genus,count)
        vertices = labeled_vertex(genus,rows,exe)
        value = tuple(sum(w*x for w,x in zip(rows.values(),col)) for col in zip(*vertices))
        out.append(value)
        print('genus',genus,'quartic marks',count,'graph columns',value,flush=True)
    return out


def ward(genus, count, potential):
    # The generator xy d_z consumes one x and y. The open z index adds one
    # z component to the symbol. The background degrees are (2g-1,2g-1,2g+1).
    quadratic,cubic,quartic = stencils()
    shear = tuple(((i,j,0),F(i*j,4)) for i,j in product((-1,1),repeat=2))
    packets_by_type = (packets(quadratic,count+1), packets(cubic,2*(genus-count)-1),
                       packets(quartic,count))
    letters = {p for terms in (quadratic,cubic,quartic,shear) for p,w in terms}
    actions = [({p: tuple(s*p[i] for s,i in zip(signs,order)) for p in letters}, signs[2])
               for order in ((0,1,2),(1,0,2)) for signs in CORNERS]
    weights = defaultdict(F)
    for (v,wv),(a,wa),(b,wb),(d,wd) in product(shear,*packets_by_type):
        fields, sign = min(((act[v], tuple(sorted(act[p] for p in a+b+d))), sign)
                           for act,sign in actions)
        weights[fields] += sign*wv*wa*wb*wd
    result = None
    for (v,normals),weight in weights.items():
        if not weight:
            continue
        fields = [(0,1,0,0,0),(0,0,*v),(1,0,0,0,0)]+[(0,0,*p) for p in normals]
        root = tuple(-sum(col) for col in zip(*fields))
        terms = []
        for phi in potential(fields+[root]):
            terms.append(weight*(root[4]*(phi[0][0]-phi[-1][0])-root[0]*(phi[1][4]-phi[-1][4])))
        if result is None:
            result = [F(0)]*len(terms)
        result = [a+b for a,b in zip(result,terms)]
    return result


def theta_check(terms):
    from bulk_eight import SHAPES, potentials
    def by_shape(points):
        return (potentials(points,shapes=(shape,)) for shape in SHAPES)
    for k in range(2):
        shapes = ward(2,k,by_shape)
        expected = (4-2*k)*terms[k][0]+(2*k+2)*terms[k+1][0]
        print('theta integrals', k, shapes, 'sum',sum(shapes),flush=True)
        assert sum(shapes) == -expected


def branch_preserving_check(two, three):
    import sympy as sp
    r = sp.Symbol('r', real=True)
    def columns(terms):
        g = len(terms)-1
        return [sp.Poly(sum(w*r**(g-i)*(1+r)**i for i,w in enumerate(col)),r)
                for col in zip(*terms)]
    divisor = columns(two)[0]*155520
    assert divisor.as_expr() == 17*r**2+44*r+30 and divisor.discriminant() == -104
    remainder = [p.rem(divisor) for p in columns(three)]
    determinant = sp.Matrix([p.all_coeffs() for p in remainder]).det()
    assert determinant == -F(1290669321259,4721004731437682688000000)
    print('independent real linear remainders at the two-loop roots',determinant,flush=True)

    # A nonlinear automorphism preserving both critical branches, composed with
    # a diagonal scaling and a linear shear into the massive coordinate.
    x,y,z,t,u,c = sp.symbols('x y z t u c')
    X,Y,Z = 2*(x+z), y/6, 3*z
    Z += u*X*Y + X**2*Y
    assert sp.Matrix([X,Y,Z]).jacobian((x,y,z)).det().expand() == 1
    potential = (Z**2+c*X**2*Y**2)/2
    jet = potential.subs({x:t*x,y:t*y,z:t**2*z}, simultaneous=True).expand().coeff(t,4)
    assert (2*jet-(3*z+u*x*y/3)**2-c*x*x*y*y/9).expand() == 0
    transverse = sp.hessian(potential,(y,z)).subs({y:0,z:0}).det()
    assert (transverse - c*x*x).expand() == 0


def rank_one_holdout(genus, exe):
    """Direct shear (x+y)^2 in the regular frame s=x+y, t=y (no discriminant division)."""
    from collections import Counter
    from critical_line_quadratic_shear import monomial_stencil
    # B = z²/2 + zs² + (s⁴+s²t²-2st³+t⁴)/2.
    quadratic = stencils()[0]
    cubic = tuple(((i,0,j),F(j,2)*(1 if i else -2)) for i,j in product((-1,0,1),(-1,1)))
    quartic = {0:1, 1:-2, 2:1, 4:1}  # index is the s degree
    weights = defaultdict(F)
    for count in range(genus+1):
        for degrees in combinations_with_replacement(quartic,count):
            if sum(4-i for i in degrees) != 2*genus:
                continue
            selections = [(quadratic,count),(cubic,2*(genus-count))]
            for i,n in Counter(degrees).items():
                selections.append((monomial_stencil(i),n))
            coeff = prod(quartic[i] for i in degrees)
            for event in product(*(packets(stencil,n) for stencil,n in selections)):
                waves = tuple(sorted(p for group,_ in event for p in group))
                weights[waves] += coeff*prod(w for _,w in event)
    letters = set(p for row in weights for p in row)
    actions = [{p:tuple(s*p[j] for s,j in zip(signs,perm)) for p in letters}
               for perm in permutations(range(3)) for signs in CORNERS]
    grouped_waves = defaultdict(F)
    for waves,w in weights.items():
        key = min(tuple(sorted(action[p] for p in waves)) for action in actions)
        grouped_waves[key] += w
    rows = {p:w for p,w in grouped_waves.items() if w}
    values = labeled_vertex(genus,rows,exe)
    actual = tuple(sum(w*x for w,x in zip(rows.values(),col)) for col in zip(*values))
    expected = {2:(F(49,25920),0), 3:(-F(13414400183,60526642176000),F(511429,2918916000))}
    assert actual == expected[genus]
    print('rank-one shear in a regular frame',genus,actual,flush=True)


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--hexagon',action='store_true',help='integrate the separated three-loop hexagon')
    if parser.parse_args().hexagon:
        from bulk_three_hexagon import potentials
        expected = [-F(1,14515200),F(1,2721600),-F(1,1451520)]
        for k,value in enumerate(expected):
            observed = ward(3,k,lambda p:(potentials(p),))
            assert observed == [value]
            print('three-loop hexagon, quartic donors',k,value,flush=True)
        raise SystemExit(0)
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        calculated = []
        for g in (2,3):
            terms = polynomial(g,exe)
            calculated.append(terms)
            rank_one_holdout(g,exe)
            if g == 2:
                assert terms == [(F(1,51840),0),(-F(1,9720),0),(F(1,5184),0)]
                theta_check(terms)
            else:
                from critical_line_binary import cubic_jet_polynomial
                assert terms[-1] == tuple(8*x for x in cubic_jet_polynomial()[2,2,2])
        branch_preserving_check(*calculated)
