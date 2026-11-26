"""Finite interpolation of seven-planar four-loop normal potential differences."""
from fractions import Fraction
from functools import cache
from itertools import combinations, product
from math import gcd
from random import Random
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from flint import fmpq_mat, fmpz_mpoly_ctx
from four_loop_planar import restricted_symbols
from planar_trace import matchings


PAIRS = tuple(combinations(range(5),2))
FRAMES = ((0,1,4,9,12),(3,-5,0,1,2))


def det3(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def normal_frame(z):
    # The signed complementary triple minors equal z_i-z_j.
    gap = z[3]-z[4]
    q = [(1,0,0),(0,1,0),(0,0,gap)]
    q.append(tuple(-sum(Fraction(int(z[i]-z[4]),gap)*q[i][d] for i in range(3))
                   for d in range(3)))
    q.append(tuple(-sum(col) for col in zip(*q)))
    return q


def gale(q):
    if len(q) != 5 or any(map(sum,zip(*q))):
        raise ValueError("Expected five conserved three-dimensional covectors")
    minors = {pair:(-1)**(sum(pair)-1)*det3(*(p for k,p in enumerate(q) if k not in pair))
              for pair in PAIRS}
    z = tuple(minors[i,4] for i in range(4))+(0,)
    assert all(v == z[i]-z[j] for (i,j),v in minors.items())
    return z


@cache
def monomials(degree):
    return [(a,b,c,d) for a in range(degree+1) for b in range(degree//2+1)
            for c in range(degree//2+1) for d in range(degree//3+1)
            if a+2*b+2*c+3*d == degree]


def pair_words(z,pair,lift_degree):
    i,j = pair
    a,b = z[i],z[j]
    three = [t for k,t in enumerate(z) if k not in pair]
    p,q,r = [3*t-sum(three) for t in three]
    variables = (3*(a+b)-2*sum(three), (a-b)**2, p*p+q*q+r*r, p*q*r)
    factor = ((-1)**(i+j-1)*(a-b)*(three[0]-three[1])
              *(three[0]-three[2])*(three[1]-three[2]))
    return [factor*variables[0]**e*variables[1]**f*variables[2]**g*variables[3]**h
            for e,f,g,h in monomials(5-lift_degree)]


def as_matrix(rows):
    return fmpq_mat([[str(x) for x in row] for row in rows])


@cache
def reciprocal_kernel(lift_degree):
    # Impose the row-sum and both-sided momentum-relation identities as POLYNOMIAL identities.
    ctx = fmpz_mpoly_ctx.get(('z0','z1','z2','z3'))
    z = (*ctx.gens(),ctx.constant(0))
    pairs = {ij:pair_words(z,ij,lift_degree) for ij in PAIRS}
    equations = [[sum(w) for w in zip(*(v for ij,v in pairs.items() if i in ij))]
                 for i in range(5)]
    equations.append([sum(z[i]*z[j]*words[k] for (i,j),words in pairs.items())
                      for k in range(len(monomials(5-lift_degree)))])
    rows = set()
    for row in equations:
        polynomials = [v.to_dict() for v in row]
        for exponent in set().union(*polynomials):
            coeff = tuple(int(p.get(exponent,0)) for p in polynomials)
            divisor = gcd(*coeff)
            rows.add(tuple(x//divisor for x in coeff))
    reduced, rank = as_matrix(sorted(rows)).rref()
    pivots = [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]
    free = [j for j in range(reduced.ncols()) if j not in pivots]
    basis = []
    for col in free:
        v = [Fraction(0)]*reduced.ncols()
        v[col] = Fraction(1)
        for i,pivot in enumerate(pivots):
            v[pivot] = -Fraction(str(reduced[i,col]))
        basis.append(v)
    assert (rank,len(basis)) == ((3,6),(3,4),(1,3))[lift_degree]
    return basis


def normal_row(z,pair,degree):
    words = pair_words(z,pair,degree)
    return [sum(v*c for v,c in zip(words,vector)) for vector in reciprocal_kernel(degree)]


def plane_rows(plane,lift):
    total = tuple(map(sum,zip(*plane)))
    plus = tuple(a+b for a,b in zip(total,lift))
    squared = matchings([*plane,lift])
    linear = Fraction(matchings([*plane,plus])-matchings([*plane,total])-squared,2)
    return restricted_symbols(plane),(linear,),(squared,)


def evaluation_row(plane,normal,outside_lifts,pair):
    if len(plane)!=7 or len(outside_lifts)!=5 or any(map(sum,zip(*plane,*outside_lifts))):
        raise ValueError("Expected seven planar waves and five lifts balancing their sum")
    z = gale(normal)
    mean = Fraction(sum(z),5)
    lift = tuple(sum((t-mean)*v for t,v in zip(z,components)) for components in zip(*outside_lifts))
    result = []
    for d,planar in enumerate(plane_rows(plane,lift)):
        normal_part = normal_row(z,pair,d)
        result.extend(x*y for x,y in product(planar,normal_part))
    assert len(result)==25
    return result


def standard_lifts(plane,z,lift):
    # A rational section of the quotient by normal-to-planar shears.
    center = tuple(-Fraction(sum(c),5) for c in zip(*plane))
    shift = tuple(Fraction(a,z[0]-z[1]) for a in lift)
    return [tuple(a+sign*b for a,b in zip(center,shift)) for sign in (1,-1,0,0,0)]


@cache
def interpolation():
    rng = Random(41259)
    nodes, rows = [], []
    for frame in range(2):
        z = FRAMES[frame]
        normal = normal_frame(z)
        for _ in range(32):
            plane = tuple((rng.randrange(-3,4),rng.randrange(1,4)) for i in range(7))
            lift = (rng.randrange(-2,3),rng.randrange(-2,3))
            actual = standard_lifts(plane,z,lift)
            for pair in PAIRS:
                row = evaluation_row(plane,normal,actual,pair)
                if as_matrix(rows+[row]).rank()>len(rows):
                    rows.append(row)
                    nodes.append((plane,normal,actual,pair))
            if len(rows)==(18 if frame==0 else 25):
                break
        if frame==0:
            # At this frame the three normal row spaces have ranks 4,3,3.
            assert len(rows)==18
    matrix = as_matrix(rows)
    assert matrix.nrows()==25 and matrix.det()
    return nodes,matrix.inv()


def reconstruct(values):
    column = list(values)
    if len(column)!=25:
        raise ValueError("Expected one exact value at each of the 25 interpolation nodes")
    return interpolation()[1]*as_matrix([[v] for v in column])


def trace_control(plane,normal,lifts,pair):
    # A full scalar boundary, evaluated by fusion, without using the interpolation row.
    i,j = pair
    weights = [det3(normal[i],normal[j],p) for p in normal]
    a,b,d = restricted_symbols(plane)
    def minor(r,s):
        return det3(*(p for k,p in enumerate(normal) if k not in (r,s)))
    return (2*a-3*b+5*d)*sum(weights[s]*(minor(i,s)**8-minor(j,s)**8)
                              for s in range(5) if s not in pair)


def cross(p,q):
    return tuple(p[(i+1)%3]*q[(i+2)%3]-p[(i+2)%3]*q[(i+1)%3] for i in range(3))


NORMAL_COEFFICIENTS = (
    ('35/54', '-71/544320', '-567557/3772137600000'),
    ('7/4', '1/20160', '-238061/1676505600000'),
    ('-14/81', '-181/1632960', '-14001731/203695430400000'),
    ('35/486', '-13/4898880', '-169219/19096446600000'),
    ('35/54', '13/136080', '42409/5867769600000'),
    ('7/324', '-1/1632960', '-8220757/2851736025600000'),
)


def normal_factors(q,pair):
    # Fused eleven-field trace, integrated genus-three hexagon, integrated genus-two theta.
    row = normal_row(gale(q),pair,0)
    return tuple(sum(x*Fraction(y) for x,y in zip(row,col)) for col in zip(*NORMAL_COEFFICIENTS))


def integrated_normal_factors(q,pair):
    from quartic_seven_plane_hexagon import normal_component as hexagon
    from quartic_seven_plane_theta import normal_component as theta
    i,j = pair
    v = cross(q[i],q[j])
    z = gale(q)
    fused = sum(det3(q[i],q[j],p)*((z[s]-z[i])**8-(z[s]-z[j])**8)
                for s,p in enumerate(q))
    return fused,hexagon(q,i,j,direction=v),theta(q,target=pair,direction=v)


PLANAR_COEFFICIENTS = (
    ('-22599229/3026332108800','143249/23351328000','1/15552'),
    ('2104429/8826801984000','-21479/108972864000','0'),
    ('-16229009/4035109478400','53/16016000','1/38880'),
)


def required_bare(plane,normal,pair):
    """Necessary bare coefficients in 1, zeta(3), zeta(5) and the undetermined trace."""
    from quartic_scalar_response import total_channels,BIPARTITE
    a,b,d = restricted_symbols(plane)
    w,h,t = normal_factors(normal,pair)
    result = [w*(Fraction(b,576)*k+Fraction(8*a-b,9)*j) for k,j in total_channels(BIPARTITE)]
    hx0,hx3,th = (sum(x*Fraction(y) for x,y in zip((a,b,d),column))
                  for column in zip(*PLANAR_COEFFICIENTS))
    result[0] -= h*hx0+t*th
    result[1] -= h*hx3
    return (*result,w*(22*a-5*b+18*d))


def check_integrals():
    from quartic_seven_plane_hexagon import restrict,planar_factors
    from quartic_seven_plane_theta import NORMAL,planar_factor
    from quartic_seven_plane_target import PLANES

    probes = (((0,1),(0,2),(0,3),(1,2)),((0,1),(0,2)))
    rows,integrals = [],[]
    for frame,pairs in zip(FRAMES,probes):
        q = normal_frame(frame)
        for pair in pairs:
            rows.append(normal_row(frame,pair,0))
            integrals.append(integrated_normal_factors(q,pair))
    assert as_matrix(rows).inv()*as_matrix(integrals) == as_matrix(NORMAL_COEFFICIENTS)
    assert as_matrix(integrals).rank() == 3
    # Two held-out geometries, including the twelve-wave transfer fixture.
    for q,pair in ((NORMAL,(0,1)),(normal_frame((-2,1,-1,2,5)),(1,4))):
        assert normal_factors(q,pair) == integrated_normal_factors(q,pair)
    with TemporaryDirectory() as tmp:
        native = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(Path(__file__).with_name('nine_sextics.cpp')),
             '-o',str(native)],check=True)
        coefficients = restrict(native)
    test_planes = [*PLANES,((1,0),)*3+((0,1),)*3+((-3,-3),)]
    original = [(*planar_factors(p,coefficients),planar_factor(p)) for p in test_planes]
    assert as_matrix([restricted_symbols(p) for p in test_planes]).inv()*as_matrix(original) == as_matrix(PLANAR_COEFFICIENTS)
    other = ((1,3),(2,3),(1,5),(-5,1),(-1,2),(3,1),(2,1))
    assert as_matrix([restricted_symbols(other)])*as_matrix(PLANAR_COEFFICIENTS) == as_matrix([
        (*planar_factors(other,coefficients),planar_factor(other))])
    for plane,expected in zip(PLANES,(Fraction(-81473804382978468857648,2698059375),
                                      Fraction(-96109434930082264435136,1618835625))):
        _,_,weight_five,ambiguity = required_bare(plane,NORMAL,(0,1))
        assert weight_five == -7*expected and not ambiguity
    print('All three unrestricted normal factors regenerated and checked on two additional configurations.')


def check():
    for f in FRAMES:
        assert gale(normal_frame(f)) == tuple(t-f[4] for t in f)
    nodes,_ = interpolation()
    coeff = reconstruct(trace_control(*node) for node in nodes)
    assert all(coeff[j,0]==0 for j in range(18,25))
    rng = Random(7142)
    for _ in range(10):
        z = tuple(rng.sample(range(-9,10),5))
        q = normal_frame(z)
        plane = tuple((rng.randrange(-2,3),rng.randrange(1,5)) for i in range(7))
        ell = [tuple(rng.randrange(-2,3) for i in range(2)) for j in range(4)]
        ell.append(tuple(-sum(col) for col in zip(*plane,*ell)))
        pair = rng.choice(PAIRS)
        row = evaluation_row(plane,q,ell,pair)
        shear = ((1,2,0),(-1,0,3))
        shifted = [tuple(e+sum(a*b for a,b in zip(p,h)) for e,h in zip(x,shear)) for p,x in zip(q,ell)]
        assert row == evaluation_row(plane,q,shifted,pair)
        assert Fraction(str((as_matrix([row])*coeff)[0,0])) == trace_control(plane,q,ell,pair)
    print('25 independent scalar readings. Planar multiplicities 3,1,1; normal dimensions 6,4,3.')


if __name__=='__main__':
    check()
    check_integrals()
