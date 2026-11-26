"""Primary planar BCOV brackets in a ternary cubic normal potential."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
import json

from elementary_scalar import ROOT
from longitudinal_theta import primary_coefficient
from three_loop_soft_jet import DOUBLE_SOFT, orbit_integers

UNIT = [(1,0,0),(0,1,0),(0,0,1)]
CORNERS = tuple(product((-1,1),repeat=3))
HARD = [(1,0,0,0,0),(0,1,0,0,0)]
FERMAT = (F(1163706689,25609494822912000),-F(74897,1976041267200))


def symmetry_tables():
    sites = tuple((*CORNERS,*UNIT,*(tuple(-v for v in p) for p in UNIT)))
    tables = []
    for order in permutations(range(3)):
        for signs in CORNERS:
            tables.append({p:tuple(s*p[i] for i,s in zip(order,signs)) for p in sites})
    return tables


def cubic_orbits(genus):
    # In (x³+y³+z³)/6 + b xyz exactly 2g derivatives are used on each axis.
    # If r labels use xyz, the rest use n=(2g-r)/3 pure cubic marks on each axis.
    tables = symmetry_tables()
    for r in range(2*genus%3,2*genus+1,3):
        n = (2*genus-r)//3
        weights = Counter()
        for corners in combinations_with_replacement(CORNERS,r):
            cube_weight = F(prod(prod(p) for p in corners),
                            8**r*prod(factorial(k) for k in Counter(corners).values()))
            for plus in product(range(n+1),repeat=3):
                points = list(corners)
                weight = cube_weight
                for p,k in zip(UNIT,plus):
                    points.extend([p]*k+[tuple(-x for x in p)]*(n-k))
                    weight *= F((-1)**(n-k),2**n*factorial(k)*factorial(n-k))
                # Every normal coordinate has even total degree. Reflections and
                # permutations of the three axes preserve the full vertex.
                orbit = min(tuple(sorted(table[p] for p in points)) for table in tables)
                weights[orbit] += weight
        yield r,{orbit:w for orbit,w in weights.items() if w}


def labeled_vertex(genus, orbits, exe, hard=HARD):
    rows = [[(0,0,*p) for p in normal]+list(hard) for normal in orbits]
    return labeled_planar_rows(genus,rows,exe)


def labeled_planar_rows(genus,rows,exe):
    """One tangent vector in the last slot, polarized along the first axis."""
    if genus == 2:
        return [(primary_coefficient(row[:-1],row[-1],HARD[0]),F(0)) for row in rows]
    data = json.loads((ROOT/'operator'/'three_loop_coefficients.json').read_text())
    assert tuple(map(tuple,data['basis_words'])) == DOUBLE_SOFT
    coefficients = [list(map(F,col)) for col in zip(*data['rational_and_zeta'])]
    fixtures = [row+[tuple(-sum(c) for c in zip(*row))] for row in rows]
    return [tuple(points[-1][0]*sum(x*y for x,y in zip(row,col)) for col in coefficients)
            for points,row in zip(fixtures,orbit_integers(exe,fixtures))]


def hard_normal_shears(exe):
    # The first two tangent momenta stay independent, and the vector
    # polarization still pairs to -1 with the balancing root.
    shifted = [(1,0,2,-1,3),(0,1,-4,2,1)]
    for genus in (2,3):
        _,points = next(cubic_orbits(genus))
        bare = labeled_vertex(genus,points,exe)
        assert any(any(value) for value in bare)
        assert labeled_vertex(genus,points,exe,shifted) == bare


def hesse(genus,exe):
    polynomial = {}
    for power,orbits in cubic_orbits(genus):
        vertices = labeled_vertex(genus,orbits,exe)
        polynomial[power] = tuple(sum(w*a for w,a in zip(orbits.values(),col))
                                  for col in zip(*vertices))
    return polynomial


def hessian_pencil(cubic, normal):
    """The reduced two-loop amplitude and the three-loop amplitude in Fermat units."""
    import sympy as sp
    parameter = sp.Dummy('t')
    if len(normal)!=3:
        raise ValueError('Supply the three normal coordinate variables.')
    form = sp.Poly(cubic,*normal)
    if form.is_zero:
        return sp.S.Zero,sp.S.Zero
    if any(sum(degrees)!=3 for degrees,_ in form.terms()):
        raise ValueError('The normal potential must be a homogeneous cubic.')
    determinant = sp.hessian(cubic,normal).det()
    pencil = sp.Poly(sp.hessian(cubic+parameter*determinant,normal).det(),parameter)
    monomial, coefficient = form.terms()[0]
    linear = sp.Poly(pencil.nth(1),*normal)
    two = sp.cancel(linear.coeff_monomial(monomial)/(12*51840*coefficient))
    assert linear.as_expr().expand() == sp.expand(12*51840*two*cubic)
    second = sp.Poly(pencil.nth(2)+12*51840*two*determinant,*normal)
    ratio = sp.cancel(-second.coeff_monomial(monomial)/(6*coefficient))
    assert sp.expand(second.as_expr()+6*ratio*cubic) == 0
    assert sp.expand(pencil.nth(3)-48*(51840*two)**2*cubic-2*ratio*determinant) == 0
    return two,ratio


def planar_coefficients(cubic,normal):
    """Return reduced Q2 and the rational, zeta(3) columns of reduced Q3."""
    two,ratio=hessian_pencil(cubic,normal)
    return two,tuple(ratio*c for c in FERMAT)


def weierstrass_graph(genus,exe):
    # y²z-x³-xz² at g=2, y²z-x³-z³ at g=3. Each nonzero allocation is forced.
    def mixed(i,j):  # 2(cosh(w_i)-1)sinh(w_j) has cubic jet w_i²w_j
        rows={}
        for s in (-1,1):
            q=tuple(s*x for x in UNIT[j])
            rows[q]=F(-s)
            for t in (-1,1):
                rows[tuple(x+t*y for x,y in zip(q,UNIT[i]))]=F(s,2)
        return rows
    def negative_cube(i):
        return {tuple(s*x for x in UNIT[i]):F(-3*s) for s in (-1,1)}
    def assignments(rows,count):
        for points in combinations_with_replacement(rows,count):
            yield points,prod(rows[p] for p in points)/prod(factorial(n) for n in Counter(points).values())
    extra={p:-a for p,a in mixed(2,0).items()} if genus==2 else negative_cube(2)
    sectors=[assignments(mixed(1,2),genus),assignments(negative_cube(0),genus-1),assignments(extra,1)]
    orbits=Counter()
    for terms in product(*sectors):
        pts=tuple(p for points,_ in terms for p in points)
        key=min(tuple(sorted(tuple(i*j for i,j in zip(p,s)) for p in pts)) for s in CORNERS)
        orbits[key]+=prod(w for _,w in terms)
    values=labeled_vertex(genus,orbits,exe)
    return tuple(sum(w*x for w,x in zip(orbits.values(),col)) for col in zip(*values))


def four_line_probes(multiplicity):
    """Polarize a pure cubic on each of four lines whose momenta sum to zero."""
    # Each pure cubic can be represented by sinh of its momentum. Divide by
    # the factorial of each repeated exponential before grouping the fixtures.
    lines = (*UNIT,(-1,-1,-1))
    weights = Counter()
    for counts in product(range(multiplicity+1),repeat=4):
        reflected = tuple(multiplicity-n for n in counts)
        key = min(tuple(sorted(counts)),tuple(sorted(reflected)))
        weights[key] += prod(F((-1)**(multiplicity-n),
                               2**multiplicity*factorial(n)*factorial(multiplicity-n))
                             for n in counts)
    for counts, weight in sorted(weights.items()):
        normals = []
        for n,p in zip(counts,lines):
            normals.extend([p]*n+[tuple(-q for q in p)]*(multiplicity-n))
        yield counts, weight, tuple(normals)


def four_line_check():
    @cache
    def two(points):
        value, = labeled_vertex(2,[points],None)
        return value[0]

    assert sum(w*two(pts) for _,w,pts in four_line_probes(1)) == F(1,51840)
    probes = list(four_line_probes(2))
    assert len(probes) == 9
    # An eight-linear control with known polarization: distribute the eight
    # normal slots over two copies of the independently measured two-loop vertex.
    control = F(0)
    for _,weight,pts in probes:
        for left in combinations(range(8),4):
            control += weight*two(tuple(sorted(pts[i] for i in left)))*two(
                tuple(sorted(p for i,p in enumerate(pts) if i not in left)))
    assert control == F(1,51840**2)
    print("Four-loop normal probe weights (times 4096):",
          [(counts,4096*w) for counts,w,_ in probes],flush=True)


def elliptic_check():
    import sympy as sp
    x,y,z,b=sp.symbols('x y z b')
    two,three=hessian_pencil((x**3+y**3+z**3)/6+b*x*y*z,(x,y,z))
    s,t=sp.expand(51840*two),sp.expand(three)
    assert (s,t) == (b**4-b,1-20*b**3-8*b**6)
    assert sp.expand(t*t-64*s**3-(1+8*b**3)**3) == 0
    # Move the flex (1:-1:0) and its tangent to Weierstrass infinity.
    k=-2*b
    A,C,D,E=4*(k**3-1)/3,-3*k*k,2*k,-sp.Rational(1,3)
    xp=A*x+C/3
    assert sp.expand(A*A*(A*x**3+C*x*x+D*x+E)-xp**3
                     +16*s*xp/3+16*t/27) == 0
    # Test a nongeneric singular cubic, and generic integral shears away from Hesse form.
    assert hessian_pencil(x*x*y,(x,y,z)) == (0,0)
    cubic = x**3+3*x*y*z+2*x*y*y+7*z**3-y*z*z+2*x*x*y
    transformed=cubic.subs({x:x+y,y:2*y+z,z:x-2*z},simultaneous=True)
    det=sp.Matrix([[1,1,0],[0,2,1],[1,0,-2]]).det()
    original=hessian_pencil(cubic,(x,y,z))
    assert hessian_pencil(transformed,(x,y,z)) == (det**4*original[0],det**6*original[1])
    # Dehomogenize the first two doubled Hesse gradients at z=1.
    # Multiplication by 1+2bxy on the quotient with basis 1,x,y,xy has this determinant.
    sym=(x,y)
    ideal=sp.groebner([x*x+2*b*y,y*y+2*b*x],*sym,order='grlex',domain=sp.QQ.frac_field(b))
    basis=(1,x,y,x*y)
    cols=[sp.Poly(ideal.reduce((1+2*b*x*y)*m)[1],*sym) for m in basis]
    multiplication=sp.Matrix([[col.coeff_monomial(m) for col in cols] for m in basis])
    assert sp.expand(multiplication.det()-t*t+64*s**3) == 0
    P,R=sp.symbols('P R')
    result=hessian_pencil(y*y*z-x**3-P*x*z*z-R*z**3,(x,y,z))
    assert tuple(map(sp.expand,(51840*result[0],result[1]))) == (-48*P,-6912*R)
    # The fourth normal line is the negative of the sum of the first three.
    frame=hessian_pencil((x**3+y**3+z**3-b*(x+y+z)**3)/6,(x,y,z))
    assert tuple(map(sp.expand,(51840*frame[0],frame[1]))) == (b,1-6*b-3*b**2)


def jacobi_graphs():
    from bulk_eight import SHAPES, potentials
    from bulk_three_hexagon import potentials as hexagon
    # [[p,q],r]+cyclic on two Hamiltonian vectors p,q and a scalar r.
    p,q=HARD
    r=(1,2,0,0,0)
    def project(phi):
        return 2*(phi[0][1]-phi[-1][1])-3*(phi[1][0]-phi[-1][0])
    for power,orbits in cubic_orbits(2):
        shapes = [F(0)]*len(SHAPES)
        for normals,weight in orbits.items():
            inputs = [p,q,r]+[(0,0,*n) for n in normals]
            points = inputs+[tuple(-sum(c) for c in zip(*inputs))]
            for i,shape in enumerate(SHAPES):
                shapes[i] += weight*project(potentials(points,shapes=(shape,)))
        print('Two-loop potential by theta shape at b power',power,shapes,flush=True)
        sign={1:-1,4:1}[power]
        assert shapes == [sign*x for x in (F(41,90720),F(11,11340),F(1,2016),F(1,2520))]
        assert sum(shapes) == sign*F(1,432)
    # The Fermat chamber starts at three loops. Compute its original lower vertex insertion.
    power,orbits = next(cubic_orbits(3))
    assert power == 0
    value = F(0)
    for normals,weight in orbits.items():
        inputs = [p,q,r]+[(0,0,*n) for n in normals]
        points = inputs+[tuple(-sum(c) for c in zip(*inputs))]
        value += weight*project(hexagon(points))
    print('Three-loop Fermat hexagon projection',value,flush=True)
    assert value == -F(1,580608)
    return value


def analytic_binary(exe):
    from worldsheet_primary import fourfold_primary
    # Normals ordered (x,y,z). The z²/2 marks remove one target coordinate.
    for g in (2,3):
        total=[F(0),F(0)]
        for n in range(1,g):
            plane = [(0,0,1,0)]*n+[(0,0,0,1)]*(g-n)+[p[:4] for p in HARD]
            root = tuple(-sum(c) for c in zip(*plane))
            lower = fourfold_primary(g,[*plane,root])
            points = tuple([UNIT[2]]*g+[UNIT[0]]*n+[UNIT[1]]*(g-n))
            lifted, = labeled_vertex(g,[points],exe)
            assert tuple(x/factorial(g) for x in lifted) == lower
            for i,value in enumerate(lower):
                total[i] += value/(factorial(n)*factorial(g-n))
        print('Analytic primary binary',g,total,flush=True)
        assert tuple(total) == ((F(1,62208),F(0)) if g==2 else
                               (F(49188779,16295634432000),-F(7783,3143448000)))

if __name__ == '__main__':
    import sympy as sp
    x,y,z=sp.symbols('x y z')
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(ROOT/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        assert hesse(2,exe) == {1:(-F(1,51840),F(0)),4:(F(1,51840),F(0))}
        assert hesse(3,exe) == {r:tuple(k*x for x in FERMAT) for r,k in ((0,1),(3,-20),(6,-8))}
        assert weierstrass_graph(2,exe) == (planar_coefficients(y*y*z-x**3-x*z*z,(x,y,z))[0],0) == (-F(1,1080),0)
        assert weierstrass_graph(3,exe) == planar_coefficients(y*y*z-x**3-z**3,(x,y,z))[1]
        analytic_binary(exe)
        hard_normal_shears(exe)
    elliptic_check()
    four_line_check()
    jacobi_graphs()
