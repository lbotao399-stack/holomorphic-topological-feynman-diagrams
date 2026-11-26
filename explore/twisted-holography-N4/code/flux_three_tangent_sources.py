"""Maximal normal pole of the primary D3 contact with three tangent-dependent scalars."""
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import factorial, prod
from pathlib import Path
from random import Random
from subprocess import run
from tempfile import TemporaryDirectory

from flux_primary import BULK, IIB, quartics, orbit_integers
from worldsheet_soft_rank import sharp_soft
from worldsheet_primary import fourfold_primary
from worldsheet_fourfold_positivity import fourfold_coefficients
from flint import fmpq_mat, fmpz_mat


def momentum_frame(genus, shifts, normal=(0,0,1), tangent=(1,1), third=(1,1)):
    # Coordinates dw1,dw2,dz1,dz2,dw3. Tangent exponents are z1,z2,z1+z2.
    d,e = tangent
    rows = [(1,0,0,0,0)]*(genus-1)+[(0,1,0,0,0)]*(genus-1)
    rows += [(i,j,x,y,0) for (i,j),(x,y) in zip(shifts,((d,0),(0,e),(third[0]*d,third[1]*e)))]
    rows += [(normal[0],normal[1],0,0,normal[2])]
    return [*rows,tuple(-sum(c) for c in zip(*rows))]


def reduced_traces(genus, frames, exe):
    if genus==2:
        return [(sum(x*y for x,y in zip(BULK,quartics(frame))),) for frame in frames]
    assert genus==3
    data = json.loads((IIB/'operator'/'three_loop_coefficients.json').read_text())
    return [tuple(sum(x*Fraction(y) for x,y in zip(row,coeff)) for coeff in
                  zip(*data['rational_and_zeta'])) for row in orbit_integers(exe, frames)]


def project_leading(genus, exe, third):
    grid = (-1,0,1)
    points = list(product(grid,repeat=2))
    frames = [momentum_frame(genus,((0,0),(0,0),(a,b)),third=third) for a,b in points]
    probes = [(((-2,4),(3,1),(3,2)),(4,-1,2),(2,-1)),
              (((4,1),(0,3),(-1,1)),(1,2,-1),(1,2))]
    frames.extend(momentum_frame(genus,*args,third=third) for args in probes)
    traces = reduced_traces(genus,frames,exe)
    # Every slice is at most quadratic. Interpolate both independent reduced channels.
    def basis(x):
        return (Fraction(x*(x-1),2), 1-x*x, Fraction(x*(x+1),2))
    leading = (Fraction(1,2),-1,Fraction(1,2))
    def evaluate(x,y):
        weights = [u*v for u,v in product(x,y)]
        return tuple(sum(w*c for w,c in zip(weights,col)) for col in zip(*traces[:9]))
    for (shifts,q,te),expected in zip(probes,traces[9:]):
        a,b = (z-third[0]*x-third[1]*y for x,y,z in zip(*shifts))
        scale = (te[0]*te[1]*q[2])**(2*genus)
        assert tuple(scale*c for c in evaluate(basis(a),basis(b)))==expected
    # In shifts (s,0),(0,t),(0,u), take [s²tu] from (-A s,u-B t).
    return tuple(-2*third[0]**2*third[1]*c for c in evaluate(leading,leading))


def repeated_soft_fourfold(genus, direction):
    count = genus-2
    frame = [tuple(int(i==j) for j in range(4)) for i in range(4)]
    def evaluate(t):
        soft = [tuple(t*x for x in direction)]*count
        rows = frame+soft
        return fourfold_primary(genus,rows+[tuple(-sum(col) for col in zip(*rows))])
    if not count:
        return evaluate(0)
    # The added waves and balancing shift are all parallel; each minor is affine in t.
    nodes = range(-genus,genus+1)
    powers = fmpq_mat([[x**j for j in range(2*genus+1)] for x in nodes])
    coefficients = powers.solve(fmpq_mat([[str(y) for y in evaluate(x)] for x in nodes])).tolist()
    assert all(not x for row in coefficients[:2*count] for x in row)
    return tuple(Fraction(str(v)) for v in coefficients[2*count])


def rank_at_infinity(genus):
    rows = [(0,0,1,0,0),(0,0,0,1,0),(2,3,1,1,0),(0,0,0,0,1)]
    rows += [tuple(-sum(c) for c in zip(*rows))]
    soft = [(1,0,0,0,0)]*(genus-1)+[(0,1,0,0,0)]*(genus-1)
    assert sharp_soft(genus,rows,soft)


def genus_four_weights(direction):
    """Resolve the repeated-soft quartic on the two simplex rays e1 and e1+e2."""
    mean = Fraction(sum(direction),5)
    centered = [x-mean for x in direction]
    second = sum(x*x for x in centered)
    fourth = sum(x**4 for x in centered)
    return (30*fourth-7*second**2)/8,(13*second**2-20*fourth)/12


def separate_second_insertion():
    """A quadratically soft orbit vanishing on the whole triple-collinear locus."""
    seed = ((0,1,2,3),(0,4,5,6),(0,1,2,4),(0,1,2,4))
    orbit = {tuple(sorted(tuple(sorted(order[i] for i in minor)) for minor in seed))
             for order in permutations(range(7))}
    assert len(orbit) == 1260
    hard = [tuple(int(i==j) for j in range(4)) for i in range(4)]+[(-1,)*4]
    leading = [word for word in orbit if all(sum(s in minor for minor in word)==1 for s in (5,6))]
    assert len(leading) == 180
    for v in ((1,0,0,0),(1,1,0,0),(1,-1,-1,0),(2,-1,4,3)):
        rows = [*hard,v,v]
        minors = {m:int(fmpz_mat([rows[i] for i in m]).det())**2 for m in combinations(range(7),4)}
        coefficient = sum(prod(minors[m] for m in word) for word in leading)
        assert coefficient == 48*genus_four_weights([*v,0])[1]
    print('Triple-collinear orbit: zero on the axial ray, 48 on the mixed ray.')


def check_five_mark_projection():
    """Integrate the two longitudinal soft columns by identity covariance at genus four."""
    rng = Random(29411)
    genus = 4
    for direction in ((1,0,0,0),(1,1,0,0)):
        # Four Abel differences from hard marks 1,...,4 to the balancing mark 5.
        hard = [fmpz_mat([[rng.randrange(-5,6) for _ in range(14)] for _ in range(genus)])
                for _ in range(4)]
        # In the mixed case use normal coordinates (1 minus 2, 3, 4)
        # and longitudinal coordinate 2. The axial case has longitudinal coordinate 1.
        mixed = direction[1]
        normal = [hard[0]-hard[1] if mixed else hard[1],hard[2],hard[3]]
        longitudinal = hard[mixed]
        projected = sum(int(fmpz_mat([*normal[0].tolist(),*normal[1].tolist(),*normal[2].tolist(),
                                      longitudinal.tolist()[i],longitudinal.tolist()[j]]).det())**2
                        for i,j in combinations(range(genus),2))
        direct = 0
        for i,j in product(range(genus),repeat=2):
            rows = [row+[direction[k]*int(h==i),direction[k]*int(h==j)]
                    for k,block in enumerate(hard) for h,row in enumerate(block.tolist())]
            direct += int(fmpz_mat(rows).det())**2
        assert direct == 2*projected and projected>0


def check_genus_four():
    separate_second_insertion()
    check_five_mark_projection()
    assert genus_four_weights([1,0,0,0,0])==(1,0)
    assert genus_four_weights([1,1,0,0,0])==(0,1)
    assert genus_four_weights([1,-1,0,0,0])==(4,1)
    assert genus_four_weights([1,-1,-1,0,0])==(4,Fraction(11,3))
    for d in product(range(-2,3),repeat=4):
        original = genus_four_weights([*d,0])
        assert original == genus_four_weights([3]+[x+3 for x in d])
        assert min(original)>=0
    # Determine the full degree-eight collinear restriction before taking [t^4].
    nodes = range(-4,5)
    rows = [[str(c) for c in fourfold_coefficients((1,t,t))] for t in nodes]
    coefficients = fmpq_mat([[t**i for i in range(9)] for t in nodes]).solve(fmpq_mat(rows))
    assert all(coefficients[i,j]==0 for i in range(4) for j in range(3))
    assert tuple(72*Fraction(str(coefficients[4,j])) for j in range(3)) == (
        -Fraction(11,16),-36,-36)


def check_complex_tangent_directions():
    import sympy as sp

    x, y, z = sp.symbols('x y z')
    mean = (1+x+y)/5
    centered = [entry-mean for entry in (1,x,y,0,0)]
    second, fourth = [sum(entry**degree for entry in centered) for degree in (2,4)]
    # The two quartic directions have a common zero exactly when p2=p4=0.
    assert sp.Matrix([[sp.Rational(30,8),-sp.Rational(7,8)],
                      [-sp.Rational(20,12),sp.Rational(13,12)]]).det() == sp.Rational(125,48)
    eliminant = sp.Poly(sp.resultant(second,fourth,y),x).monic()
    assert eliminant.degree() == 8 and eliminant.gcd(eliminant.diff()).degree() == 0
    assert (eliminant % sp.Poly(x**2+x+1,x)).is_zero
    assert sp.factor(2*(z+1)**5-5*(z+1)**2+3) == z**2*(2*z**3+10*z**2+20*z+15)
    for moment in (second,fourth):
        assert sp.rem(sp.expand(moment.subs(y,x**2)),x**2+x+1,x) == 0
    # On A+B=1 the three nonzero simplex coordinates sum to zero.
    assert sp.expand((2*fourth-second**2).subs(y,-1-x)) == 0
    assert sp.expand(second.subs(y,-1-x)-2*(1+x+x*x)) == 0
    for a in (Fraction(1,2),Fraction(1,3),Fraction(2),Fraction(-3,2)):
        q = a*a-a+1
        assert genus_four_weights([1,-a,a-1,0,0]) == (4*q*q,q*q)

    # Interpolate the two separate genus-three period channels from the trace.
    points = [(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1)]
    basis = sp.Matrix([[1,x,y,x*x,x*y,y*y]])
    matrix = sp.Matrix.vstack(*(basis.subs({x:i,y:j}) for i,j in points))
    readings = sp.Matrix([repeated_soft_fourfold(3,(1,i,j,0)) for i,j in points])
    fitted = basis * matrix.inv() * readings
    axis = repeated_soft_fourfold(3,(1,0,0,0))
    for observed,reference in zip(fitted,axis):
        assert sp.expand(observed - sp.Rational(5,4)*second*reference) == 0
    assert tuple(fitted.subs({x:-2,y:3})) == repeated_soft_fourfold(3,(1,-2,3,0))
    assert repeated_soft_fourfold(2,(1,0,0,0)) == (Fraction(1,62208),0)
    print('Complex D3 tangent momenta: eight common zeros; both genus-three channels vanish.')


def check():
    check_genus_four()
    check_complex_tangent_directions()
    for n in range(2,10):
        rank_at_infinity(n)
    with TemporaryDirectory() as tmp:
        exe = Path(tmp)/'nine'
        run(['c++','-std=c++17','-O3',str(IIB/'operator'/'nine_sextics.cpp'),'-o',str(exe)],check=True)
        for genus in (2,3):
            for a,b in ((1,1),(-2,1),(1,3),(0,2)):
                coefficient = project_leading(genus,exe,(a,b))
                lower = repeated_soft_fourfold(genus,(1,-a,-b,0))
                multiplier = 4*a*a*b*factorial(genus-1)**2//factorial(genus-2)
                assert coefficient == tuple(multiplier*x for x in lower[:len(coefficient)])
                if (a,b)==(1,1):
                    expected = ((Fraction(1,15552),) if genus==2 else
                                (Fraction(740466719,7565830272000),-Fraction(469153,5837832000)))
                    assert coefficient == expected
                    print("genus",genus,"maximal contact before divergence and flux derivatives:",*coefficient)
                if genus==3:
                    axis = repeated_soft_fourfold(genus,(1,0,0,0))
                    ratio = Fraction(5*(1+a*a+b*b)-(1-a-b)**2,4)
                    assert lower == tuple(ratio*x for x in axis)


if __name__=='__main__':
    check()
