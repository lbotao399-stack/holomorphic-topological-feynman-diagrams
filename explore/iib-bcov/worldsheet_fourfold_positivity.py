"""Full-rank genus-four Betti tests on the curve y^3=x^5-1."""
from fractions import Fraction
from itertools import combinations, product
from math import prod

from flint import fmpq_poly, nmod_mat

from worldsheet_primary import restricted_symbols


PRIME = 211
POINTS = ((13,51),(13,79),(13,81),(19,37),(19,78),(19,96),(20,45),(20,169),
          (20,208),(25,49),(25,53),(25,109),(30,45),(30,169),(32,10),(38,36))
HOLOMORPHIC_SHIFT = ((2,0,-1,-1),(0,1,-2,-1),(2,1,-2,-2),(-2,1,1,-2))


def determinant(strengths):
    def inv(a):
        return pow(a, -1, PRIME)
    roots_x = [a for a in range(PRIME) if pow(a,5,PRIME) == 1]
    roots_y = [b for b in range(PRIME) if pow(b,3,PRIME) == PRIME-1]
    # Mark three x-branch points, infinity, two other x-branch points and a y-branch point.
    slots = (4,1,0,8,2,3,7)
    momenta = [(a,0,0,0) for a in strengths]
    momenta += [(0,1,0,0),(0,0,1,0),(0,0,0,1),(-sum(strengths),-1,-1,-1)]
    rows = []
    for x,y in POINTS:
        assert (y**3-x**5+1) % PRIME == 0
        dy = 5*x**4*inv(3*y*y)
        holomorphic = (inv(y), inv(y*y), x*inv(y*y), x*x*inv(y*y))
        # Coefficients relative to dx. Each finite cusp differential has residues +1,-1.
        logarithmic = [inv(3*(x-a)) for a in roots_x]
        logarithmic += [dy*inv(5*(y-b)) for b in roots_y] + [0]
        differentials = [sum(a*logarithmic[i] for i,a in zip(slots, column))
                         + sum(a*h for a,h in zip(shift, holomorphic))
                         for column, shift in zip(zip(*momenta), HOLOMORPHIC_SHIFT)]
        rows.append([f*h % PRIME for f,h in product(holomorphic, differentials)])
    return int(nmod_mat(rows, PRIME).det())


def fourfold_coefficients(strengths):
    """Coefficients of K4, J4 and V(1,1,1) after taking the fourfold quotient."""
    u = (1,0)
    a,b,d = restricted_symbols([u]*4+[(0,s) for s in strengths])
    ambiguity = Fraction(22*a-5*b+18*d, 288)
    return tuple(-x/24 for x in (Fraction(b,576)+3*ambiguity/8,
                                  Fraction(8*a-b,9)+8*ambiguity, 24*ambiguity))


def check_double_normal_split():
    # Six line marks split to distinct slopes; the seventh planar mark is e^2.
    t = fmpq_poly([0,1])
    strengths = (1,-2,1,3,2,-1)
    slopes = (-3,-1,0,1,2,5)
    plane = [(a,a*s*t) for a,s in zip(strengths,slopes)] + [(0,1)]
    polynomials = restricted_symbols(plane)
    assert all(not p[i] for p in polynomials for i in range(4))
    a,b,d = (Fraction(str(p[4])) for p in polynomials)
    assert (a,b,d) == (931392,1862784,0)
    ambiguity = (22*a-5*b+18*d)/288
    direct = (-(b/576+3*ambiguity/8), -((8*a-b)/9+8*ambiguity), -24*ambiguity)

    pairs = combinations(range(6),2)
    matchings = [(i,j,k,l) for (i,j),(k,l) in combinations(pairs,2) if len({i,j,k,l})==4]
    differences = sum((slopes[i]-slopes[j])**2*(slopes[k]-slopes[l])**2
                      for i,j,k,l in matchings)
    assert len(matchings) == 45 and differences == 6468
    # V4(tu,tv,1) / t^4 at t=0 gives u²v² times these three coefficients,
    # from the axial and centered triple-strength components in check().
    norm = (Fraction(-11,1152), -Fraction(1,2), -Fraction(1,2))
    assert direct == tuple(2*prod(a*a for a in strengths)*differences*c for c in norm)
    assert direct == (-17787,-931392,-931392)
    # One first-order outlier; five other first jets coincide and split at t².
    transverse = (t,0,t*t,2*t*t,3*t*t,4*t*t)
    split = restricted_symbols([(1,h) for h in transverse]+[(0,1)])
    matched = sum(((transverse[i]-transverse[j])*(transverse[k]-transverse[l]))**2
                  for i,j,k,l in matchings)
    assert all(not p[h] for p in (*split,matched) for h in range(6))
    assert tuple(int(p[6]) for p in (*split,matched)) == (150,300,0,150)
    soft = restricted_symbols([(t,t*transverse[0])]+[(1,h) for h in transverse[1:]]+[(0,1)])
    assert all(not p[h] for p in soft for h in range(8))
    assert tuple(int(p[8]) for p in soft) == (150,300,0)


def check():
    check_double_normal_split()
    assert [determinant(s) for s in ((1,1,1),(1,1,-2))] == [172,129]
    axis = (0,0,1)
    centered = (-Fraction(11,32), -18, -22)
    assert fourfold_coefficients((1,1,1)) == axis
    assert fourfold_coefficients((1,1,-2)) == centered
    for strengths in product(range(-2,3), repeat=3):
        mean_part = Fraction(sum(strengths)**2, 9)
        contrast = (sum(a*a for a in strengths)-3*mean_part) / 24
        assert fourfold_coefficients(strengths) == tuple(
            prod(strengths)**2*(mean_part*a + contrast*c) for a,c in zip(axis,centered))
    print("Fourfold Betti multiplication determinants modulo 211:", 172, 129)
    print("Centered volume = -11 K4/32 - 18 J4 - 22 V(1,1,1); both volumes strictly positive.")


if __name__ == '__main__':
    check()
