"""Positive-time angular-test pentagon: exact bridge and spherical coefficients."""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, prod

import sympy as sp


def check():
    # Root the cycle at vertex 0. Factor each contracted five-dimensional
    # volume into four forms in the relative odd coordinates and take their
    # 20 by 20 determinant; this works for the heat and propagator odd ends.
    def odd_cycle(currents):
        rows = []
        factor = sp.S.One
        for vertex, current in enumerate(currents):
            pivot = next(i for i in range(5) if current[i])
            factor *= sp.Rational((-1)**pivot, current[pivot]**3)
            for i in range(5):
                if i == pivot:
                    continue
                row = [0]*20
                for site, sign in ((vertex,1),((vertex+1)%5,-1)):
                    if site:
                        row[5*(site-1)+i] += sign*current[pivot]
                        row[5*(site-1)+pivot] -= sign*current[i]
                rows.append(row)
        return factor*sp.Matrix(rows).det()

    ordinary = [sp.eye(5)[:,i]+sp.eye(5)[:,(i+1)%5] for i in range(4)]
    ordinary_end = -sum(ordinary,sp.zeros(5,1))
    assert odd_cycle(ordinary+[ordinary_end]) == 0
    transverse_end = sp.eye(5)[:,4]
    nonzero = odd_cycle(ordinary+[transverse_end])
    assert abs(nonzero) == 1
    assert odd_cycle(ordinary+[transverse_end+ordinary_end]) == nonzero

    # At genus one the distinguished level-one quartic fixes the root. A
    # ring with a labeled orientation has five choices for its partner and
    # 4! orders on the remaining vertices. Each line has two orientations.
    root_counts = {"zero scalar": 0, "horizontal": 0, "annulus": 0, "interior": 0}
    for partner, edge, orientation in product(range(5), range(5), (0,1)):
        angular = (edge+orientation) % 5
        annular = (edge+1-orientation) % 5
        category = ("zero scalar" if partner < 4 else "horizontal" if angular == 0
                    else "annulus" if annular == 0 else "interior")
        root_counts[category] += factorial(4)
    assert list(root_counts.values()) == [960,48,48,144]

    # Raise the source to level two and add exactly one axial scalar. Delete
    # its leaf at the root and retain the complete labeled cyclic ordering.
    mixed = tuple(range(4))
    flux, axial = 4, 5
    matched = set()
    zero_roots = 0
    for pair in combinations((*mixed,flux,axial),2):
        for order in permutations(i for i in (*mixed,flux,axial) if i not in pair):
            for edge, orientation in product(range(5),(0,1)):
                if any(i in mixed for i in pair):
                    zero_roots += 1
                else:
                    assert pair == (flux,axial)
                    matched.add((tuple(order),edge,orientation))
    assert zero_roots == 3360
    assert matched == set(product(permutations(mixed),range(5),(0,1)))

    # Genus-one Euler derivatives E_1 and E_2 on mixed and axial Hessians.
    # Dolbeault and polyvector factors are spectators in this scalar check.
    c, w1, w2, w3, z1, z2 = sp.symbols('c w1 w2 w3 z1 z2')
    coords = (w3,z1,w1,z2,w2)
    background = z1*w1 + z2*w2 + c*w3**2/2
    laurent = 1/(w1*w2*w3)
    hessian = sp.hessian(background, coords)
    first = sp.trace(hessian.adjugate() * sp.hessian(laurent, coords)) / 24
    second = (laurent*hessian.det()/24 + background*first).expand().coeff(c)
    assert sp.simplify(first*w3**2/laurent) == sp.Rational(1,12)
    assert sp.simplify(second/laurent) == sp.Rational(1,12)
    assert sp.simplify(second - w3**2*first/2 - laurent/24) == 0

    # Both annular weights are constant. Move the extra scalar to the flux
    # leaf of the quartic pentagon. Two independent Taylor variables make
    # the quintic and quartic contact coefficients directly comparable.
    extra, strength = sp.symbols('extra strength')
    mixed_background = z1*w1 + z2*w2
    for pole in (1, 2, 4):
        flux_mode = laurent / w3**(pole-1)
        profile = w3**(pole+1) + w1*w3 + w2**2
        primary = mixed_background + extra*profile + strength*flux_mode
        quintic = sp.expand(primary * sp.hessian(primary, coords).det() / 24)
        quintic = quintic.coeff(extra, 1).coeff(strength, 1)
        fused = mixed_background + strength*profile*flux_mode
        quartic = sp.expand(sp.hessian(fused, coords).det() / 24).coeff(strength, 1)
        difference = sp.expand(quintic-quartic)
        assert sp.simplify(difference + sp.diff(profile,w3)*sp.diff(flux_mode,w3)/12) == 0
        residue = difference.coeff(w1,-1).coeff(w2,-1).coeff(w3,-1)
        assert residue == sp.Rational(pole*(pole+1),12)

    # Smooth a spherical harmonic H_l(x) / |x|^(2 nu) on the angular
    # Gaussian arm. Compare a sphere moment followed by the radial gamma
    # integral with the beta moment from the radial Schwinger parameter.
    for ell in range(25):
        for nu, j in product(range(1,ell+3), range(9)):
            sphere = Fraction(2, factorial(j)*factorial(ell+j+2))
            radial = Fraction(factorial(ell+2-nu+j),2)
            beta = Fraction(factorial(nu-1)*factorial(ell+2-nu+j),
                            factorial(ell+2+j))
            assert sphere*radial == beta/(factorial(nu-1)*factorial(j))

    # Invert the grounded heat precision before taking differences of the
    # internal position covariances; the short edge has length 1/100000.
    for lengths in ((1, 2, 3), (Fraction(1,100000), 3, 5)):
        grounded = sp.zeros(2)
        for i, length in enumerate(lengths):
            direction = sp.Matrix([int(j == i+1)-int(j == i) for j in (1,2)])
            grounded += direction*direction.T/(4*sp.Rational(length))
        cov = sp.zeros(4)
        cov[1:3,1:3] = grounded.inv()
        for edge, vertex in product((1,2,3), range(4)):
            observed = (cov[edge,vertex]-cov[edge-1,vertex])/lengths[edge-1]
            elapsed = sp.Rational(sum(lengths[:vertex]))/sum(lengths)
            assert observed == 4*(int(edge <= vertex)-elapsed)

    # Schur complement of the ordinary path with a normal mass at the flux.
    for lengths in ((1, 2, 3, 4), (11, 2, 17, 3)):
        for site, mass in product(range(1, 4), (0, Fraction(4, 3), 53)):
            laplace = sp.zeros(5)
            for i, length in enumerate(lengths):
                weight = sp.Rational(1, 4*length)
                laplace[i,i] += weight
                laplace[i+1,i+1] += weight
                laplace[i,i+1] = laplace[i+1,i] = -weight
            laplace[site,site] += mass
            ends, interior = (0,4), (1,2,3)
            bridge = (laplace.extract(ends,ends)
                      - laplace.extract(ends,interior)*laplace.extract(interior,interior).inv()
                      * laplace.extract(interior,ends))
            alpha = sp.Rational(1, 4*sum(lengths[:site]))
            beta = sp.Rational(1, 4*sum(lengths[site:]))
            denom = alpha+beta+mass
            cross = alpha*beta/denom
            assert bridge == sp.Matrix([[cross+alpha*mass/denom, -cross],
                                        [-cross, cross+beta*mass/denom]])

    # Integrate the harmonic z_1^ell at the north pole of the real unit S^5.
    # The expansion of exp(u(z_1+bar z_1)) selects p=j and q=j+ell.
    for ell, j in product(range(30), range(18)):
        beta_moment = Fraction(2, (j+ell+1)*(j+ell+2))
        measured = beta_moment / (factorial(j)*factorial(j+ell))
        assert measured == Fraction(2, factorial(j)*factorial(j+ell+2))

    # Independent normal integration of f^n bar(f)^n / |f|^6. Strip pi^3.
    # A source derivative with respect to J supplies an extra bar(f).
    for index, powers in product(range(3), product(range(5), repeat=3)):
        if powers[index] == 0:
            continue
        order = sum(powers)
        h = 7
        angular = Fraction(2*prod(map(factorial, powers)), factorial(order+2))
        radial = Fraction(factorial(order-1), 2*h**order)
        anti = list(powers)
        anti[index] -= 1
        series = Fraction(factorial(order-1), h**order * factorial(order+2)
                          * prod(map(factorial, anti)))
        assert series*prod(map(factorial, powers+tuple(anti))) == angular*radial

    # Include both harmonic degrees, a bounded sphere-degree shift and a
    # possible zero cross precision. The identity involves no inverse precision.
    for a, k, shift, j in product(range(1,6), range(12), range(-4,5), range(5)):
        ell = a+k+shift
        if ell < 0:
            continue
        radius, annulus = Fraction(2,9), Fraction(3,2)
        for cross in (0, Fraction(31,7)):
            x, y = cross*radius**2, cross*annulus**2
            left = (radius/annulus)**(a-k)*(cross*radius*annulus)**(ell+2*j)
            right = ((radius/annulus)**(2*a)*(annulus/radius)**(-shift)
                     * y**(ell+j)*x**j)
            assert left == right
    print("All interior flux sites, 300 flux moments, 540 spherical coefficients and both-degree shifts pass.")


if __name__ == "__main__":
    check()
