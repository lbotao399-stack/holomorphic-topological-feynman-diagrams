"""Closed point insertions from unrestricted BCOV sources on minimal fivefold fields."""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial, prod

from elliptic_descendants import vertex as elliptic_vertex
from genus_one_graph_ward import graph_vertex
from quadratic_loop_selection import regular_elliptic
from scalar_ancestor_response import scalar_power
from two_loop_single_descendant import single_vertex


def one_axis(genus):
    """Unlabeled Fourier weights of (cosh z - 1)^genus/genus!, modulo z -> -z."""
    states={}
    for minus in range(genus+1):
        for plus in range(minus,genus-minus+1):
            zeros=genus-minus-plus
            if not (minus or plus):
                continue  # The fivefold momentum rank drops on this axis.
            weight=Fraction((-1)**zeros,2**(plus+minus)*factorial(minus)*factorial(zeros)*factorial(plus))
            state=(-1,)*minus+(0,)*zeros+(1,)*plus
            states[state]=weight*(1 if minus==plus else 2)
    return states


def ward_sum(genus, momenta):
    """Fivefold graph calculation after string at the zero-charge primary and constant fields."""
    states=one_axis(genus)
    total=Fraction()
    active=0
    for orbit in combinations_with_replacement(states,4):
        waves=[(p,0,0,0,0) for p in momenta]
        for direction,entries in enumerate(orbit,1):
            waves.extend(tuple(charge*int(i==direction) for i in range(5)) for charge in entries)
        root=tuple(-sum(column) for column in zip(*waves))
        if genus==1:
            value=graph_vertex((root,*waves),(len(waves)-4,)+(0,)*len(waves))
        elif genus==2:
            value=single_vertex((root,*waves),0,full_source=True)
        else:
            raise ValueError('The fivefold Ward replay covers genera one and two.')
        multiplicity=24//prod(factorial(v) for v in Counter(orbit).values())
        total+=multiplicity*prod(states[p] for p in orbit)*value
        active+=int(bool(value))
    return total,active


def one_mark(genus,momenta):
    series=[Fraction(1)]+[Fraction() for _ in range(genus)]
    for p in momenta:
        factor=[Fraction(p**(2*j),4**j*factorial(2*j+1)) for j in range(genus+1)]
        series=[sum(series[n-j]*factor[j] for j in range(n+1)) for n in range(genus+1)]
    return series[genus]


def mixed_two_loop_jets():
    """Reconstruct the mixed fourth-derivative coefficients from fivefold graphs."""
    single = ward_sum(2, (1, 0))[0]
    aligned = ward_sum(2, (1, 1))[0]
    opposite = ward_sum(2, (1, -1))[0]
    # Symmetry gives one coefficient each for the 4+0, 3+1 and 2+2 splittings.
    odd = (aligned - opposite) / 4
    paired = (aligned + opposite) / 2 - 2 * single
    assert (single, odd, paired) == (Fraction(1, 1920), 0, Fraction(1, 576))
    assert opposite == Fraction(1, 360)
    print('Reconstructed two-input coefficients at 4+0, 3+1, 2+2:', single, odd, paired, flush=True)


def protected_level_one():
    """One fixed closed level-one source on the level-three scalar profile."""
    for genus in range(1, 4):
        charges = (0, 1, -1, 0) + (0,) * genus
        levels = (1, 0, 0, 0) + (3,) * genus
        expected = Fraction(2 * factorial(3*genus+1),
                            6**genus * factorial(genus) * factorial(2*genus+2))
        assert scalar_power(genus, genus, charges, levels) / factorial(genus) == expected
        assert expected * 6**genus == Fraction(2 * comb(3*genus+2, genus), 3*genus+2)

    for genus, source, momenta, hard in (
            (1, 1, (1, 2, -3, 0), (4,)),
            (1, 2, (1, 2, -3, 0), (2, 2)),
            (2, 2, (1, -1, 0), (4,)),
            (2, 2, (1, -1, 0), (2, 3)),
            (2, 4, (1, -1, 0), (1, 2))):
        labels = (0,) + momenta + (0,) * len(hard)
        levels = (source,) + (0,) * len(momenta) + hard
        degree = 2*genus + len(momenta) + len(hard) - 2
        assert sum(levels) == degree
        labeled = (one_mark(genus, momenta) * factorial(degree) /
                   (factorial(source) * prod(factorial(j) for j in hard)))
        assert scalar_power(genus, genus, labels, levels) == labeled

    # Each of the four independent normal axes supplies one donor. Reflection
    # identifies the two nonzero cosine charges, and a zero drops the target rank.
    normal = [tuple(int(i == j) for i in range(5)) for j in range(1, 5)]
    zero = (0,) * 5
    tangent = (1, 0, 0, 0, 0)
    recoil = (0, -1, -1, -1, -1)
    terms = []
    for attached in range(4):
        transfer = tuple(sum(p[i] for p in normal[:attached]) for i in range(5))
        waves = (recoil, transfer, tangent, tuple(-t for t in tangent), zero, *normal[attached:])
        levels = (1, 3-attached) + (0,) * (len(waves)-2)
        direct = regular_elliptic(waves, levels)
        assert direct == elliptic_vertex(waves, levels)
        terms.append((-1)**attached * comb(4, attached) * direct)
    assert sum(terms) == Fraction(1, 3)
    print('Closed level-one insertion, one-loop fivefold fusion terms:', *terms, flush=True)


def flux_descendant():
    """Ward graphs on B=z1*w1+z2*w2 and a wave with general normal flux momentum."""
    from flint import fmpz_mat

    # Coordinates are (w3,z1,w1,z2,w2). Four signed cosines give the quadratic
    # jet of B, with Hessian weights (+1/2,-1/2,+1/2,-1/2).
    normals = ((0,1,1,0,0), (0,1,-1,0,0), (0,0,0,1,1), (0,0,0,1,-1))
    flux = (2,0,3,0,-2)

    def amplitude(orbit, genus, moving_flux):
        waves = [moving_flux]
        for entries, axis in zip(orbit, normals):
            waves.extend(tuple(x*n for n in axis) for x in entries)
        root = tuple(-sum(row) for row in zip(*waves))
        if genus == 1:
            return graph_vertex((root,*waves), (1,) + (0,)*len(waves))
        return single_vertex((root,*waves), full_source=True)

    elliptic = sum(amplitude(tuple((s,) for s in signs),1,flux)
                   for signs in product((-1,1), repeat=4)) / 16**2
    assert elliptic == Fraction(flux[0]**2, 24)

    states = one_axis(2)
    result = Fraction()
    for orbit in combinations_with_replacement(states,4):
        aligned = amplitude(orbit, 2, (1,0,0,0,0))
        # Transfer two normal components onto the sole w3-dependent input.
        # Its squared axial charge occurs twice in the worldsheet symbol.
        assert amplitude(orbit, 2, flux) == flux[0]**4 * aligned
        copies = 24 // prod(factorial(n) for n in Counter(orbit).values())
        result += copies * prod(states[s] for s in orbit) * aligned / 16**2
    assert result == Fraction(1, 1920)
    normal_hessian = ((1,2,-1),(2,4,0),(-1,0,3))
    for left, right in (((1,2,-1),(0,3,1)), ((2,1,0),(4,2,0))):
        hessian = [(0,0,*left), (0,0,*right)]
        hessian.extend((left[i],right[i],*row) for i,row in enumerate(normal_hessian))
        direction = tuple(left[(i+1)%3] * right[(i+2)%3] - left[(i+2)%3] * right[(i+1)%3]
                          for i in range(3))
        for i,j in product(range(3), repeat=2):
            minor = fmpz_mat([[x for col,x in enumerate(row) if col != i+2]
                              for row_index,row in enumerate(hessian) if row_index != j+2])
            assert (-1)**(i+j) * int(minor.det()) == direction[i] * direction[j]
    for genus in range(1, 9):
        vertex = Fraction(1, 4**genus * factorial(2*genus+1))
        residue = vertex * factorial(2*genus)
        assert residue == Fraction(1, 4**genus * (2*genus+1))
        if genus <= 3:
            assert scalar_power(genus, genus, (-1,1), (2*genus-1,0)) == vertex
    print('One-flux normal three-form, weighted genus-one and genus-two residues:',
          elliptic * factorial(2) / flux[0]**2, result * factorial(4), flush=True)


def check_integrated_sources():
    from elliptic_descendants import cut_recursion

    basis = tuple(tuple(int(i == j) for j in range(5)) for i in range(5))
    zero = (0,)*5
    hard = (*basis, (-1,)*5)
    levels = (2,) + (0,)*6
    assert elliptic_vertex((zero,*hard), levels) == cut_recursion((zero,*hard), levels) == Fraction(1,4)
    first_descendant = (1,) + (0,)*5
    assert elliptic_vertex(hard, first_descendant) == cut_recursion(hard, first_descendant) == Fraction(1,24)
    for source_level, expected in ((0,Fraction(1,24)), (1,Fraction(1,4))):
        degrees = (source_level, 2-source_level) + (0,)*5
        assert elliptic_vertex((zero,*hard), degrees) == cut_recursion((zero,*hard), degrees) == expected
    # Interpolate the Ward rows to a constant source. Each of the two
    # determinant factors is affine along this common direction of recoil.
    samples = []
    for step in (1,2,3):
        source = (0,step,-step,0,0)
        recoil = tuple(x-y for x,y in zip(hard[0],source))
        samples.append(graph_vertex((source,recoil,*hard[1:]), levels))
    assert samples == [Fraction(1,3), Fraction(7,12), 1]
    assert 3*samples[0] - 3*samples[1] + samples[2] == Fraction(1,4)
    for p, t in ((-2,(3,1,-1,2)), (3,(2,1,0,-1))):
        waves = tuple(tuple(a*x for x in row) for a,row in zip((p,*t),basis))
        balance = tuple(-sum(column) for column in zip(*waves))
        assert elliptic_vertex((zero,*waves,balance), levels) == Fraction((p*prod(t))**2,4)
    # One nonzero all-primary vertex from the fivefold genus-two reconstruction.
    hard = (*basis, (1,1,0,0,0), (-2,-2,-1,-1,-1))
    primary = single_vertex(hard)
    assert primary == Fraction(17,155520)
    assert single_vertex((zero,*hard)) == 9*primary
    assert single_vertex((zero,*hard), marked=1) == primary
    print("Integrated sources: elliptic partitions, cuts and hexagon Ward rows; genus-two string and dilaton pass.")


def check(full_source_test=False):
    check_integrated_sources()
    for genus in range(1,6):
        degrees=(2*genus,0,0)
        expected=Fraction(1,4**genus*factorial(2*genus+1))
        assert scalar_power(genus,genus,(-1,1,0),degrees)==expected
    probes=[(1,(-2,)),(1,(1,2,3,4,5)),(2,(1,)),(2,(1,2))]
    if full_source_test:
        probes.append((2,(1,1,1)))
    for genus,charges in probes:
        value,count=ward_sum(genus,charges)
        assert value==one_mark(genus,charges)
        print('Fivefold Ward sum, genus, tangent momenta, nonzero Fourier orbits:',value,genus,charges,count,flush=True)
    if full_source_test:
        mixed_two_loop_jets()
        protected_level_one()
        flux_descendant()
    # Test a level genuinely outside the minimal source range before applying
    # any string equations. Zero Fourier terms on a normal axis have deficient rank.
    basis=tuple(tuple(int(i==j) for i in range(5)) for j in range(5))
    root=tuple(-1 for _ in range(5))
    zero=(0,)*5
    for r in (2,5,6):
        waves=(root,*basis,*([zero]*(r-1)))
        v=elliptic_vertex(waves,(r,)+(0,)*(len(waves)-1))
        m=r-2
        assert v/factorial(m)==Fraction(1,24*factorial(m))
    print('Scalar stable graphs through genus five and unrestricted elliptic source levels 2, 5, 6 passed.')


if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser(description=__doc__)
    parser.add_argument('--full-source',action='store_true',help='also sum the level-five twelve-wave system and reconstruct the two-loop mixed jets')
    parser.add_argument('--integrated',action='store_true',help='run only the conserved integrated source checks')
    args=parser.parse_args()
    if args.integrated:
        check_integrated_sources()
    else:
        check(args.full_source)
