"""Fixed quadratic background and constant scalar descendants in the flat fivefold."""
from fractions import Fraction
from math import factorial, prod

from descendant_background import coefficient
from genus_one_graph_ward import graph_vertex
from scalar_ancestor_response import homogeneous_dr, scalar_power


def response(genus, levels):
    """Heat-normalized coefficient for one tangent scalar, one constant tangent vector."""
    if genus < 0 or any(k not in (2, 3, 4) for k in levels):
        raise ValueError("Use a nonnegative genus and constant scalar levels two to four.")
    if sum(k-1 for k in levels) != 2*genus:
        return Fraction()
    denominator = factorial(2*genus+2) * prod(factorial(k) for k in levels)
    denominator *= prod(factorial(levels.count(k)) for k in set(levels))
    return Fraction((-1)**genus * 2*factorial(2*genus+len(levels)), denominator)


def regular_elliptic(momenta, levels):
    # Move all waves in a single target coordinate and balance the last.
    # Its two Betti periods give degree at most two in this deformation.
    stride=1+max(abs(p[0]) for p in momenta)
    def sample(step):
        step*=stride
        shifted = [(p[0]+step, *p[1:]) for p in momenta[:-1]]
        last = momenta[-1]
        shifted.append((last[0]-len(shifted)*step, *last[1:]))
        return graph_vertex(tuple(shifted), tuple(levels))
    return 6*sample(2)-8*sample(3)+3*sample(4)


def check():
    from math import comb
    from itertools import product
    from flint import fmpq_poly
    for genus, levels in ((1,(3,)), (1,(2,2)), (2,(3,3)), (2,(2,4)),
                          (3,(4,4)), (3,(3,3,3)), (4,(3,3,3,3))):
        charges = (1,-1,0) + (0,)*len(levels)
        degrees = (0,0,0) + levels
        repeats = prod(factorial(levels.count(k)) for k in set(levels))
        labeled = (-1)**genus * repeats * response(genus,levels)
        assert scalar_power(genus,genus,charges,degrees) == labeled
        if genus <= 3:
            assert Fraction(str(homogeneous_dr(genus,charges).to_dict()[degrees])) == labeled
    # The unit coefficient for four independent transverse exponentials tests
    # the polarized Hessian by the actual hexagon and its tree Ward terms.
    for p in (1,-2,3):
        action=coefficient(regular_elliptic, ((0,),(p,),(0,),(-p,)),
                           (3,0,0,0), (1,1,1,1))
        assert -p*action == -Fraction(p**3,12)
    print('Scalar target stable graphs through genus four, full marked formula through genus three,'
          ' and the genus-one fivefold Ward graphs agree.')
    ternary = [Fraction(comb(3*g,g),2*g+1) for g in range(17)]
    for g in range(1,len(ternary)):
        assert ternary[g] == sum(ternary[i]*ternary[j]*ternary[g-1-i-j]
                                for i in range(g) for j in range(g-i))
    for g, coeff in enumerate(ternary):
        square = sum(ternary[i]*ternary[g-i] for i in range(g+1))
        assert (3*coeff-square)/2 == coeff/(g+1) == (-6)**g*response(g,(3,)*g)
    assert response(3,(4,4)) == -Fraction(1,576)
    # Revert the classical scalar equation degree by degree, with all three
    # background amplitudes active, independently of the marked-point formula.
    for strengths in ((1,-2,3),(-3,1,2)):
        s=fmpq_poly([0,1])
        implicit=s
        for degree in range(2,14):
            rhs=s+sum((amplitude*implicit**k/factorial(k)
                       for k,amplitude in enumerate(strengths,2)),fmpq_poly())
            implicit += fmpq_poly([0]*degree+[rhs[degree]])
        for g in range(7):
            terms=Fraction()
            for counts in product(range(2*g+1),range(g+1),range(2*g//3+1)):
                if sum((k-1)*m for k,m in enumerate(counts,2))!=2*g:
                    continue
                ks=tuple(k for k,m in enumerate(counts,2) for _ in range(m))
                terms+=response(g,ks)*prod(a**m for a,m in zip(strengths,counts))
            assert terms==Fraction(str(implicit[2*g+1]))*(-1)**g/(g+1)
    print('Level three, heat units:',*(response(g,(3,)*g) for g in range(6)))


if __name__ == '__main__':
    check()
