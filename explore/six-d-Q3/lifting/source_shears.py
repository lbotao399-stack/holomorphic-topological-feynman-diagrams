"""Classical mixed source generators in the repository's component jet convention."""
from fractions import Fraction
from itertools import product

from type_d import shear
from ops import L, Theory
from salg import EPS, add, gen, lderiv, letters_of, mul


def opposite(poly, flavor, axis):
    # Sign gives the standard action of the two odd generators on jet coefficients.
    result = {}
    for letter in letters_of(poly):
        kind, color, *powers = letter
        n = powers[axis]
        if n == 0:
            continue
        if kind == 1:
            target, sign = 2 + flavor, -1
        elif (kind, flavor) == (4, 1):
            target, sign = 2, 1
        elif (kind, flavor) == (3, 2):
            target, sign = 2, -1
        else:
            continue
        powers[axis] -= 1
        result = add(result, mul(gen((target,color,*powers)), lderiv(poly,letter)), sign*n)
    return result


def check():
    theory = Theory(3,{key:Fraction(value) for key,value in EPS.items()})
    count = 0
    for kind,flavor,axis in product(('c','b','g1','g2'),(1,2),range(3)):
        powers = (1,2,1)
        word = L(kind,0,powers)
        expected = powers[axis] + (kind == 'b' or kind == 'g'+str(flavor))
        anticommutator = add(opposite(shear(word,flavor,axis),flavor,axis),
                             shear(opposite(word,flavor,axis),flavor,axis))
        assert not add(anticommutator,word,-expected)
        probe = L(kind,0,tuple(int(k==axis) for k in range(3)))
        assert not add(theory.Q0(opposite(probe,flavor,axis)),
                       opposite(theory.Q0(probe),flavor,axis))
        count += 2
    return count


if __name__ == '__main__':
    print("Mixed Cartan and classical Ward relations:",check())
