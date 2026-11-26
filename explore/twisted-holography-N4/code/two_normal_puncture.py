"""Primary two-normal contact on D3 puncture jets and the regular brane unit."""

from fractions import Fraction
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import sys

import sympy as sp

from quadratic_d3_factorization import exterior_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "iib-bcov"))
from longitudinal_theta import primary_coefficient


def two_loop_vertex():
    # Two labeled hard scalars, one quartic normal jet, two quadratic jets on
    # the other normal axis, and one hard vector. The two equal donors cost 2!.
    # Homogeneity fixes four derivatives on the first normal and four on the
    # second. Double softness leaves exactly two on each second-normal donor.
    normal = [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (0, 0, 0, 1, 0)]
    powers = (4, 3, 2, 1, 1, 2)
    coefficient = Fraction(0)
    for values in product(*(range(degree+1) for degree in powers)):
        a, b, c, s, t, u = values
        hard = [(a, 0, 0, 0, s), (0, b, 0, 0, t)]
        vector = (c, c, 0, 0, u)
        result = primary_coefficient([*hard, *normal], vector, (1, -1, 0, 0, 0))/2
        # Check the whole momentum polynomial before extracting the jet.
        expected = Fraction((a-b)*(a*b*u-b*c*s-a*c*t)**4, 62208)
        assert result == expected
        weight = prod(Fraction((-1)**(d-n)*comb(d,n), factorial(d))
                      for d,n in zip(powers, values))
        coefficient += weight*result
    assert coefficient == Fraction(1, 5184)
    print("Two-loop reduced puncture scalar from the seven-field vertex:", coefficient)


def regular_brane():
    z1, z2, w1, w2, w3 = sp.symbols("z1 z2 w1 w2 w3")
    normal = (w1, w2, w3)
    wedge, contraction = exterior_matrices()
    potential = (w1**2+w2**2)/2 + w1**4/24
    differential = sum((w*h for w,h in zip(normal, contraction)), sp.zeros(8))
    differential += (w1/2+w1**3/24)*wedge[0] + w2*wedge[1]/2
    assert (differential**2-potential*sp.eye(8)).applyfunc(sp.expand) == sp.zeros(8)
    assert differential.subs(dict.fromkeys(normal, 0)) == sp.zeros(8)
    for genus in range(2, 9):
        degree = 2*genus
        for coefficient in (z1**degree/factorial(degree), z2**(degree-1)/factorial(degree-1)):
            primitive = coefficient*wedge[2]
            assert differential*primitive+primitive*differential == w3*coefficient*sp.eye(8)
        axial = w3**(degree-2)/factorial(degree-2)
        vector = axial*(z1+z2)**2/2
        assert sp.diff(vector,z1)-sp.diff(vector,z2) == 0
    print("Two scalar input multiplications are exact; the closed unit survives restriction.")


if __name__ == "__main__":
    two_loop_vertex()
    regular_brane()
