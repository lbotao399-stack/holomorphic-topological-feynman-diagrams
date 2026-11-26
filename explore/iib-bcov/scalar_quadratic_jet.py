"""The full mixed-wave projection of Q on a fixed scalar point function and quadratic background."""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod

from flint import fmpq_mat
from scalar_fan_waves import reduced


def quadratic_coefficient(n, hessian):
    """Coefficient for the scalar exponential, with the universal heat and string factors stripped."""
    if n<0 or len(hessian)!=3 or any(len(row)!=3 for row in hessian):
        raise ValueError("A nonnegative genus and a three-dimensional Hessian are required.")
    matrix = fmpq_mat([[str(x) for x in row] for row in hessian])
    if matrix!=matrix.transpose():
        raise ValueError("The Hessian must be symmetric.")
    determinant = F(str(matrix.det()))
    return determinant**n/(2**n*factorial(2*n+1)*factorial(n+1)*factorial(n))


def polarized_coefficient(n, directions, strengths):
    retained = (0,0,0,1,0)
    vector = (0,0,0,0,1)
    total = F(0)
    # A direction occurring more than n times cannot occupy distinct triples.
    for counts in product(range(n+1),repeat=len(directions)):
        if sum(counts)!=3*n:
            continue
        waves = [(*p,0,0) for p,k in zip(directions,counts) for _ in range(k)]
        weight = prod(a**k/factorial(k) for a,k in zip(strengths,counts))
        total += weight*reduced(waves,retained,vector,retained)
    return total


def check():
    directions = ((1,0,0),(0,1,0),(0,0,1),(1,2,-1))
    strengths = (F(1),F(-2),F(3),F(1,2))
    hessian = [[sum(a*p[i]*p[j] for a,p in zip(strengths,directions))
                for j in range(3)] for i in range(3)]
    assert fmpq_mat([[str(x) for x in row] for row in hessian]).det()==-4
    for n in range(1,5):
        assert quadratic_coefficient(n,hessian)==polarized_coefficient(n,directions,strengths)
    identity = [[int(i==j) for j in range(3)] for i in range(3)]
    assert [quadratic_coefficient(n,identity) for n in range(4)]==[
        F(1),F(1,24),F(1,5760),F(1,5806080)]
    from longitudinal_theta import primary_coefficient
    plane = [(*p,0,0) for p in (*directions,(1,-1,1))]
    scalar,vector = (0,0,0,1,0),(0,0,0,0,1)
    for count,expected in ((4,F(11,17280)),(5,F(31,2160))):
        for scale in (1,2,-1):
            waves = [tuple(scale*x for x in p) for p in plane[:count]]
            assert primary_coefficient([*waves,scalar],vector,scalar)==expected*scale**12
    print('Quadratic-jet projection verified through four loops and both lower-order two-loop symbols checked.')


if __name__=='__main__':
    check()
