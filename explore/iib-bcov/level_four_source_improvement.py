"""An exact improvement of the rank-five level-four source on the genus-two probe.

The primitive has a trivector T of momentum k, a vector v, and seven scalar waves.
Its scalar-symmetrized Fourier symbol is
  sum_sigma det(k,s1,s2,s3,s4)^3 T(s2,s3,s4) k(v),
where sigma permutes all seven scalar momenta, including the three spectators.
Determinants use the fivefold volume form; this degree-19 polynomial is a
translation- and SL(5)-equivariant finite-jet expression. Put the odd trivector
coefficient before the odd vector coefficient. Its nine primary input coordinates
have total degree 2*9-3-1=14. Multiplication by lambda^4 gives degree six, and
the quantum differential produces degree seven, as for the rank-five E4 source.
"""
from fractions import Fraction
from itertools import permutations
from math import factorial

from flint import fmpz_mat


def basis(i):
    return tuple(int(i==j) for j in range(5))


def det(rows):
    return int(fmpz_mat(rows).det())


def main():
    volume = tuple(basis(j) for j in range(1,5))
    vector = (1,2,0,0,0)
    scalar = [tuple(sign*x for x in p) for p in volume for sign in (1,-1)]
    assert det((*volume,vector))==1
    assert sum(x*y for x,y in zip((-2,1,0,0,0),vector))==0
    rows = []
    for i,q in enumerate(scalar):
        # Only the rank-four wave W with a scalar q can produce the required trivector.
        # iota_(e1+q) W = iota_q W. Every other fusion leaves the trivector slot empty.
        k = tuple(x+y for x,y in zip(basis(0),q))
        rest = scalar[:i]+scalar[i+1:]
        value = count = 0
        for s in permutations(rest,4):
            determinant = det((k,*s))
            trivector = det([p[1:] for p in (q,*s[1:])])
            term = determinant**3*trivector*sum(x*y for x,y in zip(k,vector))
            if term:
                count += factorial(3)
                value += factorial(3)*term
        assert count==288
        rows.append(value)
    assert rows==[-864,288,-288,-288,-288,-288,-288,-288]
    # Q0 A0=-partial_Omega(A0^2/2) on primaries, with T ordered first.
    variation = -sum(rows)
    assert variation==2304
    source = sum(map(Fraction,('91/360','347/540','3227/6480','137/1215','179/38880')))
    assert source==Fraction(19579,12960)
    assert source-(source/variation)*variation==0
    print('Unit exterior trace; reduced genus-two source',source)
    print('Classical variation of the degree-19 primitive:',variation)
    print('Added primitive coefficient in matching reduced units:',-source/variation)


if __name__=='__main__':
    main()
