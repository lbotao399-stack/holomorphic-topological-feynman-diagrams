"""Signed fan coefficient and SU(3) three-loop check without ribbon reduction."""
from functools import cache
from itertools import permutations,product

import numpy as np


def odd_sign():
    # Coordinate order is theta_1,theta_2,theta_3,beta_1,beta_2,beta_3.
    for order in permutations(range(3)):
        exterior = [3,4,5]
        coefficient = 1
        for coordinate in reversed(order):
            slot = exterior.index(3+coordinate)
            coefficient *= -(-1)**slot
            exterior.pop(slot)
            coefficient *= (-1)**sum(x<coordinate for x in exterior)
            exterior = sorted([coordinate,*exterior])
        assert exterior==[0,1,2] and coefficient==-1


def color_three():
    def unit(i,j):
        result = np.zeros((3,3),dtype=np.int64)
        result[i,j] = 1
        return result

    off_diagonal = [(i,j) for i in range(3) for j in range(3) if i!=j]
    first = np.diag([1,-1,0])
    second = np.diag([0,1,-1])
    basis = [unit(i,j) for i,j in off_diagonal]+[first,second]
    # Three times the trace-dual Chevalley basis. Only traceless generators
    # enter. This is a direct adjoint calculation with integral matrices.
    dual = [3*unit(j,i) for i,j in off_diagonal]+[2*first+second,first+2*second]
    assert [[int(np.trace(x@y)) for y in dual] for x in basis]==(3*np.eye(8,dtype=int)).tolist()
    d = e = unit(0,1)+unit(1,0)
    f = np.diag([1,1,-2])

    @cache
    def symmetrized_input(indices):
        value = 0
        for a,b,c,k in permutations(indices):
            value += int(np.trace(basis[a]@basis[b]@basis[c]@basis[k]@f))
        return value

    total = 0
    for indices in product(range(8),repeat=4):
        source = symmetrized_input(tuple(sorted(indices)))
        if source:
            commutator = d
            for i in indices:
                commutator = dual[i]@commutator-commutator@dual[i]
            total += source*int(np.trace(commutator@e))
    assert total==20412
    coefficient = total//3**4
    assert coefficient==252==(3**3+12*3)*int(np.trace(d@e@f+d@f@e))
    # With D=E=lambda_1/2 divide by four. For the standard Hermitian
    # generators tr(t_A t_B)=delta_AB/2 supplies tau^3, leaving 63/8.
    return coefficient


if __name__=='__main__':
    odd_sign()
    for n in range(2,20):
        source_phase = (-1)**n*1j**(n+1)
        component_phase = -(-1j)**(n-2)*1j**2
        structure_phase = 1j**(-n-1)
        assert source_phase*component_phase*structure_phase==-1j**n
    print('All six Grassmann orders give −theta_1 theta_2 theta_3.')
    print('SU(3) direct nested-commutator integral-basis sum:',color_three())
