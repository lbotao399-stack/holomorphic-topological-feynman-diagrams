"""Lower bounds on the transverse bulk boundary rank from the cubic-soft subspace."""
from itertools import combinations_with_replacement
from math import prod

import numpy as np
from flint import fmpz

from soft_trace_dimensions import multiply_by_slot, extract


def trace_dimensions(genus, prime=10000019):
    slots = 2*genus+3
    shape = (2*genus,)*5
    size = tuple(weight+5-i for i,weight in enumerate(shape))
    if not fmpz(prime).is_prime() or prime<=slots or prime**2>=2**63:
        raise ValueError('Need a prime above the Newton recursion length with signed 64-bit residue products.')
    if prod(size)*prime*(slots+22)>=2**63:
        raise ValueError('Slot character sums exceed the finite-array arithmetic bound.')
    unit = np.zeros(size,dtype=np.int64)
    unit[(0,)*5] = 1
    answer = []
    quadratics = [tuple((i==k)+(j==k) for k in range(5))
                  for i,j in combinations_with_replacement(range(5),2)]
    for cutoff in (2,3):
        powers = [unit]
        for m in range(1,slots+1):
            total = np.zeros_like(unit)
            for stride in range(1,m+1):
                previous = powers[m-stride]
                term = multiply_by_slot(previous,stride,2)
                if cutoff==3:
                    for quadratic in quadratics:
                        shift = tuple(stride*k for k in quadratic)
                        if any(v>=limit for v,limit in zip(shift,size)):
                            continue
                        high = tuple(slice(v,None) for v in shift)
                        low = tuple(slice(None,-v or None) for v in shift)
                        term[high] -= previous[low]
                total += term
            powers.append((total%prime)*pow(m,-1,prime)%prime)
        answer.append(extract(powers[slots],shape,prime))
    # A component with exactly h quadratic slots is the symmetric product
    # of h quadratic monomials times the cubic-soft part at n-h slots.
    if slots>=8:
        for h in range(1,min(3,slots-6)+1):
            total = np.zeros_like(unit)
            for factors in combinations_with_replacement(quadratics,h):
                shift = tuple(map(sum,zip(*factors)))
                high = tuple(slice(v,None) for v in shift)
                low = tuple(slice(None,-v or None) for v in shift)
                total[high] += powers[slots-h][low]
            answer.append((answer[-1]+extract(total%prime,shape,prime))%prime)
    return tuple(answer)


if __name__=='__main__':
    # Each second multiplicity is at most the corresponding first multiplicity,
    # independently computed in characteristic zero by soft_target.genus_counts.
    for g,expected in ((2,(1,0)),(3,(18,2,6,12,17)),(4,(2398,122,649,1527,2149)),
                       (5,(1005099,31671,198576,525868,827474))):
        result = trace_dimensions(g)
        assert result==expected
        print(g,*result)
