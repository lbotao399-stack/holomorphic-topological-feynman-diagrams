"""Exterior multiplication of the classical rank-five scalar-source point classes."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import factorial


def exterior(left,right):
    result = Counter()
    for a,x in left.items():
        for b,y in right.items():
            if a & b:
                continue
            sign = 1
            remaining = a
            while remaining:
                lowest = remaining & -remaining
                if (b & (lowest-1)).bit_count()%2:
                    sign = -sign
                remaining -= lowest
            result[a|b] += sign*x*y
    return {mask:x for mask,x in result.items() if x}


def mask(*bits):
    return sum(1<<i for i in bits)


def divided_power_traces(field,weight):
    # Bits 0..4 are the polyvector variables. The remaining bits carry odd
    # field coefficients. The five theta variables stand to their left.
    power = field
    out = []
    for number in range(2,6):
        power = exterior(power,field)
        terms = exterior(weight,power)
        top = {k>>5:x for k,x in terms.items() if k&31==31}
        assert all(x%factorial(number)==0 for x in top.values())
        out.append({k:x//factorial(number) for k,x in top.items()})
    return out


def certificate():
    # A rank-four bivector has this normal form over the reals or complexes.
    # Its orbit is dense. The 5 vector components are independent odd generators.
    field = {mask(i,5+i):1 for i in range(5)}
    field.update({mask(0,1):1,mask(2,3):1})
    # Weight scalar, bivector and seven copies of the positive-rank field
    # already account for at least nine of the ten slots in a nonzero product.
    # Test separately every admissible trivector and every weight vector.
    for trivector,weight_vector in [((),())]+[(s,()) for s in combinations(range(5),3)]+[
            ((),(i,)) for i in range(5)]:
        source = field.copy()
        weight = {0:1}
        if trivector:
            source[mask(*trivector,10)] = 1
        if weight_vector:
            weight[mask(*weight_vector,10)] = 1
        first,second,third,fourth = divided_power_traces(source,weight)
        assert exterior(first,fourth)==exterior(second,third)
    # A witness for all the products O_0 O_s. Its positive-rank primary
    # field has only one vector component, so terms with three or more
    # copies of the vector disappear.
    witness = {mask(0,1):1,mask(2,3):1,mask(4,5):1,mask(2,3,4,6):1}
    first,second,third,fourth = divided_power_traces(witness,{0:1})
    assert len(first)==len(second)==1 and not third and not fourth
    orientation = exterior(first,second)[3]
    assert abs(orientation)==1
    # On this field the trace generating series is exp(phi*z)*(first+z*second).
    def trace_at(index, phi):
        trace = Counter({k:Fraction(phi**index,factorial(index))*v for k,v in first.items()})
        if index:
            for k,v in second.items():
                trace[k] += Fraction(phi**(index-1),factorial(index-1))*v
        return trace

    factors = {}
    for r,s in combinations(range(5),2):
        m = r+s
        factors[r,s] = Fraction((s-r)*factorial(m-1),factorial(r)*factorial(s))
        for phi in range(8):
            value = exterior(trace_at(r,phi),trace_at(s,phi)).get(3,0)
            expected = orientation*Fraction(phi**(m-1),factorial(m-1))
            assert value==factors[r,s]*expected
    assert {r+s for r,s in factors}==set(range(1,8))
    assert [factors[pair] for pair in [(1,4),(2,4),(3,4)]]==[3,5,5]
    print('Sixteen exterior tests and ten protected pairs passed; the seven constant symbols are independent.')


if __name__=='__main__':
    certificate()
