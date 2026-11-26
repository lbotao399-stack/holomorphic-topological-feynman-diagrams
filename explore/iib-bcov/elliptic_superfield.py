"""Full exterior-algebra tests of the all-arity genus-one bulk Ward identity."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import factorial
from random import Random

from elliptic_hessian import add, determinant, momenta_by_subset, scale
from genus_one_graph_ward import plus


@lru_cache(maxsize=16384)
def exterior_sign(left, right):
    parity = 0
    while left:
        bit = left & -left
        parity += (right & (bit-1)).bit_count()
        left ^= bit
    return (-1)**parity


def multiply(left,right):
    result = {}
    for (a,m),x in left.items():
        for (b,n),y in right.items():
            if a&b or m&n:
                continue
            key = (a|b,m|n)
            result[key] = result.get(key,0)+x*y*exterior_sign(m,n)
    return {k:v for k,v in result.items() if v}


def ward(momenta, levels, ranks, rng):
    n = len(momenta)
    unit = {(0,0):1}
    charges = momenta_by_subset(momenta)
    def dx(field,j):
        return {key:v*charges[key[0]][j] for key,v in field.items() if charges[key[0]][j]}
    def odd(field,j):
        bit = 1<<(n+j)
        return {(w,m^bit):v*exterior_sign(bit,m^bit) for (w,m),v in field.items() if m&bit}
    def delta(field):
        return add(*(dx(odd(field,j),j) for j in range(5)))
    def bracket(even,field):
        return add(delta(multiply(even,field)),scale(multiply(delta(even),field),-1),
                   scale(multiply(even,delta(field)),-1))
    def determinant_rows(rows):
        return determinant(rows,multiply,unit)
    def hessian(field):
        return [[dx(dx(field,j),i) for j in range(5)] for i in range(5)]
    def directional_det(matrix,direction):
        return add(*(determinant_rows([direction[j] if j==i else row for j,row in enumerate(matrix)])
                     for i in range(5)))
    def integrated(field):
        # Conserved total momentum and top theta degree select the integrated n-wave trace.
        top = ((1<<5)-1)<<n
        return {m^top:v for (w,m),v in field.items() if w==(1<<n)-1 and m&top==top}

    fields = [{} for j in range(5)]
    for i,(k,r) in enumerate(zip(levels,ranks)):
        assert 0<=r<=4-k
        coefficient = (1<<i) if r%2 else 0
        for directions in combinations(range(5),r):
            mask = coefficient | sum(1<<(n+j) for j in directions)
            fields[k][1<<i,mask] = rng.choice((-2,-1,1,2))
    def derivative(u,order):
        out = {}
        for j in reversed(range(order,5)):
            out = add(multiply(out,u),scale(fields[j],Fraction(1,factorial(j-order))))
        return out
    u = {}
    for _ in range(n):
        u = derivative(u,0)
    prime = derivative(u,1)
    h, log_h, power = unit, {}, unit
    for k in range(1,n+1):
        power = multiply(power,prime)
        if not power:
            break
        h, log_h = add(h,power),add(log_h,scale(power,Fraction(1,k)))
    powers = [unit]
    for _ in range(5):
        powers.append(multiply(powers[-1],u))
    f0 = add(scale(powers[2],Fraction(-1,2)),
             *(scale(multiply(field,powers[j+1]),Fraction(1,factorial(j+1)))
               for j,field in enumerate(fields)))
    qa = [scale(delta(f0),-1),*(delta(field) for field in fields[:-1])]
    qu = multiply(h,add(*(scale(multiply(powers[j],qa[j]),Fraction(1,factorial(j)))
                          for j in range(5))))
    qlog = multiply(h,add(multiply(derivative(u,2),qu),
                           *(scale(multiply(powers[j-1],qa[j]),Fraction(1,factorial(j-1)))
                             for j in range(1,5))))
    assert qu==scale(bracket(u,u),Fraction(-1,2))
    assert qlog==add(delta(u),scale(bracket(u,log_h),-1))
    matrix = hessian(u)
    det = determinant_rows(matrix)
    dx_u = [dx(u,i) for i in range(5)]
    pair = [[bracket(dx_u[i],dx_u[j]) for j in range(5)] for i in range(5)]
    assert directional_det(matrix,pair)==scale(multiply(det,delta(u)),2)
    density = multiply(log_h,det)
    assert integrated(bracket(u,density))==integrated(scale(multiply(density,delta(u)),-2))
    variation = add(multiply(qlog,det),multiply(log_h,directional_det(matrix,hessian(qu))))
    value = integrated(variation)
    assert value==integrated(multiply(delta(u),det))
    return value


def check():
    rng = Random(151127)
    cases = [((0,)*6,(1,)*6), ((0,)*6,(0,0,1,1,2,2)),
             ((0,0,0,0,0,0,1),(0,0,0,1,1,2,2)),
             ((0,0,0,0,0,0,0,2),(0,0,0,1,0,2,1,2)),
             ((0,0,0,0,0,0,1,1),(0,0,1,1,1,0,1,2)),
             ((0,0,0,0,0,0,0,0,3),(0,0,0,0,0,0,4,1,1)),
             ((0,0,0,0,0,0,0,0,0,4),(0,0,0,0,0,0,0,4,2,0))]
    for levels,ranks in cases:
        waves = tuple(tuple(rng.randrange(-2,3) for j in range(5)) for i in levels[1:])
        waves += (tuple(-x for x in plus(waves)),)
        value = ward(waves,levels,ranks,rng)
        assert value
        print('levels',levels,'polyvector ranks',ranks,'signed trace',list(value.values()),flush=True)


if __name__=='__main__':
    check()
