"""Exact SL(5) multiplicities of positive divergence-free jets with a current or scalar output."""
from collections import Counter
from itertools import combinations
from math import comb
import numpy as np


def nodes(n, order):
    """Rotation representatives of unordered distinct roots, including 1."""
    for others in combinations(range(1, order), n - 1):
        indices = (0,) + others
        rotations = {tuple(sorted((i - j) % order for i in indices)) for j in indices}
        if indices == min(rotations):
            # If the stabilizer is nontrivial this orbit is shorter.
            yield indices, order * len(rotations) // n


def exterior(characters, highest_degree, max_arity, inv_arity, reduce):
    # Newton identities in exterior arity, graded by the positive Euler degree.
    coefficients = [[0] * (highest_degree + 1) for _ in range(max_arity + 1)]
    coefficients[0][0] = 1
    for q in range(1, max_arity + 1):
        for degree in range(q, highest_degree + 1):
            total = 0
            for r in range(1, q + 1):
                for h in range(1, (degree - (q-r)) // r + 1):
                    prev = coefficients[q-r][degree-r*h]
                    total = reduce(total + (-1)**(r-1) * prev * characters[r][h])
            coefficients[q][degree] = inv_arity(total, q)
    return coefficients


def bounds(n, degree, max_arity):
    # Dimensions at the identity bound the multiplicities integrally.
    dimensions = [n*comb(n+h, n-1)-comb(n+h-1, n-1) for h in range(degree+1)]
    def divide(x,q):
        assert x % q == 0
        return x//q
    return exterior([dimensions]*(max_arity+1),degree,max_arity,divide,lambda x:x)


def torus_data(n, order, qmax, scalar_degree, determinant_power, prime, generator):
    assert (prime - 1) % order == 0
    primitive = pow(generator, (prime-1)//order, prime)
    assert len({pow(primitive,i,prime) for i in range(order)}) == order
    base = np.array([pow(primitive,i,prime) for i in range(order)], dtype=np.int64)
    points, sizes = zip(*nodes(n,order))
    positions = np.array(points)
    # The full Weyl average is a sum over unordered distinct n-tuples.
    factors = np.array(sizes,dtype=np.int64) * pow(order**n,-1,prime) % prime
    factors = factors * base[(-determinant_power*positions.sum(axis=1)) % order] % prime
    for i in range(n):
        for j in range(n):
            if i != j:
                factors = factors * (1-base[(positions[:,i]-positions[:,j])%order]) % prime
    characters = [None]
    for r in range(1,qmax+1):
        roots = base[(r*positions)%order]
        duals = base[(-r*positions)%order]
        complete = [1]+[0]*(scalar_degree+1)
        for z in duals.T:
            for j in range(1, len(complete)):
                complete[j] = (complete[j] + z*complete[j-1]) % prime
        vector = roots.sum(axis=1) % prime
        characters.append([(vector*complete[h+1]-complete[h])%prime
                           for h in range(scalar_degree+1)])
    covector = base[(-positions)%order].sum(axis=1) % prime
    return factors, covector, characters


def multiplicities(n, order, qmax, scalar_degree, determinant_power, prime, generator):
    factors, covector, characters = torus_data(n, order, qmax, scalar_degree, determinant_power, prime, generator)
    coeff = exterior(characters,scalar_degree,qmax,
                     lambda x,q:x*pow(q,-1,prime)%prime,lambda x:x%prime)
    scalar = [int(((factors*coeff[q][scalar_degree])%prime).sum() % prime)
              for q in range(qmax+1)]
    current = [int((((factors*covector%prime)*coeff[q][scalar_degree-1])%prime).sum()%prime)
               for q in range(qmax+1)]
    return current, scalar


def dimensions(n, order, qmax, determinant_power):
    degree = -n*determinant_power
    assert degree>=1
    # Exponents of a monomial after the determinant and Weyl shifts lie
    # between -(degree+qmax)+(1-n)-determinant_power and qmax+(n-1)-determinant_power.
    assert order > max(degree+qmax+n-1+determinant_power,qmax+n-1-determinant_power)
    dim = bounds(n,degree,qmax)
    biggest = max(max(row[degree],n*row[degree-1]) for row in dim)
    mod_product=1
    answer=[[0]*(qmax+1) for _ in range(2)]
    for prime,gen in [(998244353,3),(469762049,3),(167772161,3),(754974721,11),
                      (2013265921,31),(1811939329,13)]:
        new=multiplicities(n,order,qmax,degree,determinant_power,prime,gen)
        for prev,cur in zip(answer,new):
            for q,v in enumerate(cur):
                prev[q]+=mod_product*((v-prev[q])*pow(mod_product,-1,prime)%prime)
        mod_product*=prime
        if mod_product>biggest:
            break
    assert mod_product>biggest
    for q in range(qmax+1):
        assert answer[0][q]<=n*dim[q][degree-1] and answer[1][q]<=dim[q][degree]
    return answer


def profiles(q, total, lower=1):
    if q == 0:
        if total == 0:
            yield ()
        return
    for h in range(lower,total//q+1):
        for rest in profiles(q-1,total-h,h):
            yield (h,)+rest


def profile_dimensions(total_dimensions):
    # A single residue is exact since all profile dimensions are nonnegative
    # and their sums are already certified below the prime.
    prime = 998244353
    factors,covector,characters = torus_data(5,32,9,15,-3,prime,3)
    exteriors={}
    for h in range(1,16):
        e=[1]
        for k in range(1,min(9,15//h)+1):
            acc=0
            for r in range(1,k+1):
                acc=(acc+(-1)**(r-1)*characters[r][h]*e[k-r])%prime
            e.append(acc*pow(k,-1,prime)%prime)
        exteriors[h]=e
    for q,total,measure,certified in [(9,14,factors*covector%prime,total_dimensions[0][9]),
                                      (8,15,factors,total_dimensions[1][8])]:
        assert certified < prime
        result={}
        for profile in profiles(q,total):
            character=measure.copy()
            for h,k in Counter(profile).items():
                character=character*exteriors[h][k]%prime
            result[tuple(x+1 for x in profile)]=int(character.sum()%prime)
        assert sum(result.values())==certified
        yield result


if __name__=='__main__':
    current,scalar=dimensions(5,32,15,-3)
    for q,(j,k) in enumerate(zip(current,scalar)):
        print(f'{q:2} {j:8} {k:8}')
    for row in (current,scalar):
        print('Euler characteristic',sum((-1)**q*v for q,v in enumerate(row)))
    for result in profile_dimensions((current,scalar)):
        for profile,size in result.items():
            print(profile,size)
