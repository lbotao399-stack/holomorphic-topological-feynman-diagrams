"""S9 decomposition of products of three five-row minors covering all nine slots."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import factorial, prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from genus_two_probes import det, partitions


@cache
def unrestricted(cycles):
    # Complete symmetric characters on the conserved permutation representation.
    h = [1,-1,0,0,0,0,0,0]
    for length in cycles:
        for j in range(length,len(h)):
            h[j] += h[j-length]
    return det([[h[3-i+j] if 3-i+j>=0 else 0 for j in range(5)] for i in range(5)])


def soft_character(cycles):
    # Each omitted permutation cycle used only the degree-zero slot.
    return sum((-1)**(len(cycles)-k)*sum(unrestricted(c) for c in combinations(cycles,k))
               for k in range(len(cycles)+1))


@cache
def alternant_character(beta,cycles):
    if not cycles:
        return int(beta==tuple(reversed(range(len(beta)))))
    total = 0
    for i,value in enumerate(beta):
        moved = value-cycles[0]
        if moved<0 or moved in beta:
            continue
        sign = (-1)**sum(moved<x<value for x in beta)
        remaining = tuple(sorted((*beta[:i],*beta[i+1:],moved),reverse=True))
        total += sign*alternant_character(remaining,cycles[1:])
    return total


def irreducible(shape,cycles):
    return alternant_character(tuple(v+len(shape)-1-i for i,v in enumerate(shape)),cycles)


def blocks():
    types = list(map(tuple,partitions(9)))
    orders = [prod(j**k*factorial(k) for j,k in Counter(c).items()) for c in types]
    values = [soft_character(c) for c in types]
    with TemporaryDirectory() as temp:
        exe = Path(temp)/'soft-character'
        run(['c++','-O3','-std=c++17',str(Path(__file__).with_name('soft_target_char.cpp')),
             '-o',str(exe)],check=True)
        queries = [(3,0,1,len(c),*(j for length in c for j in (length,1))) for c in types]
        data = '\n'.join(' '.join(map(str,q)) for q in queries)
        replay = run([exe],input=data,text=True,capture_output=True,check=True)
        assert list(map(int,replay.stdout.split()))==values
    result = {}
    for shape in map(lambda c:tuple(reversed(c)),types):
        value = sum(Fraction(t*irreducible(shape,c),order)
                    for t,c,order in zip(values,types,orders))
        assert value.denominator==1 and value>=0
        if value:
            result[shape] = int(value)
    assert all(sum(m*irreducible(shape,c) for shape,m in result.items())==t
               for c,t in zip(types,values))
    assert values[0]==5460 and sum(m*(m+1)//2 for m in result.values())==68
    return result


if __name__=='__main__':
    for shape,multiplicity in blocks().items():
        print(''.join(map(str,shape)),irreducible(shape,(1,)*9),multiplicity)
