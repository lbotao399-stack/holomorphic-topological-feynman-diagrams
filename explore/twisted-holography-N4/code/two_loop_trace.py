"""Two-loop LLC single-trace fan. Exact U(N) ribbons and an independent SU(2) check."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import prod


def contract(words,letter):
    (i,j),(k,l) = [(slot,pos) for slot,word in enumerate(words)
                    for pos,field in enumerate(word) if field==letter]
    word = words[i]
    if i==k:
        return [other for at,other in enumerate(words) if at!=i]+[
                word[j+1:l],word[l+1:]+word[:j]]
    other = words[k]
    return [w for at,w in enumerate(words) if at not in (i,k)]+[
            word[j+1:]+word[:j]+other[l+1:]+other[:l]]


def cyclic(word):
    return min(word[i:]+word[:i] for i in range(len(word))) if word else 'N'


def ribbon_sum(insertion):
    result = Counter()
    for spoke in permutations('ABC'):
        for side in product((1,-1),repeat=3):
            word = 'D'
            for label,orientation in zip(spoke,side):
                word = label+word if orientation==1 else word+label
            words = [insertion,word+'E']
            for letter in 'ABC':
                words = contract(words,letter)
            result[tuple(sorted(map(cyclic,words)))] += prod(side)
    return {words:c for words,c in result.items() if c}


def epsilon(i,j,k):
    if len({i,j,k})<3:
        return 0
    return (-1)**sum(x>y for x,y in ((i,j),(i,k),(j,k)))


def pauli_check():
    # For t_i=sigma_i/2 the trace-pairing coefficient is 1/2.
    # Use D=1,E=2,F=3, for which the structure constant is one.
    total = 0
    for a,b,c in product(range(3),repeat=3):
        trace_times_eight = (int(a==b and c==2)-int(a==c and b==2)+int(a==2 and b==c))
        for x,y,z in permutations((a,b,c)):
            total += trace_times_eight*sum(epsilon(x,0,r)*epsilon(y,r,s)*epsilon(z,s,1)
                                           for r,s in product(range(3),repeat=2))
    assert Fraction(total,8)==-Fraction(5,4)==-Fraction(1,2)**3*(2**2+6)


if __name__=='__main__':
    from single_particle import bdy_B,bulk_coker
    terms = ribbon_sum('ABCF')
    assert terms=={('DEF','N','N'):1,('DFE','N','N'):-1,('DEF',):6,('DFE',):-6}
    assert not ribbon_sum('ABC')
    pauli_check()
    assert bulk_coker((-1,-1,-3,-1,-1))==(1,1,2,0)
    assert bdy_B((0,0,2,0,0),2)==1
    print('Two-loop handle coefficients 1,6. Cubic trace vanishes. SU(2) component −5/4.')
