"""Noncommutative descent identities and stable integrated primitive dimensions."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import factorial
from sympy import SparseMatrix


def cyclic(word):
    # C is odd; F=dC is even. Canonicalize the ordinary graded trace word.
    sign=1
    cycle={}
    while word not in cycle:
        cycle[word]=sign
        if word[0]=='C' and word[1:].count('C')%2:
            sign=-sign
        word=word[1:]+word[:1]
    key=min(cycle)
    return key,cycle[key] if cycle[word]==sign else 0


def descended(length):
    # Coefficients of L int_0^1 tr(C (t F + t(1-t) C^2)^(L-1)) dt.
    from itertools import product
    counted=defaultdict(Fraction)
    for factors in product(('F','CC'),repeat=length-1):
        word,sign=cyclic('C'+''.join(factors))
        pairs=factors.count('CC')
        counted[word]+=sign*Fraction(factorial(length)*factorial(pairs),factorial(length+pairs))
    return {key:v for key,v in counted.items() if v}


def total_derivative(words):
    output=defaultdict(Fraction)
    for word,value in words.items():
        for j,ch in enumerate(word):
            sign=(-1)**word[:j].count('C')
            # (d+Qhat) C = F+CC and (d+Qhat) F = [C,F].
            terms=(('F',1),('CC',1)) if ch=='C' else (('CF',1),('FC',-1))
            for replacement,s in terms:
                target,t=cyclic(word[:j]+replacement+word[j+1:])
                output[target]+=t*s*sign*value
    return {key:v for key,v in output.items() if v}


def forms(k):
    for number in range(min(k,3)+1):
        for spatial in combinations(range(3),number):
            for first in range(k-number+1):
                for theta in range(4):
                    yield spatial,first,k-number-first,theta


def derivative_rank(k,odd_weight):
    # At zero spatial Taylor order only the two odd-coordinate derivatives remain.
    # dz is odd; dtheta is even and can have arbitrary powers.
    basis=lambda m:[t for t in forms(m) if (t[1]+t[2]+t[3].bit_count())%2==odd_weight]
    source,target=basis(k),{x:i for i,x in enumerate(basis(k+1))}
    entries={}
    for j,(sp,a,b,theta) in enumerate(source):
        for bit in range(2):
            if theta & (1<<bit):
                sign=(-1)**((theta & ((1<<bit)-1)).bit_count()+len(sp))
                row=target[sp,a+(bit==0),b+(bit==1),theta^(1<<bit)]
                entries[row,j]=sign
    return len(source),SparseMatrix(len(target),len(source),entries).rank()


def integrated_primitive(p):
    result=[0,0]
    for length in range((p+2)//2,p+1):
        k=2*length-1-p
        if not 0<=k<length:
            continue
        for weight in (0,1):
            dim,out=derivative_rank(k,weight)
            prev=derivative_rank(k-1,weight)[1] if k else 0
            # Incoming at the top of the truncated complex is absent.
            result[(p+weight)%2]+=dim-prev-(out if k<length-1 else 0)
    return tuple(result)


def direct_translation_ranks(length, weight):
    """Cyclic chains dual to densities at total momentum zero lie in ker(total d/dz)."""
    from stable import boundary, cyclic_words, representative
    from sympy.polys.matrices import DomainMatrix
    words = list(cyclic_words(length,weight))
    joint, derivatives = [], []
    for word in words:
        column=defaultdict(int)
        for slot,mon in enumerate(word):
            for j in range(3):
                if mon[j]:
                    reduced=list(mon)
                    reduced[j]-=1
                    out,sign=representative(word[:slot]+(tuple(reduced),)+word[slot+1:])
                    column[j,out]+=mon[j]*sign
        derivatives.append(dict(column))
        column.update({(3,out):v for out,v in boundary(word).items()})
        joint.append(column)

    def rank(columns):
        rows={}
        matrix={}
        for j,column in enumerate(columns):
            for key,value in column.items():
                if value:
                    matrix[rows.setdefault(key,len(rows)),j]=value
        return DomainMatrix.from_Matrix(SparseMatrix(len(rows),len(words),matrix)).rank()

    r=rank(derivatives)
    return len(words)-r,rank(joint)-r


def weighted_model(p,weight):
    if any(k>1 for k in weight[:3]):
        return 0
    def count(k):
        return sum(tuple(int(i in sp) for i in range(3))==weight[:3]
                   and (a+(theta&1),b+(theta>>1))==weight[-2:]
                   for sp,a,b,theta in forms(k))
    def rank(k):
        if not any(weight[-2:]):
            return 0
        r=0
        for j in range(k+1):
            r=count(j)-r
        return r

    dimension=0
    for length in range((p+2)//2,p+1):
        k=2*length-1-p
        if 0<=k<length:
            dimension+=count(k)-rank(k-1)-(rank(k) if k<length-1 else 0)
    return dimension


def check_total_momentum():
    from functools import cache
    ranks=cache(direct_translation_ranks)
    spatial=[(0,0,0),(1,0,0),(1,1,0),(2,0,0),(1,1,1)]
    odd=[(0,0),(1,0),(1,1),(2,0)]
    for ev in spatial:
        for od in odd:
            w=ev+od
            for p in (2,3):
                dim,rank=ranks(p,w)
                assert dim-rank-ranks(p+1,w)[1]==weighted_model(p,w),(p,w)
    return len(spatial)*len(odd)*2


if __name__=='__main__':
    for length in range(1,9):
        assert total_derivative(descended(length))=={'F'*length:1}
    assert integrated_primitive(2)==(9,8)
    print(check_total_momentum(),'independent zero-total-momentum cyclic ranks verified')
    cubic=integrated_primitive(3)
    assert sum(cubic)==35
    print('descent identities through degree eight; integrated quadratic and cubic parity dimensions',
          integrated_primitive(2),cubic)
