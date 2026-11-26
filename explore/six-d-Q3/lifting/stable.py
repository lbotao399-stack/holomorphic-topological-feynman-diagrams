"""Check stable primitive multiplicities by exact ranks in the graded cyclic complex."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import product
from math import comb


def multiplicity(length, weight):
    support = sum(k>0 for k in weight)
    if not support:
        return length%2
    column = length-1-sum(max(k-1,0) for k in weight[-2:])
    return comb(support-1,column) if 0<=column<support else 0


def local_dimension(length,weight,traceless=True):
    """Coefficient of the graded symmetric algebra on stable primitive classes (N >= length)."""
    zero=(0,)*(len(weight)+1)
    target=(length,*weight)
    counts={zero:1}
    for p in range(2 if traceless else 1,length+1):
        for w in product(*(range(k+1) for k in weight)):
            m=multiplicity(p,w)
            if not m:
                continue
            degree=(p,*w)
            odd=(p+sum(w[-2:]))%2
            updated=counts.copy()
            for old,value in counts.items():
                limit=min((t-j)//i for i,j,t in zip(degree,old,target) if i)
                if odd:
                    limit=min(limit,m)
                for k in range(1,limit+1):
                    new=tuple(j+k*i for j,i in zip(old,degree))
                    ways=comb(m,k) if odd else comb(m+k-1,k)
                    updated[new]=updated.get(new,0)+value*ways
            counts=updated
    return counts.get(target,0)


@cache
def factors(k, length, odd):
    if length == 1:
        return ((k,),) if not odd or k<=1 else ()
    return tuple((a,)+rest for a in range(min(k,1 if odd else k)+1)
                 for rest in factors(k-a,length-1,odd))


def odd_degree(monomial):
    return sum(monomial[-2:])%2


def rotate(word):
    p = odd_degree(word[-1])
    sign = (-1)**(len(word)-1+p*sum(map(odd_degree,word[:-1])))
    return (word[-1],)+word[:-1],sign


@cache
def representative(word):
    rotated = word
    sign = 1
    orbit = {}
    while rotated not in orbit:
        orbit[rotated] = sign
        rotated, s = rotate(rotated)
        sign *= s
    canonical = min(orbit)
    return canonical, orbit[canonical] if sign==1 else 0


def cyclic_words(length, weight):
    possibilities = [factors(k,length,i>=3) for i,k in enumerate(weight)]
    for parts in product(*possibilities):
        word = tuple(zip(*parts))
        canonical,sign = representative(word)
        if word==canonical and sign:
            yield word


def multiply(a,b):
    if any(x and y for x,y in zip(a[-2:],b[-2:])):
        return (),0
    return tuple(x+y for x,y in zip(a,b)),(-1)**(a[-1]*b[-2])


def boundary(word):
    output = Counter()
    if len(word)==1:
        return output
    for i in range(len(word)):
        if i+1==len(word):
            ordered,sign = rotate(word)
            j=0
        else:
            ordered,sign,j = word,(-1)**i,i
        factor,local = multiply(ordered[j],ordered[j+1])
        if local:
            target,s = representative(ordered[:j]+(factor,)+ordered[j+2:])
            output[target] += sign*local*s
    return {k:v for k,v in output.items() if v}


def chain_ranks(length,weight):
    pivots={}
    dimension=0
    for word in cyclic_words(length,weight):
        dimension+=1
        col=boundary(word)
        square=Counter()
        for target,v in col.items():
            for other,m in boundary(target).items():
                square[other]+=v*m
        assert not any(square.values())
        while col:
            pivot=min(col)
            if pivot not in pivots:
                pivots[pivot]={k:Fraction(v,col[pivot]) for k,v in col.items()}
                break
            scalar=col[pivot]
            for k,v in pivots[pivot].items():
                value=col.get(k,0)-scalar*v
                if value:
                    col[k]=value
                else:
                    col.pop(k,None)
    return dimension,len(pivots)


if __name__=='__main__':
    weights = sorted({(0,0,0,a,b) for a in range(4) for b in range(4)}
                     | {(1,0,0,a,b) for a in range(3) for b in range(3)}
                     | {(1,1,0,0,0),(1,1,1,0,0),(2,1,0,1,0),(1,0,0,1,2)})
    for w in weights:
        ranks=[chain_ranks(d,w) for d in range(1,6)]
        for d in range(1,5):
            dimension,rank=ranks[d-1]
            assert dimension-rank-ranks[d][1] == multiplicity(d,w), (d,w)
    for w in weights:
        for length in (2,3):
            assert local_dimension(length,w)==multiplicity(length,w)
    # At zero weight the absolute traceless primitive ghosts have odd lengths 3,5,7,...
    assert local_dimension(8,(0,0,0,0,0))==1
    print(len(weights)*4, 'exact cyclic homology dimensions and all corresponding boundary squares verified')
