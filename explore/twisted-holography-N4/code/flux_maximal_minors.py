"""Independent genus-four D3 contact from the two global eleven-slot determinant words."""
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial
from flint import fmpz_mat


@cache
def words():
    subsets=[s for s in combinations(range(11),5) if sum(i<7 for i in s)<=2]
    masks=[sum(1<<i for i in s) for s in subsets]
    first,second=[],[]
    for i,j,k,l in combinations_with_replacement(range(len(masks)),4):
        a,b,c,d=masks[i],masks[j],masks[k],masks[l]
        if (a|b|c|d)!=2047:continue
        twice=(a&b)|(a&c)|(a&d)|(b&c)|(b&d)|(c&d)
        thrice=(a&b&c)|(a&b&d)|(a&c&d)|(b&c&d)
        if thrice.bit_count()!=4:continue
        multiplicity=twice.bit_count(),(a&b&c&d).bit_count()
        if multiplicity==(5,0):first.append((i,j,k,l))
        if multiplicity==(4,1):second.append((i,j,k,l))
    return subsets,first,second


def matrices(slopes,tangents,left,right,normal):
    waves=[(s,d,0,0,0) for s,d in zip(slopes,tangents)]
    waves.extend(((0,0,left,1,0),(0,0,right,1,0),(normal[0],0,normal[1],0,normal[2])))
    waves.append(tuple(-sum(col) for col in zip(*waves)))
    return waves


def totals(slopes,left,right,normal):
    subsets,first,second=words()
    waves=matrices(slopes,(2,1,1,1,1,1,1),left,right,normal)
    minors={s:int(fmpz_mat([waves[k] for k in s]).det()) for s in subsets}
    sq=[minors[s]**2 for s in subsets]
    values=[sum(sq[i]*sq[j]*sq[k]*sq[l] for i,j,k,l in pool) for pool in (first,second)]
    # The signed third restriction is reconstructed with full five-row minors.
    def pair(a,b):
        if a>b:return -pair(b,a)
        return minors[a,b,7,8,9]
    third=0
    for center in range(7):
        for a,b in combinations((i for i in range(7) if i!=center),2):
            rest=[i for i in range(7) if i not in (a,b,center)]
            for j in range(1,4):
                k,l=[rest[i] for i in range(1,4) if i!=j]
                third+=pair(a,b)**2*pair(a,center)*pair(b,center)*pair(rest[0],rest[j])**2*pair(k,l)**2
    return *values,third


def check():
    assert tuple(map(len,words()[1:]))==(7560,2520)
    row=[0,0,0]
    # Inclusion-exclusion at 192 integer points extracts the homogeneous degree-eight coefficient.
    for s in product(range(3),*(range(2) for _ in range(6))):
        scale=(-1)**sum(s)*comb(2,s[0])
        for j,c in enumerate(totals(s,1,0,(2,-1,1))):row[j]+=scale*c
    assert all(v%2==0 for v in row)
    first,second,third=(v//2 for v in row)
    assert (first,second,third)==(51840,259200,360)
    assert 22*first-5*second+24*18*third==0
    # Homogeneity extracts both fourth-degree source jets by finite differences at 25 points.
    s=(1,3,-1,2,5,7,-3);reference=totals(s,1,0,(2,-1,1))
    cross=[0,0,0]
    for t,u in product(range(5),repeat=2):
        weight=(-1)**(t+u)*comb(4,t)*comb(4,u)
        for j,c in enumerate(totals(s,t,u,(-3,7,2))):cross[j]+=weight*c
    assert tuple(cross)==tuple(c*70*2**8*factorial(4)**2 for c in reference)
    k=Fraction(70*second,13824)
    j=Fraction(70*first,27)-Fraction(70*second,216)
    assert (k,j)==(Fraction(2625,2),50400)
    return k,j


if __name__=='__main__':
    print('Direct eleven-slot trace coefficients',*check())
