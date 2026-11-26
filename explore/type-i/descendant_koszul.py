"""Polynomial contraction of the one-closed-field differential q wedge on holomorphic jets."""
from collections import Counter
from itertools import combinations


def sign(slots,i):
    return (-1)**sum(j<i for j in slots)


def differential(source):
    out=Counter()
    for (exponents,slots),c in source.items():
        for i in range(5):
            if i not in slots:
                powers=list(exponents);powers[i]+=1
                out[tuple(powers),tuple(sorted((*slots,i)))]+=sign(slots,i)*c
    return out


def homotopy(source):
    out=Counter()
    for (exponents,slots),c in source.items():
        for i in slots:
            if exponents[i]:
                powers=list(exponents);powers[i]-=1
                out[tuple(powers),tuple(j for j in slots if j!=i)]+=sign(slots,i)*exponents[i]*c
    return out


if __name__=='__main__':
    count=0
    for k in range(6):
        for indices in combinations(range(5),k):
            for powers in ((0,0,0,0,0),(1,0,0,0,0),(0,1,2,0,3),(5,3,1,4,2)):
                key=powers,indices
                value={key:1}
                result=differential(homotopy(value))
                result.update(homotopy(differential(value)))
                assert {a:b for a,b in result.items() if b}==({key:sum(powers)+5-k} if sum(powers)+5-k else {})
                count+=1
    print('polynomial Koszul homotopy identities',count)
