"""Independent cubature for the four polarized marked BCOV theta periods."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import permutations
from math import factorial,prod
from flint import fmpz_mat
from marked_theta import MOMENTA,CURRENT,partitions

UNIT=[[1,0,0],[0,1,0],[0,0,1],[-1,-1,-1]]


def compositions(total,size):
    if size==1:
        yield (total,)
    else:
        for x in range(total+1):
            for tail in compositions(total-x,size-1):
                yield (x,*tail)


def rule():
    # Uniform four-simplex. Verify all degree-seven barycentric moments.
    nodes=[]
    for j in range(4):
        d=11-2*j
        weight=F((-1)**j*d**7,64*factorial(j)*factorial(11-j))
        for b in compositions(3-j,5):
            nodes.append((tuple(2*x+1 for x in b),d,weight))
    for powers in compositions(7,5):
        assert sum(w*F(prod(i**k for i,k in zip(b,powers)),d**7)
                   for b,d,w in nodes)==F(prod(factorial(k) for k in powers),factorial(11))
    return nodes


@cache
def cube():
    nodes=Counter()
    for gaps,den,w in rule():
        times=tuple(sum(gaps[:i]) for i in range(1,5))
        for assignment in permutations(times):
            nodes[assignment,den]+=w
    assert sum(nodes.values())==1
    return nodes


@cache
def integrals(paths):
    # Repeated p1 and p2 span both cycle components in their directions.
    # Removing these rows and columns leaves the actual six-by-six minor.
    zero=[0,0,0]
    columns=[(p if h!=2 else zero)+(p if h!=0 else zero) for p,h in zip(UNIT,paths)]
    moments=Counter()
    for (times,d),weight in cube().items():
        v=[times[i]-times[3] for i in range(3)]
        square=int(fmpz_mat([v+zero,zero+v,*columns]).det())**2
        if not square:
            continue
        # Factor 2*d^2 clears both averaged Green terms. Last four labels
        # remain ordered. The four repeated-axis times are independent uniforms.
        for at,sign in ((2,1),(3,-1)):
            x=times[at]
            for j,y in enumerate(times):
                if CURRENT[j+4]:
                    xy=2*x*y
                    bridge=2*d*min(x,y)-xy if paths[at]==paths[j] else 0
                    moments[(d,weight,0)]+=sign*square*xy*CURRENT[j+4]
                    moments[(d,weight,1+paths[at])]+=sign*square*bridge*CURRENT[j+4]
            moments[d,weight,0]+=sign*square*4*x*d
            moments[d,weight,4+paths[at]]+=sign*square*x*(d-x)
    result=[F() for _ in range(7)]
    for (d,w,kind),total in moments.items():
        result[kind]+=total*w/(2*d**6)
    return result


def run(shape):
    total=F()
    for groups in partitions(shape):
        path={i:h for h,grp in enumerate(groups) for i in grp}
        if path[0]==path[1] or path[2]==path[3]:
            continue
        vals=integrals(tuple(path[i] for i in range(4,8)))
        # Integrate the three conductances by their independent simplex moments.
        powers=[6-p for p in shape]
        moment=F(prod(factorial(k) for k in powers),factorial(9))
        total+=moment*vals[0]
        for j,p in enumerate(powers):
            number=sum(path[k]==j for k in range(4))
            total+=moment*F(9,p)*(vals[1+j]+number*vals[4+j])
    return total/(2*prod(factorial(n) for n in Counter(shape).values()))


if __name__=='__main__':
    expected=(F(1,2430),F(17,17010),F(1,1890),F(1,1890))
    actual=tuple(run(shape) for shape in ((1,5,5),(2,4,5),(3,3,5),(3,4,4)))
    assert actual==expected and sum(actual)==F(1,405)
    print('Four independent ordered-region cubatures:',actual,'sum',sum(actual))
