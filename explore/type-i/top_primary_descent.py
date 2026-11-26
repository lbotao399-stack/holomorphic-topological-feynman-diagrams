"""Primary rank-five disk: all single and double traces on nine conserved gauge waves."""
from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations,combinations_with_replacement
from one_primary_descent import add,coordinates,project
from external_closed import stable_colour,relabel,parity


@cache
def minor_coordinates(subset):
    # The ninth open momentum is minus the sum of the other eight.
    if 8 not in subset:
        return {variables:value for (_,variables),value in coordinates(subset).items()}
    result=Counter()
    for i in range(8):
        if i not in subset[:-1]:
            choices=(*subset[:-1],i)
            sign=-parity(choices)
            for monomial,value in minor_coordinates(tuple(sorted(choices))).items():
                result[monomial]+=sign*value
    return {m:v for m,v in result.items() if v}


@cache
def chart_product(a,b):
    out=Counter()
    for u,x in minor_coordinates(a).items():
        for v,y in minor_coordinates(b).items():
            out[tuple(sorted(u+v))]+=x*y
    return {m:v for m,v in out.items() if v}


def chart(symbol):
    out=Counter()
    for (a,b),x in symbol.items():
        if x:
            for m,v in chart_product(a,b).items():
                out[m]+=x*v
    return {m:v for m,v in out.items() if v}


def wheel(labels,colours):
    by_rank=[Counter(),Counter()]
    for (power,words),weight in colours.items():
        for s,sign in relabel(words,labels):
            a=tuple(sorted(s[:5]))
            factor=weight*parity(s[:5])*sign
            for i in range(5):
                b=(*s[5:],s[i])
                key=tuple(sorted((a,tuple(sorted(b)))))
                by_rank[power][key]+=(5-i)*factor*parity(b)
    return [chart(symbol) for symbol in by_rank]


def differential(left,right):
    out=Counter()
    for j in range(9):
        order=[(j+i)%9 for i in range(9)]
        momenta=[tuple(order[:2])]+[(i,) for i in order[2:]]
        for direction in [-1,1]:
            for start in range(8):
                s=[(start+direction*i)%8 for i in range(8)]
                weight=-parity(order)*parity(s)
                add(out,tuple(momenta[s[i]] for i in left),tuple(momenta[s[i]] for i in right),weight)
    return chart(out)


def closed_sewing(representative):
    out=Counter()
    for (internal,external),value in representative.items():
        left=[(0,),(1,)]+[(i+2,) for i in internal]
        right=[(0,1) if i==-2 else (i+2,) for i in external]
        add(out,tuple(left),tuple(right),value)
    return chart(out)


def run():
    traces=[(tuple(range(9)),),((0,1),tuple(range(2,9))),((0,1,2),tuple(range(3,9))),
            ((0,1,2,3),tuple(range(4,9)))]
    colours=stable_colour(5)
    values=[wheel(labels,colours) for labels in traces]
    print('chart sizes',[[len(v) for v in row] for row in values],flush=True)
    assert not any(values[2]+values[3])
    target=Counter(values[0][0])
    for k,c in values[0][1].items():target[k]+=32*c
    subsets=combinations(range(7),5)
    cols=((key,differential(*key)) for key in combinations_with_replacement(subsets,2))
    rank,[(rep,rem),(generic,rrem)]=project(cols,[target,values[0][1]])
    print('single-trace CE rank',rank,'representative',rep,'remainder size',len(rem),
          'slope remainder',len(rrem),flush=True)
    assert rank==10 and not rem and not rrem
    from primary_disk import solve
    _,_,disk=solve(3,5,(-2,))
    tree=closed_sewing(disk)
    variations=[]
    for free in [(1,0),(0,1)]:
        _,_,disk=solve(3,5,(-2,),free)
        delta=Counter(closed_sewing(disk))
        delta.subtract(tree)
        variations.append((free,delta))
    for multiple in (-48,48):
        remainder=Counter(values[1][0])
        for k,v in tree.items():remainder[k]-=multiple*v
        rank2,[(correction,residual)]=project(variations,[remainder])
        assert rank2==2
        if multiple==-48:
            assert not residual
            free=tuple(correction[i]/multiple for i in [(1,0),(0,1)])
            assert free==(F(81,280),-F(11,30))
            _,_,selected=solve(3,5,(-2,),free)
            assert {key:multiple*v for key,v in closed_sewing(selected).items()}==values[1][0]
            print('unique rank-(3,5) disk coefficients',free,'ordered sewing multiplier',multiple)
        else:
            assert residual
    return values


if __name__=='__main__':
    run()
