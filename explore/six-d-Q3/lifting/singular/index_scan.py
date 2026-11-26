"""Alternating singlet index by Cartan weights, without building CE matrices."""

from collections import defaultdict
from itertools import product


def index(bound):
    zero = (0,)*5
    mon = {(zero,0): 1}
    order = sorted((w for w in product(*(range(n+1) for n in bound[:3]),range(2),range(2))
                    if w!=zero and all(a<=b for a,b in zip(w,bound))), key=lambda w:(sum(w),w))
    for w in order:
        odd = sum(w[3:]) != 1
        for color in (1,0,-1):
            new = defaultdict(int)
            for (deg,charge), val in mon.items():
                new[deg,charge] += val
                maxdeg = min((b-a)//c for a,b,c in zip(deg,bound,w) if c)
                if odd:
                    maxdeg = min(1,maxdeg)
                for k in range(1,maxdeg+1):
                    dest = tuple(a+k*c for a,c in zip(deg,w))
                    new[dest,charge+k*color] += val * (-1)**k
            mon = {key: val for key,val in new.items() if val}
    return {deg:mon.get((deg,0),0)-mon.get((deg,1),0)
            for deg in product(*(range(n+1) for n in bound))}


def degree_totals(limit):
    """Sum over labeled source multidegrees using body monomial multiplicities."""
    from math import comb
    mon = {(0,0,0):1}
    for h in range(limit+1):
        for odd_weight in (0,1,2):
            if not (0<h+odd_weight<=limit):
                continue
            mult = comb(h+2,2)*comb(2,odd_weight)
            exterior = odd_weight != 1
            for color in (1,0,-1):
                new = defaultdict(int)
                for (hh,tt,c),n in mon.items():
                    cap = (limit-hh-tt)//(h+odd_weight)
                    if exterior:
                        cap = min(mult,cap)
                    for k in range(cap+1):
                        copies = comb(mult,k) if exterior else comb(mult+k-1,k)
                        new[hh+k*h,tt+k*odd_weight,c+k*color] += n * (-1)**k * copies
                mon = {key:val for key,val in new.items() if val}
    return {(h,t):mon.get((h,t,0),0)-mon.get((h,t,1),0)
            for h in range(limit+1) for t in range(limit-h+1)}


if __name__ == "__main__":
    import sys
    bound = tuple(map(int,sys.argv[1:])) if len(sys.argv)>1 else (2,2,1,4,4)
    values = index(bound)
    negative = sorted(((w, val) for w,val in values.items() if val<0), key=lambda item:sum(item[0]))
    print(f"{len(values)} weights, {len(negative)} negative, first: {negative[:5]}")
