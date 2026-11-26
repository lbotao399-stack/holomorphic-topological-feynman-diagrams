"""Exhaustive sl2 finite-weight comparison in three even and two odd coordinates."""

from itertools import combinations_with_replacement
from math import factorial, prod
from collections import Counter

from relative import Sector


def orbit(entries):
    return factorial(len(entries)) // prod(factorial(n) for n in Counter(entries).values())


def weights(limit):
    for hol in combinations_with_replacement(range(limit+1),3):
        if not 2<=sum(hol)<=limit:
            continue
        for odd in combinations_with_replacement(range(limit-sum(hol)+1),2):
            if sum(hol)+sum(odd)<=limit:
                yield tuple(reversed(hol))+tuple(reversed(odd)), orbit(hol)*orbit(odd)


def scan(limit):
    from index_scan import index
    from diagonal_products import count
    data=[]
    for w,mult in weights(limit):
        sector=Sector(w)
        dims=dict(sector.analyze())[None]
        assert all(a==b for a,b in dims.values()), (w,dims)
        if sum(w)<=6:
            assert {k:a for k,(a,b) in dims.items()} == sector.adjacent_ranks(), w
        assert sum((-1)**k*a for k,(a,b) in dims.items()) == index(w)[w]
        assert sum(a for a,b in dims.values()) == count(w).get(w,0), w
        data.append((w,mult,dims))
    return data


if __name__=="__main__":
    import sys
    limit=int(sys.argv[1]) if len(sys.argv)>1 else 8
    rows=scan(limit)
    total=sum(mult*sum(a for a,b in dims.values()) for _,mult,dims in rows)
    print(f"degree through {limit}: {len(rows)} permutation orbits, {sum(m for _,m,_ in rows)} weights, {total} relative classes, 0 in restriction kernel")
