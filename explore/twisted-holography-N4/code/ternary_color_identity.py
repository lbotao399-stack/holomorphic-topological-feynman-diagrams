"""Universal metric-Lie-algebra color Jacobi check of the two composed fan terms."""

from collections import defaultdict
from itertools import permutations, product


def rooted_words(vertices):
    """Expand the invariant pairing of external color 5 with five free Lie letters.

    A trivalent vertex written (a,b,c) denotes <[a,b],c>. Cyclically rotate
    the incoming parent edge to the final position, then commute the children.
    External edges are 0,...,5. All others occur exactly twice.
    """
    def expand(edge, previous):
        if edge < 6 and previous is not None:
            return {(edge,): 1}
        index, = [i for i, vertex in enumerate(vertices) if i != previous and edge in vertex]
        v = vertices[index]
        start = v.index(edge)
        left = expand(v[(start+1)%3], index)
        right = expand(v[(start+2)%3], index)
        words = defaultdict(int)
        for a, x in left.items():
            for b, y in right.items():
                words[a+b] += x*y
                words[b+a] -= x*y
        return words
    return expand(5, None)


def channels():
    """Keys record the spinors of p^3 p^4 u^5 and an ordered five-color word.

    Common units are i g^3 tau^3/(48*32^2*pi^4). Each incoming color is
    fixed by a different undifferentiated scalar trace spectator.
    """
    two, twice_one = defaultdict(int), defaultdict(int)
    spin = ((0, 1, 1), (1, 0, -1))
    for ordered in permutations((0, 1, 2)):
        i,j,k = ordered
        for first, second in ((3, 4), (4, 3)):
            sign = 1 if first == 3 else -1
            # Q0 hits either odd jet of Q2.  Bring the remaining second jet
            # past the two first jets, then put the latter in external color order.
            trees = (
                ((i,6,7),(j,7,8),(k,8,5),(6,first,second)),
                ((i,5,7),(j,7,8),(k,8,6),(6,first,second)),
            )
            weights = [rooted_words(tree) for tree in trees]
            for a,b,ab in spin:
                for h,n,hn in spin:
                    for (pa,ph,ua,ub), words, side in zip(
                            ((a,h,b,n),(b,n,a,h)), weights, (1,-1)):
                        p3,p4 = (pa,ph) if first == 3 else (ph,pa)
                        key = p3,p4,min(ua,ub),max(ua,ub)
                        for word, coefficient in words.items():
                            two[key,word] += sign*ab*hn*side*coefficient

        # The first fan uses the symmetrized pair i,j; all six ordered
        # triples above exhaust the three cyclic choices and this symmetrization.
        for d, x, y in permutations((3, 4, 5)):
            tree = ((i,6,d),(j,6,7),(7,8,x),(k,8,y))
            words = rooted_words(tree)
            sign = (-1)**sum(a>b for l,a in enumerate((d,x,y)) for b in (d,x,y)[l+1:])
            for a,b,ab in spin:
                for h,n,hn in spin:
                    jets = {d:(a,),x:(h,),y:(n,)}
                    # The input derivative lands on endpoint X or Y. Their
                    # triangle weights are 2/3 and 1/3; output 5 is the second jet.
                    if d == 5:
                        continue
                    linear = 2 if x == 5 else 1
                    saved = jets[5][0]
                    key = jets[3][0],jets[4][0],min(b,saved),max(b,saved)
                    for word, coefficient in words.items():
                        twice_one[key,word] -= 4*linear*sign*ab*hn*coefficient
    strip = lambda terms: {k:v for k,v in terms.items() if v}
    return strip(two), strip(twice_one)


if __name__ == "__main__":
    two, twice_one = channels()
    assert two == {k:-v for k,v in twice_one.items()}
    all_spins = list(product(range(2),range(2),((0,0),(0,1),(1,1))))
    per_spin = [sum(k==(a,b,*pair) for k,_ in two) for a,b,pair in all_spins]
    print("Coefficients match in each of 12 spinor projections, nonzero words:", per_spin)
    print("Distinct free words across all channels:", len({word for _,word in two}))
