"""Two true Betti factors in the fourfold D3 bound, with stable-graph multiplicities."""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations_with_replacement, product
from math import comb, factorial, prod

from soft_quartic_psi_bound import intersection, cut_intersection, upper_bound

LINKS = ((0,1,4),(1,4,2),(2,4,4),(3,4,4))


def split(pair, part):
    return (pair[0] in part) != (pair[1] in part)


def components(genus_sum, slots, last):
    """A component dimension fixes its genus from the cotangent exponents."""
    genera = []
    for degrees in slots:
        deficit = sum(degrees) + 3 - len(degrees)
        if deficit < 0 or deficit % 3:
            return 0
        genera.append(deficit//3)
    if sum(genera) != genus_sum:
        return 0
    return prod(intersection(g, ds, last) for g,ds in zip(genera,slots))


@cache
def graph_terms(genus, degrees, first, second, last):
    """Common one-edge normal term, ordered two-edge trees, and parallel-edge graph."""
    n = len(degrees)
    normal = tree = parallel = Fraction(0)
    for colors in product(range(2),repeat=n):
        left = {i for i,label in enumerate(colors) if label==0}
        if first[0] not in left or not split(first,left) or not split(second,left):
            continue
        a,b = ([degrees[i] for i in range(n) if colors[i]==label] for label in (0,1))
        normal -= components(genus,(a+[1],b+[0]),last) + components(genus,(a+[0],b+[1]),last)
        parallel += Fraction(1,2)*components(genus-1,(a+[0,0],b+[0,0]),last)
    for colors in product(range(3),repeat=n):
        left,middle,right = ({i for i,c in enumerate(colors) if c==v} for v in range(3))
        if split(first,left) and split(second,right):
            tree += components(genus,([degrees[i] for i in left]+[0],
                                      [degrees[i] for i in middle]+[0,0],
                                      [degrees[i] for i in right]+[0]),last)
    return normal,tree,parallel


def add_at(degrees,*marks):
    result = list(degrees)
    for i in marks:
        result[i] += 1
    return result


def four_times_betti_pair(genus,degrees,p,q,last):
    raw = sum(intersection(genus,add_at(degrees,i,j),last) for i in p for j in q)
    first_cut = sum(cut_intersection(genus,add_at(degrees,i),*q,last) for i in p)
    second_cut = sum(cut_intersection(genus,add_at(degrees,j),*p,last) for j in q)
    normal,tree,parallel = graph_terms(genus,tuple(degrees),p,q,last)
    return raw-first_cut-second_cut+normal+tree-parallel/3


def descendant_vertex(genus,momenta,degrees,last=False):
    """Li--Wang's full Ng=2 coefficient on conserved one- or two-dimensional momenta."""
    n = len(degrees)
    if not n or len(momenta)!=n or genus*len(momenta[0])!=2:
        raise ValueError('Expected genus two in dimension one or genus one in dimension two.')
    if any(len(row)!=len(momenta[0]) for row in momenta) or any(map(sum,zip(*momenta))):
        raise ValueError('The equal-dimensional momenta must be conserved.')
    if any(type(k) is not int or k<0 for k in degrees):
        raise ValueError('Descendant levels must be nonnegative integers.')
    if sum(degrees)!=3*genus-3+n-2:
        return Fraction(0)
    pairs = {(i,j):-sum(x*y for x,y in zip(momenta[i],momenta[j]))
             for i in range(n) for j in range(i+1,n)}
    return sum(Fraction(1 if p==q else 2,8)*pairs[p]*pairs[q]*
               four_times_betti_pair(genus,degrees,p,q,last)
               for p,q in combinations_with_replacement(pairs,2) if pairs[p] and pairs[q])


def pixton_quartic(charges, degrees):
    """Charge-degree four of the genus-two graph formula, with the 1/4 genus factor."""
    n = len(charges)
    assert len(degrees) == n and sum(charges) == 0 and sum(degrees) == n + 1
    result = sum(Fraction(x*x*y*y, 8) * intersection(2, add_at(degrees, i, j))
                 for i, x in enumerate(charges) for j, y in enumerate(charges))
    for colors in product(range(2), repeat=n):
        if colors[0]:
            continue
        charge = sum(x for x, side in zip(charges, colors) if side == 0)
        if not charge:
            continue
        sides = [[d for d, side in zip(degrees, colors) if side == vertex]
                 for vertex in (0, 1)]
        for i, x in enumerate(charges):
            legs = add_at(degrees, i)
            slots = [[d for d, side in zip(legs, colors) if side == vertex] + [0]
                     for vertex in (0, 1)]
            result -= Fraction(charge**2 * x*x, 4) * components(2, slots, False)
        for vertex in (0, 1):
            slots = [side + [int(j == vertex)] for j, side in enumerate(sides)]
            result -= Fraction(charge**4, 8) * components(2, slots, False)
        result -= Fraction(charge**4, 48) * components(1, [s + [0, 0] for s in sides], False)
    for colors in product(range(3), repeat=n):
        ends = [sum(x for x, side in zip(charges, colors) if side == vertex)
                for vertex in (0, 2)]
        if all(ends):
            slots = [[d for d, side in zip(degrees, colors) if side == vertex]
                     + [0] * (2 if vertex == 1 else 1) for vertex in range(3)]
            result += Fraction(prod(ends)**2, 8) * components(2, slots, False)
    return result


def check_pixton():
    for modulus in range(2, 24):
        for charge in range(modulus):
            weights = [u * (modulus-u) for u in range(modulus)]
            convolution = sum(w * weights[(charge-u) % modulus] for u, w in enumerate(weights))
            shift = charge * (modulus-charge)
            assert Fraction(convolution, modulus) == Fraction(modulus**4-1, 30) - Fraction(shift*(shift+2), 6)
    checked = 0
    for charges in ((1,-1), (2,-3,1), (0,2,-1,-1), (1,3,-2,-2)):
        n = len(charges)
        for degrees in product(range(n+2), repeat=n):
            if sum(degrees) != n+1:
                continue
            action = descendant_vertex(2, [(x,) for x in charges], degrees)
            assert pixton_quartic(charges, degrees) == action
            if max(degrees) == n+1:
                others = [x for x, d in zip(charges, degrees) if not d]
                expected = sum(Fraction(x**4, 1920) for x in others)
                expected += sum(Fraction(x*x*y*y, 576) for i,x in enumerate(others) for y in others[i+1:])
                assert action == expected
            checked += 1
    print('Pixton quartic matches the two-Betti vertex on', checked, 'charge and descendant fixtures.')


def bounds(links,last):
    multiplicities = Counter({(i,j):m for i,j,m in links})
    divisor = 2**14 * prod(factorial(m) for m in multiplicities.values())
    values = {}
    for p,q in combinations_with_replacement(multiplicities,2):
        powers = multiplicities.copy()
        powers.subtract((p,q))
        if min(powers.values()) < 0:
            continue
        total = 0
        for parts in product(*(range(n+1) for n in powers.values())):
            degrees = [0]*5
            weight = 1
            for ((i,j),n),k in zip(powers.items(),parts):
                degrees[i] += k
                degrees[j] += n-k
                weight *= comb(n,k)
            total += weight * four_times_betti_pair(4,degrees,p,q,last)
        values[p,q] = 2*total/divisor
    return values


def check():
    check_pixton()
    for last in (False,True):
        # On a genus-zero curve every Abel map is constant. A genus-one single-pair
        # form has complex rank one even as the curve and all its marks vary.
        for pair in ((0,1),(0,2),(3,4)):
            assert four_times_betti_pair(0,(0,)*5,(0,1),pair,last)==0
        normal,tree,parallel = graph_terms(1,(0,0),(0,1),(0,1),last)
        assert (normal,tree,parallel)==(0,0,Fraction(1,2))
        for n in (2,3,4):
            for ds in product(range(n-1),repeat=n):
                if sum(ds)==n-2:
                    assert four_times_betti_pair(1,ds,(0,1),(0,1),last)==0
        for ds in ((1,0,0),(0,1,0),(0,0,1)):
            assert descendant_vertex(1,((1,0),(0,1),(-1,-1)),ds,last)==Fraction(1,24)
        assert descendant_vertex(2,((1,),(-1,)),(3,0),last)==Fraction(1,1920)
        assert descendant_vertex(2,((1,),(-1,)),(2,1),last)==Fraction(1,640)
        result=bounds(LINKS,last)
        for (p,q),value in result.items():
            one = [upper_bound(tuple(sorted(LINKS,key=lambda link:link[:2]!=s)),last)[1]
                   for s in (p,q)]
            assert 0 < value <= min(one)
        if last:
            assert result==control
        else:
            control=result
    print('Two retained Betti factors, zero-based marks:')
    for (p,q),value in result.items():
        print(p,q,value)
    print('Minimum',min(result.values()))


if __name__=='__main__':
    check()
