"""Paired planar momenta in the difference of two rooted fivefold descendant measurements."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

IIB = Path(__file__).resolve().parents[2]/'iib-bcov'
sys.path.insert(0, str(IIB))
from elliptic_descendants import vertex as elliptic
from genus_one_graph_ward import graph_vertex


def advance(states, axis, required):
    """One labeled axis mark is left free or fused into a retained positive-level scalar."""
    result = Counter() if required else states.copy()
    for (first,second,levels),weight in states.items():
        for i,level in enumerate(levels):
            if not level:
                continue
            directions = [list(first),list(second)]
            directions[i//2][axis] += 1
            after = tuple(k-int(j==i) for j,k in enumerate(levels))
            result[(*map(tuple,directions),after)] -= weight
    return result


def canonical(first, second, levels):
    # Swap axial coordinates e3 and e4: they have the same momentum multiplicities.
    swapped = (first[1],first[0],first[2]),(second[1],second[0],second[2])
    first,second = min((first,second), swapped)
    return first,second,tuple(sorted(levels[:2])),tuple(sorted(levels[2:]))


def combinations_for(genus, first_levels, forced=False, probes=(((1,1),1),((2,0),-2))):
    zero = (0,0,0)
    all_terms = Counter()
    for last_levels,sign in probes:
        states = Counter({(zero,zero,(*first_levels,*last_levels)):1})
        if forced:
            # Include the root and vector mark in the forgettable set. Each must enter a tail.
            for axis in (0,1):
                states = advance(states,axis,True)
        for axis in range(3):
            for _ in range(genus):
                states = advance(states,axis,False)
        for (first,second,levels),weight in states.items():
            all_terms[canonical(first,second,levels)] += sign*weight
    return {key:weight for key,weight in all_terms.items() if weight}


def differences(terms):
    result = {}
    for (first,second,k,h),weight in terms.items():
        assert second==(0,0,0) and h in ((0,2),(1,1))
        if h==(1,1):
            assert terms[first,second,k,(0,2)] == -2*weight
            result[first,k] = weight
    return result


def waves(genus,first,second,levels):
    p = [(1,0,*[x-y for x,y in zip(first,(genus+1,genus+1,genus))]),
         (-1,0,0,0,0),(0,1,*second),(0,-1,0,0,0)]
    n = (genus+1,genus+1,genus)
    for i in range(3):
        p += [tuple(int(j==i+2) for j in range(5))]*(n[i]-first[i]-second[i])
    assert not any(map(sum,zip(*p))) and len(p)==2*genus+3+sum(levels)
    return tuple(p), (*levels,*(0 for _ in p[4:]))


@lru_cache(maxsize=256)
def elliptic_graph(first,second,k,h):
    p,levels = waves(1,first,second,(*k,*h))
    value = graph_vertex(p,levels)
    assert value == elliptic(p,levels)
    return value


def one_loop_check():
    expanded = {}
    # Enumerate individual donor maps before using either opposite-pair or axial symmetry.
    for levels in ((1,1,1,1),(1,1,2,0)):
        expanded[levels[2:]] = Fraction()
        for assignment in product(range(5),repeat=3):
            used = tuple(assignment.count(i) for i in range(4))
            if any(i>k for i,k in zip(used,levels)):
                continue
            first = tuple(int(site in (0,1)) for site in assignment)
            second = tuple(int(site in (2,3)) for site in assignment)
            remain = tuple(k-i for k,i in zip(levels,used))
            p,desc = waves(1,first,second,remain)
            # Undo the shear separately on both pairs. Each negative planar leg receives
            # exactly the subset of the three independent labeled donors prescribed above.
            p = [list(row) for row in p]
            for axis,i in enumerate(assignment):
                if i in (1,3):
                    p[i-1][axis+2] -= 1
                    p[i][axis+2] += 1
            value = elliptic(tuple(map(tuple,p)),desc)
            expanded[levels[2:]] += (-1)**sum(used)*value
    for forced in (False,True):
        for last,wanted in expanded.items():
            table = combinations_for(1,(1,1),forced,((last,1),))
            assert sum(w*elliptic_graph(*key) for key,w in table.items()) == wanted
    print('Genus-one direct labeled polarization and hexagon/tree Ward reduction:',expanded)


if __name__=='__main__':
    direct = differences(combinations_for(3,(2,2)))
    forced = differences(combinations_for(3,(2,2),True))
    assert len(direct)==24 and len(forced)==8
    for level in range(7):
        row = [(k,w) for k,w in direct.items() if sum(k[0])==level]
        assert sum(abs(w) for _,w in row)==(1,18,144,504,756,0,0)[level]
    print('Eight rank-reduced fivefold differences:')
    for (axes,levels),weight in sorted(forced.items()):
        print(axes,levels,weight)
    one_loop_check()


def axial_fiber_degree(vectors):
    """Three Hodge modes per axial coordinate. Integrate the nine point positions."""
    from itertools import combinations
    if len(vectors)!=9 or any(len(v)!=3 for v in vectors):
        raise ValueError('Nine three-component transverse momenta required.')
    def det(three):
        a,b,c=(vectors[i] for i in three)
        return (a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])
                +a[2]*(b[0]*c[1]-b[1]*c[0]))
    triples={frozenset(t):det(t)**2 for t in combinations(range(9),3)}
    labels=frozenset(range(9))
    return sum(v*w*triples[labels-a-b] for a,v in triples.items() if v
               for b,w in triples.items() if w and not a&b)


def elementary_source(genus, counts):
    a,b,c=counts
    assert a+b+c==2*genus+1
    axes=((0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1))
    root=(0,0,3-a,2-b,3-c)
    p=((1,0,-3,-3,-3),(-1,0,0,0,0),(0,1,0,0,0),(0,-1,0,0,0),
       *(axes[0],)*a,*(axes[1],)*(b+1),*(axes[2],)*c,root)
    assert not any(map(sum,zip(*p)))
    return p


def ordinary_polarization_check(two_loop=False):
    from itertools import combinations
    from math import factorial
    from mixed_sextic_vertex import labeled_sum,sextic
    planar=((1,0),(-1,0),(0,1),(0,-1))
    periods=[]
    for first in ((0,2),(1,1)):
        candidates=[]
        for second in ((1,1),(2,0)):
            levels=(*first,*second)
            value=sextic(3,planar,levels)
            assert value==labeled_sum(3,planar,levels)
            candidates.append(value)
        periods.append(tuple(216*(x-2*y) for x,y in zip(*candidates)))
    assert periods==[(Fraction(253,168000),0),(-Fraction(7631,2910600),Fraction(48,13475))]
    ordinary=tuple(y-2*x for x,y in zip(*periods))
    assert ordinary==(-Fraction(327949,58212000),Fraction(48,13475))
    reduced=differences(combinations_for(3,(2,2),True))
    positions=Counter()
    for (abc,levels),weight in reduced.items():
        n=tuple(sorted((4-abc[0],4-abc[1],3-abc[2]),reverse=True))
        positions[n,levels] +=weight
    assert sorted(positions.values())==[-54,2,2,54,162]
    for n in range(4):
        for b in range(4-n):
            p=elementary_source(1,(n,b,3-n-b))
            value=sum((1 if i==1 else -2)*(1 if j==1 else -2)*graph_vertex(
                p,(i,2-i,j,2-j,*(0 for _ in p[4:]))) for i,j in product((0,1),repeat=2))
            assert value==0
    if two_loop:
        from genus_two_graph_ward import vertex
        p=elementary_source(2,(2,1,2))
        values=[vertex(p,(i,2-i,j,2-j,*(0 for _ in p[4:]))) for i,j in product((0,1),repeat=2)]
        assert 4*values[0]-2*values[1]-2*values[2]+values[3]==0
        print('Individual two-loop Ward terms:',*values,flush=True)
    total=0
    for a in range(8):
        for b in range(8-a):
            c=7-a-b
            p=elementary_source(3,(a,b,c));root=p[-1]
            degree=axial_fiber_degree([w[2:] for w in p[4:]])
            assert degree==(216 if (a,b,c) in ((2,2,3),(3,1,3),(3,2,2)) else 0)
            total+=Fraction(root[2]*degree,factorial(a)*factorial(b)*factorial(c))
    assert total==9
    # Independent coefficient of det(sum t_i w_i w_i^T)^3 by binary finite differences.
    # The homogeneous polynomial has degree nine, so only its squarefree term survives.
    from flint import fmpz_mat
    probe=((1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1),(1,0,1),(2,1,1),(-1,2,1),(1,1,-2))
    prefactor=0
    for k in range(10):
        for choose in combinations(probe,k):
            gram=fmpz_mat([[sum(w[i]*w[j] for w in choose) for j in range(3)] for i in range(3)])
            prefactor+=(-1)**(9-k)*int(gram.det())**3
    assert prefactor==axial_fiber_degree(probe)
    print('Two independent relations:',periods,'double difference:',ordinary)
    print('First surviving fixed Fourier coefficient:',*(total*x/216 for x in ordinary))


if __name__=='__main__':
    ordinary_polarization_check('--two-loop' in sys.argv)
