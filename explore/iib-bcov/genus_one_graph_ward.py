"""Recover flat fivefold genus-one trace vertices from the one-loop graph Ward equation."""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, prod

from flint import fmpq_mat


def plus(momenta):
    return tuple(map(sum, zip(*momenta)))


@lru_cache(maxsize=4096)
def square(momenta):
    return Fraction(str(fmpq_mat([[str(x) for x in row] for row in momenta]).det()))**2


@lru_cache(maxsize=1)
def circle_green():
    """Integrate the grounded resistance Green matrix and sum the labeled hexagons."""
    # A grounded six-edge circle has G_ij=min(s_i,s_j)-s_i*s_j.
    # For unit-simplex lengths the means are i/6 and i*(j+1)/42, i <= j.
    orders = list(permutations(range(6)))
    normalization = 12*factorial(5)
    def mean(i, j):
        return Fraction(min(i,j)*(6-max(i,j)),42)
    diagonal = sum(mean(order[0],order[0]) for order in orders)/normalization
    off_diagonal = sum(mean(order[0],order[1]) for order in orders)/normalization
    # Momentum conservation eliminates the common off-diagonal term.
    return diagonal-off_diagonal


def circle_groups(remaining, levels, count):
    if not remaining:
        if not count:
            yield ()
        return
    if not count:
        return
    first, *rest = remaining
    for size in range(len(remaining)-count+1):
        for extra in combinations(rest, size):
            group = (first, *extra)
            if sum(levels[i] for i in group) != len(group)-1:
                continue
            others = tuple(i for i in rest if i not in extra)
            for partition in circle_groups(others,levels,count-1):
                yield (group, *partition)


def classical(levels, group):
    degrees = [levels[i] for i in group]
    return Fraction(factorial(sum(degrees)), prod(factorial(k) for k in degrees))


def graph_terms(momenta, levels):
    """The incidence potentials at individual waves, with one common origin per graph."""
    count = len(levels)
    labels = tuple(range(count))
    hexagon = [[Fraction(0)]*5 for _ in labels]
    tree = [[Fraction(0)]*5 for _ in labels]
    for groups in circle_groups(labels,levels,6):
        charges = [plus(momenta[i] for i in group) for group in groups]
        weight = (circle_green()*square(tuple(charges[:5]))
                  * prod(classical(levels,group) for group in groups))
        for group, charge in zip(groups,charges):
            for slot in group:
                for j in range(5):
                    hexagon[slot][j] += weight*charge[j]
    # The other end of the single tree edge is the root. Its voltage is zero.
    for size in range(2,count-3):
        for group in combinations(labels,size):
            if sum(levels[i] for i in group) != size-2:
                continue
            charge = plus(momenta[i] for i in group)
            outside = [i for i in labels if i not in group]
            fused = tuple(momenta[i] for i in outside)+(charge,)
            fused_levels = tuple(levels[i] for i in outside)+(0,)
            weight = classical(levels,group)*graph_vertex(fused,fused_levels)
            for slot in group:
                for j in range(5):
                    tree[slot][j] += weight*charge[j]
    return hexagon, tree


def separate_divergences(momenta, potentials):
    """Find the uniquely determined slot divergences on nonzero, jointly noncollinear waves."""
    for i, j in combinations(range(len(momenta)),2):
        for x, y in combinations(range(5),2):
            det = momenta[i][x]*momenta[j][y]-momenta[i][y]*momenta[j][x]
            if not det:
                continue
            dx, dy = (potentials[i][a]-potentials[j][a] for a in (x,y))
            value = Fraction(dx*momenta[j][y]-dy*momenta[j][x],det)
            common = [entry-value*p for entry,p in zip(potentials[i],momenta[i])]
            values = []
            for row, p in zip(potentials,momenta):
                axis = next((a for a in range(5) if p[a]),None)
                if axis is None:
                    raise ValueError('At zero momentum the divergence cannot determine this slot.')
                coefficient = Fraction(row[axis]-common[axis],p[axis])
                assert all(a-c==coefficient*b for a,b,c in zip(row,p,common))
                values.append(coefficient)
            return tuple(values)
    raise ValueError('Two nonparallel wave momenta are required to fix the common potential.')


@lru_cache(maxsize=4096)
def recovered(momenta, levels):
    if len(momenta)!=len(levels) or sum(levels)!=len(levels)-6:
        raise ValueError('The fivefold defect at n waves has total descendant level n-6.')
    circle, tree = graph_terms(momenta,levels)
    return separate_divergences(momenta,[plus(pair) for pair in zip(circle,tree)])


@lru_cache(maxsize=4096)
def graph_vertex(momenta, levels):
    if not any(levels):
        # The momentum determinant of five conserved primary waves vanishes.
        return Fraction(0)
    i = next(i for i,k in enumerate(levels) if k)
    lowered = tuple(k-int(j==i) for j,k in enumerate(levels))
    return recovered(momenta,lowered)[i]


def seven_wave_formula(momenta, d):
    pairs = {(i,j):square(tuple(p for s,p in enumerate(momenta) if s not in (i,j)))
             for i,j in combinations(range(7),2)}
    totals = [sum(b for group,b in pairs.items() if i in group) for i in range(7)]
    common = [tuple(sum(b*momenta[j if i==r else i][axis] for (i,j),b in pairs.items()
                        if r in (i,j)) for axis in range(5)) for r in range(7)]
    assert len(set(common))==1
    return tuple((totals[d] if r==d else totals[d]+totals[r]-pairs[tuple(sorted((d,r)))])/24
                 for r in range(7))


def replay():
    from elliptic_descendants import vertex
    from random import Random
    random = Random(774211)
    assert circle_green()==Fraction(1,24)
    for n in range(6,10):
        independent_waves = [tuple(random.randrange(-4,5) for _ in range(5)) for i in range(n-1)]
        momenta = tuple(independent_waves)+ (tuple(-x for x in plus(independent_waves)),)
        assert all(any(q) for q in momenta)
        cases = 0
        for slots in combinations_with_replacement(range(n),n-6):
            levels = tuple(slots.count(i) for i in range(n))
            graph = recovered(momenta,levels)
            if n==7:
                assert graph==seven_wave_formula(momenta,slots[0])
            for i,value in enumerate(graph):
                raised = tuple(k+int(j==i) for j,k in enumerate(levels))
                assert value==vertex(momenta,raised), (levels,i,value)
                cases += 1
        print(n,'waves:',cases,'labeled free-divergence coefficients verified.',flush=True)
    # Five basis waves, one generic wave, and their balancing wave.
    basis = tuple(tuple(int(i==j) for j in range(5)) for i in range(5))
    extra = (1,2,3,4,5)
    points = (*basis,extra,tuple(-x for x in plus((*basis,extra))))
    levels = (0,0,0,0,0,1,0)
    graph = recovered(points,levels)
    assert graph==seven_wave_formula(points,5)
    print('Seven-wave example with descendant at (1,2,3,4,5):',*graph)


if __name__=='__main__':
    replay()
