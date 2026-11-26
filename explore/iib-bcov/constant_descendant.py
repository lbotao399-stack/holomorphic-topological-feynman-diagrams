"""Constant first descendant on flat fivefold primary point graphs."""
from fractions import Fraction
from itertools import product
from math import comb, factorial, prod

from longitudinal_theta import primary_coefficient


def labeled_attachments(weights, count):
    """Attach each labeled zero wave once, summing over all choices of action vertex."""
    states = {tuple(0 for _ in weights): 1}
    for _ in range(count):
        following = {}
        for occupancy, multiplicity in states.items():
            for vertex, weight in enumerate(weights):
                shifted = list(occupancy)
                shifted[vertex] += 1
                key = tuple(shifted)
                following[key] = following.get(key, 0) + multiplicity*(weight+occupancy[vertex])
        states = following
    return sum(states.values())


def check():
    # The geometrical loops have 4h+1 action vertices; vertex genera are additional.
    for loops in range(3):
        vertices = 4*loops+1
        edges = 5*loops+1
        for leading in product(range(3), repeat=min(vertices,2)):
            genera = [*leading, *([0]*(vertices-len(leading)))]
            genus = loops+sum(genera)
            valences = [2*g+3 for g in genera]
            weights = [2*g-2+l for g,l in zip(genera,valences)]
            for spokes in range(1,loops+2):
                external = sum(valences)-2*edges+spokes
                assert external == 2*genus+1+spokes
                assert sum(weights) == 2*genus+external-spokes == 4*genus+1
            for count in range(5):
                assert labeled_attachments(weights,count) == factorial(count)*comb(4*genus+count,count)
    # Genus zero: (m+3)-field vertex with three primary inputs has coefficient m!.
    for count in range(12):
        assert labeled_attachments([1],count) == factorial(count)
    # The classical intersection formula also allows nonzero hard descendant levels.
    for levels in ((1,),(2,1),(1,1,1,2),(4,3,2,1)):
        total_level = sum(levels)
        denominator = prod(factorial(k) for k in levels)
        hard = Fraction(factorial(total_level),denominator)
        for count in range(8):
            full = Fraction(factorial(total_level+count),denominator)
            assert full == hard*labeled_attachments([1+total_level],count)
        # Keep at least three hard marks. Extra level-zero hard slots pad the
        # list to the stable arity 3+total_level-primary_count.
        for primary_count in range(min(3,total_level,3+total_level-len(levels))+1):
            lowered = Fraction(0)
            for drops in product(*(range(k+1) for k in levels)):
                if sum(drops) == primary_count:
                    assignments = factorial(primary_count)//prod(factorial(d) for d in drops)
                    lowered += assignments*Fraction(factorial(total_level-primary_count),
                                                    prod(factorial(k-d) for k,d in zip(levels,drops)))
            for dilatons in range(4):
                direct = Fraction(factorial(total_level+dilatons),denominator)
                assert lowered*labeled_attachments([1+total_level],dilatons) == direct
    # Ward degree with hard descendants distributed among a one-loop graph's
    # five action vertices, two of which carry genus-one local actions.
    genera = (0,1,0,1,0)
    for levels in product(range(3),repeat=5):
        weights = [4*g+1+d for g,d in zip(genera,levels)]
        total_genus = 1+sum(genera)
        assert labeled_attachments(weights,3) == 6*comb(4*total_genus+sum(levels)+3,3)
    frame = [tuple(int(i==j) for i in range(5)) for j in range(5)]
    elementary = primary_coefficient([frame[i] for i in (0,0,1,2,3)],frame[4],(0,0,0,-1,0))
    assert elementary == -Fraction(1,31104)
    # Compare two geometric loops with the single genus-two action vertex.
    for count in range(5):
        assert labeled_attachments([1]*9,count) == labeled_attachments([9],count)
    assert elementary * labeled_attachments([9],1) == -Fraction(1,3456)
    assert elementary * labeled_attachments([9],2) == -Fraction(5,1728)
    # Remove the soft primary first: it can lower the hard level-one scalar
    # (leaving nine dilaton attachments), or it can lower the soft dilaton.
    assert elementary*labeled_attachments([10],1) == elementary*(9+1)
    print('Two-loop elementary scalar: reduced coefficients with zero, one, two labeled dilatons:',
          elementary, elementary*9, elementary*90)


if __name__ == '__main__':
    check()
