"""Finite conformal integration-by-parts system for the q=6 prism weight-five coefficient."""
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, permutations, product
from math import factorial, prod
from pathlib import Path
import json

from flint import fmpq, fmpq_mat
from quartic_prism_completion import FACES, completed

Q = 6
PRIME = 1009
VERTICES = range(6)
PAIRS = tuple(combinations(VERTICES, 2))
SLOT = {pair: i for i,pair in enumerate(PAIRS)}
SLOT.update({(j,i): k for (i,j),k in tuple(SLOT.items())})
INCIDENT = [[SLOT[i,j] for j in VERTICES if j != i] for i in VERTICES]
TRIPLES = [tuple(SLOT[i,j] for i,j in combinations(t,2)) for t in combinations(VERTICES,3)]
ORDER = tuple(SLOT[i,j] for j in range(2,6) for i in range(2)) + tuple(
    SLOT[i,j] for i,j in combinations(range(2,6),2)) + (SLOT[0,1],)
UNDO_ORDER = tuple(ORDER.index(j) for j in range(15))
PERMUTATIONS = defaultdict(list)
for perm in permutations(VERTICES):
    slots = tuple(SLOT[perm[i],perm[j]] for i,j in PAIRS)
    PERMUTATIONS[slots[SLOT[0,1]]].append(tuple(slots[j] for j in ORDER))


@cache
def canonical(edges):
    # The displayed chart places a missing pair at (0,1).
    best = min(tuple(edges[j] for j in p) for i,a in enumerate(edges) if a == 0
               for p in PERMUTATIONS[i])
    return tuple(best[j] for j in UNDO_ORDER)


def converges(edges):
    return max(edges) < Q and all(sum(edges[k] for k in t) < 2*Q for t in TRIPLES)


def bubble(edges):
    return any(sum(edges[i] > 0 for i in row) <= 3 for row in INCIDENT)


def reflections(edges):
    for v in VERTICES:
        neighbors = [j for j in VERTICES if j != v and edges[SLOT[v,j]]]
        if len(neighbors) != 4 or any(edges[SLOT[v,j]] < 0 for j in neighbors):
            continue
        child = list(edges)
        ratio = Fraction(1)
        for j in neighbors:
            a = edges[SLOT[v,j]]
            child[SLOT[v,j]] = Q-a
            ratio *= Fraction(factorial(Q-a-1),factorial(a-1))
        for j,k in combinations(neighbors,2):
            child[SLOT[j,k]] += edges[SLOT[v,j]]+edges[SLOT[v,k]]-Q
        assert converges(child)
        yield canonical(tuple(child)),ratio


def component(first):
    # Entries are the factors F(first)/F(chart) of the completed cut.
    factors = {first: Fraction(1)}
    queue = [first]
    for g in queue:
        for h, ratio in reflections(g):
            weight = factors[g]*ratio
            if h not in factors:
                factors[h] = weight
                queue.append(h)
            else:
                assert factors[h] == weight
    return factors


def enumerate_charts():
    # Four equator edges meet each apex; each power is in [3-q,q-1].
    stars = [s for s in product(range(3-Q,Q),repeat=4) if sum(s)==2*Q]
    labels = {}
    representatives = []
    for a in stars:
        for b in stars:
            if list(zip(a,b)) != sorted(zip(a,b)):
                continue
            margins = [2*Q-i-j for i,j in zip(a,b)]
            for x,y in product(range(4-2*Q,Q),repeat=2):
                z = 2*Q-margins[3]-x-y
                c = (x,y,margins[0]-x-y,z,margins[1]-x-z,margins[2]-y-z)
                if max(c) >= Q or min(c) < 4-2*Q:
                    continue
                graph = (0,*a,*b,*c)
                if not converges(graph):
                    continue
                seed = canonical(graph)
                if seed in labels:
                    continue
                orbit = component(seed)
                vanishes = any(bubble(g) for g in orbit)
                index = -1 if vanishes else len(representatives)
                if not vanishes:
                    representatives.append(seed)
                for g, ratio in orbit.items():
                    assert g not in labels
                    labels[g] = index, 1/ratio
    return labels,representatives


def relations(edges, labels):
    # The two references are equator points. This keeps edge 01 absent.
    for v in VERTICES:
        finite = [j for j in range(2,6) if j != v]
        for infinity,ref in permutations(finite,2):
            candidates = [(edges,edges[SLOT[v,infinity]]-edges[SLOT[v,ref]])]
            valid = True
            for j in VERTICES:
                if j in (v,infinity,ref) or not edges[SLOT[v,j]]:
                    continue
                for u,w,sign in ((v,j,-1),(j,v,1)):
                    shifted = list(edges)
                    shifted[SLOT[v,j]] += 1
                    shifted[SLOT[ref,infinity]] += 1
                    shifted[SLOT[u,ref]] -= 1
                    shifted[SLOT[w,infinity]] -= 1
                    if not converges(shifted):
                        valid = False
                        break
                    candidates.append((tuple(shifted),sign*edges[SLOT[v,j]]))
                if not valid:
                    break
            if not valid:
                continue
            row = defaultdict(Fraction)
            for graph,weight in candidates:
                if not weight or bubble(graph):
                    continue
                idx, ratio = labels[canonical(graph)]
                if idx >= 0:
                    row[idx] += weight*ratio
            row = {i:c for i,c in row.items() if c}
            if row:
                yield row


def independent_rows(labels,count):
    pivots = {}
    selected = []
    for graph in labels:
        for row in relations(graph,labels):
            work = {i: c.numerator*pow(c.denominator,-1,PRIME) % PRIME for i,c in row.items()}
            work = {i:w for i,w in work.items() if w}
            while work:
                column = min(work)
                top = work.pop(column)
                if column not in pivots:
                    inv = pow(top,-1,PRIME)
                    pivots[column] = {j:v*inv % PRIME for j,v in work.items()}
                    selected.append(row)
                    break
                for j,v in pivots[column].items():
                    result = (work.get(j,0)-top*v) % PRIME
                    if result:
                        work[j] = result
                    else:
                        work.pop(j,None)
            if len(selected) == count-1:
                return selected
    raise ArithmeticError('Insufficient rank')


def solve():
    labels, representatives = enumerate_charts()
    size = len(representatives)
    assert (len(labels), size) == (5447,727)
    rows = independent_rows(labels,size)
    data = json.loads((Path(__file__).resolve().parents[1]/'nonSP/verify/prism_forcer/periods.json').read_text())
    measured = []
    for ex, coeff in zip(data['exponents'], data['coefficients']):
        matrix = completed(tuple(v+1 for v in ex),Q)
        g = canonical(tuple(matrix[i][j] for i,j in PAIRS))
        value = Fraction(str(coeff[2]))*factorial(Q-1)/prod(map(factorial,ex))
        measured.append((*labels[g],value))
    # One nonzero measured period supplies the scale of the homogeneous kernel.
    anchor, ratio, given = next(m for m in measured if m[2])
    matrix = fmpq_mat(size,size)
    for i,row in enumerate(rows):
        for j,c in row.items():
            matrix[i,j] = fmpq(c.numerator,c.denominator)
    matrix[size-1,anchor] = 1
    rhs = fmpq_mat(size,1)
    rhs[size-1,0] = fmpq(given.numerator*ratio.denominator,given.denominator*ratio.numerator)
    vector = matrix.solve(rhs)
    values = [Fraction(str(vector[i,0])) for i in range(size)]
    assert all((values[idx]*factor if idx >= 0 else 0)==observed for idx,factor,observed in measured)
    # Prove compatibility with every generated rational row, including unselected identities.
    for graph in labels:
        for row in relations(graph,labels):
            assert not sum(c*values[i] for i,c in row.items())
    print('5447 completed charts, 727 reflected unknowns, rank 726 mod 1009;')
    print('one exact normalization reproduces all 1481 integrated weight-five coefficients.')
    print('Nonzero components:',sum(bool(v) for v in values),'anchor:',given)


if __name__ == '__main__':
    solve()
