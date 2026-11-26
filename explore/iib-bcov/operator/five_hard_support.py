"""Exact support tests for the five-site mixed four-loop Ward integrand."""
from itertools import combinations
from random import Random

import networkx as nx
from flint import fmpz_mat, nmod_mat

from fourfold_quartic_orbits import label_images, oriented_symmetries
from five_hard_mixed_double import HARD, NORMAL, four_voltage, hard_order
from quartic_seven_plane_lift import G, cubic_cores, fixture
from quartic_seven_plane_primitive import admissible


def restrictions(core):
    """Use only cycles supported on the chosen paths, including the entire core."""
    out = []
    for size in range(2, 10):
        for subset in combinations(range(9), size):
            graph = nx.MultiGraph(core[e] for e in subset)
            if not nx.is_connected(graph) or nx.has_bridges(graph):
                continue
            incidence = fmpz_mat([[int(core[e][0] == v)-int(core[e][1] == v)
                                   for e in subset] for v in graph])
            cycles, loops = incidence.nullspace()
            columns = {e: tuple(int(cycles[i, j]) for j in range(loops))
                       for i, e in enumerate(subset)}
            exposed = set(graph) & {v for e in range(9) if e not in columns for v in core[e]}
            out.append((columns, 3*loops + 1-size, max(0, len(exposed)-1)))
    return out


def excluded(marks, order, tests):
    for columns, need, boundary in tests:
        sources = [(columns[e], NORMAL[j]) for e, j in zip(marks, order) if e in columns]
        if need > 0:
            jumps = [[a*b for a in cycle for b in wave] for cycle, wave in sources]
            if not jumps or fmpz_mat(jumps).rank() < need:
                return 'jumps'
        if boundary < 3:
            local = [wave for _, wave in sources]
            if (fmpz_mat(local).rank() if local else 0) + boundary < 3:
                return 'boundary'
    return None


def metric(core, marks, times, prime):
    edges = fixture(core, marks)
    cycles = G(11, edges).C
    incidence = nmod_mat([[int(a == i)-int(b == i) for i in range(1, 11)]
                          for a, b in edges], prime)
    weighted = nmod_mat([[int(incidence[j,i])*pow(t, -1, prime) for i in range(10)]
                         for j, t in enumerate(times)], prime)
    laplacian = incidence.transpose()*weighted
    if not laplacian.det() or not sum(times) % prime:
        return None
    green = laplacian.inv()
    transfer = weighted*green
    hard = [[int(transfer[i,j]) for j in range(5,10)] for i in range(14)]
    def resistance(i, j):
        return int(green[i-1,i-1]+green[j-1,j-1]-2*green[i-1,j-1])
    gram = nmod_mat([[sum(t*a*b for t,a,b in zip(times,c,d)) for c in cycles]
                     for d in cycles], prime)
    return cycles, hard, resistance, gram.inv()*gram.det()


def ranks_and_density(network, pair, prime):
    cycles, transfer, resistance, adjugate = network
    order = hard_order(pair)
    def block(waves):
        values = [sum(row[k]*waves[j] for k,j in enumerate(order)) % prime for row in transfer]
        return [[v*c for v,c in zip(values,cycle)] for cycle in cycles]
    quotient = [[1]*14]
    for wave in zip(*NORMAL):
        quotient.extend(block(wave))
    projected = nmod_mat(quotient, prime).rank()
    extra = block([p[1] for p in HARD])
    full = nmod_mat(quotient + extra, prime).rank()
    voltage = four_voltage(resistance, pair) % prime
    density = 0
    if full == 14 and projected == 13 and voltage:
        d = nmod_mat([[int(nmod_mat(quotient + [row], prime).det()) for row in extra]], prime)
        density = int((d*adjugate*d.transpose())[0,0])*voltage % prime
    # For all ten pairs the omitted projective factors and label multiplicity agree.
    return full, projected, voltage, density


def color_orbits(marks, group):
    left = set(combinations(range(5), 2))
    while left:
        pair = min(left)
        word = tuple(int(i in pair) for i in range(5))
        words = label_images(word, marks, group)
        pairs = {tuple(i for i, v in enumerate(w) if v) for w in words}
        left -= pairs
        yield pair, len(pairs)


def census():
    prime = 1000003
    rng = Random(13349)
    for index, core in enumerate(cubic_cores()):
        tests = restrictions(core)
        actions = tuple(oriented_symmetries(core))
        counts = {}
        shapes = 0
        labeled_shapes = 0
        unresolved = []
        for marks, size in admissible(core):
            group = [a for a in actions if tuple(sorted(a[e][0] for e in marks)) == marks]
            assert len(group)*size == len(actions)
            fixtures = [metric(core, marks, rng.sample(range(2, 233), 14), prime) for _ in range(3)]
            tests_by_pair = {pair: [ranks_and_density(n,pair,prime) for n in fixtures if n is not None]
                             for pair in combinations(range(5),2)}
            # Add the ten actual pairs at the same fixed metric, before taking shape orbits.
            assert any(sum(x[3] for x in col) % prime for col in zip(*tests_by_pair.values()))
            shapes += 1
            labeled_shapes += size
            for pair, colors in color_orbits(marks, group):
                reason = excluded(marks, hard_order(pair), tests)
                observed = tests_by_pair[pair]
                witness = any(f == 14 and q == 13 and v for f,q,v,_ in observed)
                assert not (reason and witness), (index, marks, pair, reason)
                if witness:
                    reason = 'witness'
                if not reason:
                    unresolved.append((marks, pair, observed[0]))
                n, mass = counts.get(reason,(0,0))
                counts[reason] = n+1, mass+12*colors*size
        assert not unresolved, (index, unresolved)
        assert (shapes, labeled_shapes) == [(11,69),(0,0),(36,171),(80,786),(26,1134)][index]
        assert counts.get('witness',(0,0)) == [(47,6912),(0,0),(217,17820),(675,92484),(177,132840)][index]
        assert sum(n for _,n in counts.values()) == 120*labeled_shapes
        print(index, 'hard shapes and placements', shapes, labeled_shapes,
              'color orbits and labeled orders', counts, flush=True)


if __name__ == '__main__':
    census()
