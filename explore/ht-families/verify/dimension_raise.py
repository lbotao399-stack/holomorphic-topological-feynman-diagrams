"""Signed cotree subdivision in the project's canonical graph and projective orientation."""
from fractions import Fraction
from itertools import combinations
from math import prod
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "code"))
from htk import Graph, brute, fan, polygon, theta
from htk.linalg import det


def subdivide(graph, selected):
    rank = {edge:j for j,edge in enumerate(selected)}
    root = graph.root + len(selected)

    def old_vertex(v):
        return root if v == graph.root else v

    initial = []
    appended = []
    for edge,(tail,head) in enumerate(graph.edges):
        tail,head = old_vertex(tail),old_vertex(head)
        if edge in rank:
            vertex = graph.root + rank[edge]
            initial.append((tail,vertex))
            appended.append((vertex,head))
        else:
            initial.append((tail,head))
    return Graph(root+1,initial+appended)


def split_times(times, selected, ratios):
    result = list(times)
    for edge,ratio in zip(selected,ratios):
        result[edge] = (1-ratio)*times[edge]
        result.append(ratio*times[edge])
    return result


def check(graph, holomorphic, topological, seed):
    rng = Random(seed)
    n = graph.loops
    times = [Fraction(rng.randrange(2,7),rng.randrange(2,7)) for _ in graph.edges]
    cycle = graph.cycle_matrix()
    sources = [[rng.randrange(-9,10) for _ in range(holomorphic)] for _ in range(graph.root)]
    original = brute.density(graph,times,sources,T=topological)
    assert original.coeff
    period = brute.symanzik_U(graph,times,cycle)
    lifted_sources = [row+[0] for row in sources] + [[0]*holomorphic+[1] for _ in range(n)]
    adjustment = (-1)**(n*(holomorphic+topological+1))
    total = Fraction(0)
    signed_minors = set()
    count = 0
    for selected in combinations(range(graph.ne),n):
        minor = det([[row[i] for i in selected] for row in cycle])
        if not minor:
            continue
        signed_minors.add(minor)
        raised = subdivide(graph,selected)
        tree = list(range(graph.ne))
        assert raised.sigma(tree) == (-1)**n * minor
        expected = original.scaled(adjustment*minor/period)
        for trial in range(2):
            ratios = [Fraction(rng.randrange(1,9),9) for _ in selected]
            values = split_times(times,selected,ratios)
            value = brute.density(raised,values,lifted_sources,T=topological)
            assert value.same_as(expected)
        total += adjustment*minor*(value.coeff/original.coeff)*prod(times[i] for i in selected)
        count += 1
    assert total == 1
    print(f"H={holomorphic}, T={topological}, n={n}: {count} cotrees, "
          f"minor signs {sorted(signed_minors)}, basis factor {adjustment}, signed sum = 1")


if __name__ == "__main__":
    check(polygon(3),2,0,3101)
    check(polygon(4),3,0,3102)
    check(fan(3),2,0,3103)
    check(theta(1,2,2)[0],1,1,3104)
