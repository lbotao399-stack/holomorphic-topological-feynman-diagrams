"""Signed S2 periods on reproducible integer momenta at the cubic output vertices."""
from argparse import ArgumentParser
from fractions import Fraction
from itertools import cycle
import json
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'code'))
from htk import Graph, SPReduction
from htk.brute import jacobian_density
from htk.conventions import S2_kernel
from rooted_trees import graphs
from wick_weights import automorphisms, wick_sign


def probe(graph, outputs, number):
    basis = ((1,0,0),(0,1,0),(0,0,1))
    rotating = dict(zip(outputs, cycle(basis)))
    if graph.loops < 5:
        yield rotating
        return
    times = range(1, graph.ne+1)
    rng = Random(number)
    momenta = rotating
    for _ in range(24):
        if jacobian_density(graph, times, momenta):
            yield momenta
            break
        momenta = {v: rng.choice(basis) for v in outputs}
    else:
        yield {v: (1,k,k*k) for k,v in enumerate(outputs,1)}
    # A nonzero density at one point can still integrate to zero.
    rng = Random(400+number)
    for _ in range(24):
        momenta = {v: rng.choice(basis) for v in outputs}
        if jacobian_density(graph, times, momenta):
            yield momenta


def records(n):
    if n < 1:
        raise ValueError("rooted cubic periods require positive loop order")
    for number, edges in enumerate(graphs(n), 1):
        graph = Graph(2*n+2, edges)
        outputs = [v for v in range(graph.root) if graph.degree(v) == 2]
        reduction = SPReduction(graph)
        for momenta in probe(graph, outputs, number):
            value = S2_kernel(reduction.kernel(momenta, check=1), n)
            if value:
                break
        weight = Fraction(wick_sign(edges, graph.root), automorphisms(edges, graph.root))
        yield dict(number=number, plan='tight' if reduction.tight_plan else 'generic',
                   momenta=momenta, value=str(value), wick_weight=str(weight),
                   weighted_value=str(weight*value))


def periods(n):
    for record in records(n):
        yield record['number'], record['plan'] == 'tight', Fraction(record['value'])


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('loops', type=int, nargs='?', default=4)
    parser.add_argument('--jsonl', action='store_true', help='include every non-root output momentum')
    args = parser.parse_args()
    for record in records(args.loops):
        if args.jsonl:
            print(json.dumps(record, separators=(',', ':')), flush=True)
        else:
            print(record['number'], record['plan'], record['value'], flush=True)
