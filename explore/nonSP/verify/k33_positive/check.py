"""Exact nonvanishing certificate for the four-source K3,3 kernel. Requires a
C++17 compiler (clang++ or g++), Python, NumPy, NetworkX, and SymPy."""
import json
import os
import subprocess
import sys
import tempfile
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
from graphs import K33
from vk import G, fdet

EDGES = range(9)
INJECTIONS = [(1, 0), (0, 1), (1, 1), (-2, -2), (0, 0), (0, 0)]
GRAPH = nx.Graph(K33)
EDGE_INDEX = {frozenset(pair): i for i, pair in enumerate(K33)}


def edge_permutation(automorphism):
    return tuple(EDGE_INDEX[frozenset((automorphism[a], automorphism[b]))] for a,b in K33)


def permutation_sign(perm):
    return (-1)**sum(p > q for i,p in enumerate(perm) for q in perm[i+1:])


AUTOMORPHISMS = [edge_permutation(m) for m in nx.isomorphism.GraphMatcher(
    GRAPH, GRAPH).isomorphisms_iter()]
assert len(AUTOMORPHISMS) == 72
for swap in ({0:2,2:0}, {1:3,3:1}, {0:1,1:0,2:3,3:2,4:5,5:4}):
    assert permutation_sign(edge_permutation({v:swap.get(v,v) for v in GRAPH})) == -1


def grid_input(injections=INJECTIONS, edges=K33):
    edge_index = {frozenset(pair): i for i, pair in enumerate(edges)}
    trees = []
    records = []
    for retained in combinations(EDGES, 5):
        tree = nx.Graph(edges[e] for e in retained)
        if len(tree) != 6 or not nx.is_tree(tree):
            continue
        complement = tuple(int(e not in retained) for e in EDGES)
        trees.append(complement)
        currents = np.zeros((9, 2), dtype=int)
        for v in range(5):
            path = nx.shortest_path(tree, v, 5)
            for a,b in zip(path, path[1:]):
                e = edge_index[frozenset((a,b))]
                currents[e] += (1 if edges[e] == (a,b) else -1)*np.array(injections[v])
        records.append([e for e in EDGES if complement[e]] + currents.flatten().tolist())
    if edges == K33:
        assert len(trees) == 81
    cycles = G(6, edges).C
    for first in combinations(range(1,9), 4):
        last = tuple(e for e in range(1,9) if e not in first)
        term = permutation_sign(first+last)
        for group in (first,last):
            term *= fdet([[row[e] for e in group] for row in cycles])
        if term:
            records.append([int(term),*first,*last])
    return trees, records


def numerator(records, uniform_value=2916, symmetrize=True, executable=None):
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        inp, binary, output = [directory / name for name in ('coefficients','grid','values')]
        inp.write_text(''.join(' '.join(map(str,row))+'\n' for row in records))
        if executable is None:
            subprocess.run([os.getenv('CXX','c++'),'-std=c++17','-O3',str(HERE/'grid.cpp'),
                            '-o',str(binary)],check=True)
            executable = binary
        tree_count = sum(len(row) == 22 for row in records)
        subprocess.run([str(executable),str(inp),str(output),str(tree_count)],check=True)
        values = np.fromfile(output,dtype=np.int64).astype(object).reshape([4]*9)
    inv6 = 6*sp.Matrix([[x**k for k in range(4)] for x in range(1,5)]).inv()
    transform = np.array(inv6.tolist(),dtype=np.int64).astype(object)
    for axis in EDGES:
        table = np.moveaxis(values,axis,0)
        out = (transform@table.reshape(4,-1)).reshape(table.shape)
        assert not np.any(out%6)
        values = np.moveaxis(out//6,0,axis)
    assert all(sum(index)==15 for index in zip(*values.nonzero()))
    assert values.sum()==uniform_value
    return sum(values.transpose(p) for p in AUTOMORPHISMS) if symmetrize else values


def margins_and_hepp_period(exponents, subsets, cycles, with_period):
    # The greedy spanning tree simultaneously minimizes the exponent on every
    # initial edge subset in a Hepp sector.
    partial = [Fraction(1)]
    for mask in range(1,512):
        edges = subsets[mask]
        margin = sum(exponents[j]+1 for j in edges) - 6*cycles[mask]
        if mask<511:
            assert margin>0
        if with_period:
            value = sum((partial[mask^(1<<j)] for j in edges),Fraction(0))
            partial.append(value/margin if mask<511 else value)
    if with_period:
        return partial[-1]


def main():
    trees,records = grid_input()
    symmetrized = numerator(records)
    first = {}
    sizes = {}
    def representative(a):
        if a not in first:
            orbit = {tuple(a[j] for j in p) for p in AUTOMORPHISMS}
            r = min(orbit)
            first.update(dict.fromkeys(orbit,r))
            sizes[r] = len(orbit)
        return first[a]
    numerators = Counter()
    for ex in zip(*symmetrized.nonzero()):
        key = representative(tuple(map(int,ex)))
        numerators[key] += int(symmetrized[ex])
    assert len(numerators)==175 and sum(sizes[a] for a in numerators)==10980

    cert = json.loads((HERE/'fields.json').read_text())
    denom = cert['denominator']
    # Every right hand side below is recorded in the basis of orbit *averages*.
    residual = Counter({a:denom*v for a,v in numerators.items()})
    tested = set(numerators)
    for h,edge,weight in cert['fields']:
        assert sum(h)==11 and min(h)>=0 and max(h)<=2
        for tree in trees:
            a = representative(tuple(i+j for i,j in zip(h,tree)))
            tested.add(a)
            residual[a] -= weight*(1+h[edge]-5*tree[edge])
    assert denom==10 and len(cert['fields'])==222
    assert len(+residual)==399 and min(residual.values())>=0
    assert min((+residual).values())==1103 and sum(residual.values())==1262466

    subsets=[tuple(j for j in EDGES if mask>>j&1) for mask in range(512)]
    loop_counts=[min(sum(tree[j] for j in edges) for tree in trees) for edges in subsets]
    bound = Fraction(0)
    for a in sorted(tested):
        hepp = margins_and_hepp_period(a,subsets,loop_counts,residual[a]>0)
        if hepp is not None:
            bound += residual[a]*hepp
    # Three denominators: group averaging, source polynomial = 6, rational fields.
    bound /= 72*6*denom
    assert bound == Fraction(6252279056460127,31603654656000)
    print('Cubic interpolation over all 4^9 integer time assignments passed.')
    print('The averaged numerator has 10980 monomials in 175 graph-automorphism orbits.')
    print('222 exact divergences give 399 positive orbit coefficients; minimum = 1103.')
    print('Strict lower bound for the scalar =',bound/81**6)
    print('Strict upper bound for the scalar =',bound)


if __name__=='__main__':
    main()
