"""Source basis and independent integer probes for the zero-displacement prism."""
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import prod
from pathlib import Path
from random import Random
import sys

import numpy as np
import sympy as sp

from periods import MAPS, PERMS, PRISM, support

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'k33_positive'))
from check import grid_input, numerator, permutation_sign
from vk import G, density_exact

SIGNS = tuple(map(permutation_sign, PERMS))
PAIRS = tuple(combinations(range(5),2))
GRAPH = G(6,PRISM)


def character():
    total = 0
    for g,sign in zip(MAPS,SIGNS):
        h = [1,-1,0,0,0,0]
        seen = set()
        for v in range(6):
            if v in seen:
                continue
            size = 0
            while v not in seen:
                seen.add(v)
                size += 1
                v = g[v]
            for k in range(size,6):
                h[k] += h[k-size]
        total += sign*(h[4]**2-h[3]*h[5])
    assert (total,len(MAPS)) == (480,12)
    return total//len(MAPS)


def brackets(probes):
    source = np.asarray(probes,dtype=np.int64)
    assert np.all(source.sum(axis=1)==0)
    p = source[:,[[g[v] for v in range(6)] for g in MAPS],:]
    return [p[:,:,a,0]*p[:,:,b,1] - p[:,:,a,1]*p[:,:,b,0] for a,b in PAIRS]


def column(tables, seed):
    return prod(tables[k] for k in seed) @ np.array(SIGNS)


class Pivots:
    def __init__(self):
        self.rows = {}

    def add(self, values):
        prime = 65521
        v = np.array(values,dtype=np.int64) % prime
        for k, row in sorted(self.rows.items()):
            v = (v-v[k]*row) % prime
        nonzero = np.flatnonzero(v)
        if not len(nonzero):
            return False
        k = int(nonzero[0])
        self.rows[k] = v*pow(int(v[k]), -1, prime) % prime
        return True


@cache
def sampling():
    dim = character()
    random = Random(748)
    entries = [(-1,-1,0,0,1,1),(-1,-1,-1,1,1,1),(-1,-1,0,0,0,2),
               (-2,-1,0,1,1,1),(-1,0,0,0,0,1)]
    probes = [[(1,0),(0,1),(1,1),(-2,-2),(0,0),(0,0)]]
    for _ in range(150):
        x = random.sample(random.choice(entries),6)
        y = random.sample(random.choice(entries),6)
        probes.append(list(zip(x,y)))
    tables = brackets(probes)
    seeds, cols, independent = [], [], Pivots()
    for seed in combinations_with_replacement(range(len(PAIRS)),4):
        # Standard two-row tableaux, in bracket-column order.
        if any(PAIRS[i][1]>PAIRS[j][1] for i,j in zip(seed,seed[1:])):
            continue
        col = column(tables,seed)
        if independent.add(col):
            seeds.append(seed)
            cols.append(col)
        if len(seeds)==dim:
            break
    assert len(seeds)==dim
    rows = np.array(cols).T
    independent = Pivots()
    indices = [i for i,row in enumerate(rows) if independent.add(row)]
    assert len(indices)==dim
    matrix = sp.Matrix(rows[indices].tolist())
    extra = next(i for i in range(len(rows)) if i not in indices and all(any(v) for v in probes[i])
                 and all(sp.Matrix([rows[i].tolist()])*matrix.inv()))
    return seeds, [probes[i] for i in indices+[extra]], matrix


@cache
def orbits():
    return {tuple(ex[e] for e in p):i for i,ex in enumerate(support()) for p in PERMS}


def weights(source, binary):
    trees, records = grid_input(source,PRISM)
    assert len(trees)==75
    uniform = 75**6*density_exact(GRAPH,[1]*9,source[:-1])
    poly = numerator(records,uniform,symmetrize=False,executable=binary)
    # Independent grounded Kirchhoff solution at unequal times.
    times = [2,1,1,3,1,2,1,1,2]
    u = sum(prod(t**k for k,t in zip(row,times)) for row in trees)
    power = [[t**k for k in range(4)] for t in times]
    expected = u**6*density_exact(GRAPH,times,source[:-1])
    value = 0
    result = [0]*len(support())
    for row in zip(*poly.nonzero()):
        ex = tuple(map(int,row))
        coefficient = int(poly[row])
        result[orbits()[ex]] += coefficient
        value += coefficient*prod(p[k] for p,k in zip(power,ex))
    assert value==expected
    return result
