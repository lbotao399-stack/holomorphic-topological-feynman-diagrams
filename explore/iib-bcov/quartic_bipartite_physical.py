"""Exact physical sums from certified bipartite orbit periods in matrix-edge order."""
from argparse import ArgumentParser
from gzip import open as open_gzip
from itertools import islice, permutations, product
from pathlib import Path

import numpy as np
from flint import fmpq


def matrix_orderings():
    graph_edges = [(0,1), (0,3), (0,5), (1,2), (1,4), (2,3), (2,5), (3,4), (4,5)]
    orderings = []
    for left, right in (((0,2,4), (1,3,5)), ((1,3,5), (0,2,4))):
        for rows, columns in product(permutations(left), permutations(right)):
            orderings.append([graph_edges.index(tuple(sorted((i,j)))) for i in rows for j in columns])
    return orderings


def physical_weights(bank):
    weights = {}
    shifts = np.array([16**j for j in range(8,-1,-1)], dtype=np.int64)
    orderings = matrix_orderings()
    with open_gzip(bank, "rt") as stream:
        while batch := list(islice(stream, 8192)):
            data = [v.split() for v in batch]
            exponents = np.array([v[:9] for v in data], dtype=np.int64)
            keys = np.full(len(data), 16**9, dtype=np.int64)
            for labels in orderings:
                np.minimum(keys, exponents[:, labels] @ shifts, out=keys)
            for key, row in zip(keys.tolist(), data, strict=True):
                if key in weights:
                    raise ValueError("Duplicate core orbit in the physical bank")
                weights[key] = (fmpq(row[9]), fmpq(row[10]))
    return weights


def sum_periods(directory, weights, homogeneous):
    sums = [[fmpq(0) for _ in range(3)] for _ in range(2)]
    annihilator = [fmpq(0), fmpq(0)]
    overlap = 0
    files = [open(directory / f"channel_{i}.tsv") for i in range(3)]
    try:
        for rows in zip(*files, strict=True):
            entries = [r.split() for r in rows]
            if len({a for a, b in entries}) != 1:
                raise ValueError("Component tables list different orbit keys")
            key = int(entries[0][0])
            if key not in weights:
                continue
            coefficients = weights.pop(key)
            for column, (_, value) in enumerate(entries):
                period = fmpq(value)
                for j in range(2):
                    sums[j][column] += coefficients[j] * period
            if key in homogeneous:
                overlap += 1
                for j in range(2):
                    annihilator[j] += coefficients[j] * homogeneous[key]
    finally:
        for f in files:
            f.close()
    if weights:
        raise ValueError("Some physical orbits are missing")
    return sums, annihilator, overlap


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--homogeneous", type=Path)
    args = parser.parse_args()
    bank = Path(__file__).with_name("quartic_k33.weights.gz")
    homogeneous = {}
    if args.homogeneous:
        with open(args.homogeneous) as f:
            homogeneous = {int(key):fmpq(value) for key,value in map(str.split,f)}
    sums, kernel, overlap = sum_periods(args.directory, physical_weights(bank), homogeneous)
    for label, values in zip(("K","J"), sums):
        print(label, *(str(x) for x in values))
    if homogeneous:
        assert kernel == [0,0]
        print("Homogeneous support", len(homogeneous), "physical overlap", overlap, "pairings", *kernel)
