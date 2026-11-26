"""Scalar one-loop counterterms on six divergence-free vector waves, modulo 1009."""

from itertools import permutations, product
from math import factorial, prod
from os import environ
from random import Random

import numpy as np


PRIME = int(environ.get("TWIST_COUNTERTERM_PRIME", "1009"))
if PRIME not in (251, 1009):
    raise ValueError("verified replay primes are 251 and 1009")


def cycles(p):
    rest = set(range(len(p)))
    sizes = []
    while rest:
        i = next(iter(rest))
        size = 0
        while i in rest:
            rest.remove(i)
            size += 1
            i = p[i]
        sizes.append(size)
    return tuple(sorted(sizes))


def contraction_orbits():
    slots = tuple(range(6))
    perms = list(permutations(slots))
    character = 0
    for p in perms:
        sizes = cycles(p)
        character += prod(sum(k for k in sizes if r % k == 0) - 2 for r in sizes)
    assert character == 9 * factorial(6)
    maps = []
    orbit = []
    ids = {}
    seen = set()
    for f in product(slots, repeat=6):
        if f in seen or any(i == j for i, j in enumerate(f)):
            continue
        images = set()
        for p in perms:
            image = [0] * 6
            for i in slots:
                image[p[i]] = p[f[i]]
            images.add(tuple(image))
        k = len(ids)
        ids[f] = k
        maps.extend(images)
        orbit.extend([k] * len(images))
        seen.update(images)
    assert len(seen) == 5**6
    return np.asarray(maps).T, np.asarray(orbit), ids


def det(a):
    a = np.asarray(a, dtype=np.int64).copy() % PRIME
    value = 1
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j, i]), None)
        if pivot is None:
            return 0
        if pivot != i:
            a[[pivot, i]] = a[[i, pivot]]
            value = -value
        value = value * int(a[i, i]) % PRIME
        row = a[i] * pow(int(a[i, i]), -1, PRIME) % PRIME
        a[i + 1 :] = (a[i + 1 :] - a[i + 1 :, i, None] * row) % PRIME
    return value


def pivots(rows):
    a = np.array(rows, dtype=np.int64) % PRIME
    columns = []
    for i in range(a.shape[1]):
        k = len(columns)
        p = next((j for j in range(k, len(a)) if a[j, i]), None)
        if p is None:
            continue
        a[[p, k]] = a[[k, p]]
        a[k] = a[k] * pow(int(a[k, i]), -1, PRIME) % PRIME
        a[k + 1 :] = (a[k + 1 :] - a[k + 1 :, i, None] * a[k]) % PRIME
        columns.append(i)
    return columns


def waves(n, rng):
    p = np.array([[rng.randrange(1, PRIME) for _ in range(5)] for _ in range(n - 1)])
    p = np.vstack((p, -sum(p))) % PRIME
    u = np.array([[rng.randrange(1, PRIME) for _ in range(5)] for _ in range(n)])
    for q, v in zip(p, u):
        j = next(i for i, x in enumerate(q) if x)
        v[j] = (v[j] - (v @ q) * pow(int(q[j]), -1, PRIME)) % PRIME
    assert not np.any((p * u).sum(axis=1) % PRIME)
    return p, u


def evaluate(p, u, maps, orbit):
    pairings = u @ p.T % PRIME
    factors = np.ones(len(orbit), dtype=np.int64)
    for i in range(6):
        factors = factors * pairings[i, maps[i]] % PRIME
    # Every floating-point orbit sum is an integer at most 720*(PRIME-1).
    totals = np.bincount(orbit, weights=factors).astype(np.int64)
    return det(p[:5]) * totals % PRIME


def differential(p, u, maps, orbit):
    out = np.zeros(max(orbit) + 1, dtype=np.int64)
    for i in range(7):
        for j in range(i + 1, 7):
            rest = [k for k in range(7) if k not in (i, j)]
            q = np.vstack((p[i] + p[j], p[rest])) % PRIME
            bracket = (u[i] @ p[j]) * u[j] - (u[j] @ p[i]) * u[i]
            v = np.vstack((bracket, u[rest])) % PRIME
            out += (-1) ** (i + j) * evaluate(q, v, maps, orbit)
    return out % PRIME


def verify():
    maps, orbit, ids = contraction_orbits()
    rng = Random(211103049)
    scalar = [evaluate(*waves(6, rng), maps, orbit) for _ in range(16)]
    outgoing = [differential(*waves(7, rng), maps, orbit) for _ in range(16)]
    characteristic = {cycles(f): k for f, k in ids.items() if len(set(f)) == 6}
    columns = [characteristic[c] for c in ((2, 2, 2), (3, 3), (2, 4))]
    assert not np.any(np.array(outgoing)[:, columns])
    ranks = len(pivots(scalar)), len(pivots(outgoing)), len(pivots(np.array(scalar)[:, columns]))
    assert ranks == (9, 6, 3)
    print(PRIME, "scalar rank, outgoing rank, characteristic rank", *ranks)


if __name__ == "__main__":
    verify()
