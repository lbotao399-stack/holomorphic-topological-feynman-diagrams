"""Algebra layer: Jacobi, invariance, vertex symmetry, hatQ0^2 = 0 on random jet monomials, for every
theory family and both index conventions of G.  Log: test_alg.log"""
import random
import sys

import numpy as np

import superalg as sa
import polyops as po
from wick import Algebra
import theories


def rand_mono(alg, rng, nlet, maxjet):
    lets = []
    while len(lets) < nlet:
        A = rng.randrange(alg.N)
        al = tuple(rng.randint(0, maxjet) for _ in range(alg.H))
        lets.append((A, al))
    s, m = po.canon(lets, alg.pi)
    return m if s else None


def main():
    rng = random.Random(1)
    out = []
    for name, th in theories.catalogue_small().items():
        for variant in (0, 1):
            alg = Algebra(th['alg'], th['H'], variant=variant)
            Algebra._q0cache = {}
            bad = 0
            tested = 0
            while tested < 25:
                m = rand_mono(alg, rng, rng.randint(1, 3), 1)
                if m is None:
                    continue
                tested += 1
                r = alg.q0(alg.q0({m: 1}))
                if r:
                    bad += 1
            line = '%-28s H=%d N=%3d variant %d: Q0^2 = 0 on %d/%d random monomials' % (
                name, th['H'], alg.N, variant, tested - bad, tested)
            print(line, flush=True)
            out.append(line)


if __name__ == '__main__':
    main()
