"""Float cross-check of generic K5 values (tables/K5_q<q>.json, kind G) by tropical sampling
(fiveloop/tropmc.py, Hepp sectors, Borinsky arXiv:2008.12310).
usage: python3 check_tropmc.py Q NVEC NSAMPLES [seed]"""
import json
from fractions import Fraction as F
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'cores'))
sys.path.insert(0, os.path.join(HERE, '..', 'fiveloop'))
from tropmc import Trop  # noqa: E402

Z = {'1': 1.0, 'z3': 1.2020569031595942854, 'z5': 1.0369277551433699263, 'z7': 1.0083492773819228268}


def fval(v):
    tot = 0.0
    for k, x in v.items():
        m = 1.0
        for f in k.split('*'):
            m *= Z[f]
        tot += float(F(x)) * m
    return tot


if __name__ == '__main__':
    q, nvec, ns = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    d = json.load(open(os.path.join(HERE, 'tables', 'K5_q%d.json' % q)))
    E = [tuple(e) for e in d['edges']]
    gen = [r for r in d['values'] if r[1] == 'G']
    random.Random(seed).shuffle(gen)
    for a, kd, v in gen[:nvec]:
        ex = fval(v)
        m, e = Trop(E, list(a), q).estimate(ns, seed=seed)
        print('K5 q=%d %s [G]: exact %.6e  tropical MC %.6e +- %.2e  (%.1f sigma)' % (
            q, tuple(a), ex, m, e, abs(m - ex) / e), flush=True)
