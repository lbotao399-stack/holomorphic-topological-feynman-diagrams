"""Coverage: which convergent weight orbits of a core admit a sequential plan (direct or after
completion + decompletion).   usage: python3 coverage.py Q1 [Q2 ...]"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import gutil
import seq
from validate import orbits

K33 = ((0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5))
PRISM = ((0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5))


def wheel(n):
    return tuple([(0, i) for i in range(1, n + 1)] + [(i, i % n + 1) for i in range(1, n + 1)])


def cores_list():
    out = [('K4', wheel(3)), ('W4', wheel(4)), ('prism', PRISM), ('K33', K33)]
    sk = json.load(open(os.path.join(HERE, '..', 'fiveloop', 'skel5.json')))
    for x in sk:
        out.append((x['name'], tuple(tuple(e) for e in x['edges'])))
    return out


def graph_sequential(E):
    wd = seq.merge([(u, v, 1) for u, v in E])
    return bool(seq.decompositions(wd))


if __name__ == '__main__':
    qs = [int(a) for a in sys.argv[1:]] or [3]
    for name, E in cores_list():
        line = '%-6s seq(G)=%-5s' % (name, graph_sequential(E))
        for q in qs:
            t = time.time()
            reps = orbits(E, q)
            nd = nc = 0
            for a in reps:
                c = seq.plans_for(list(E), list(a), q)
                if c:
                    if any(p['route'] == 'direct' for p in c):
                        nd += 1
                    else:
                        nc += 1
            line += ' | q=%d: %d orbits, direct %d, completion %d, none %d (%.0fs)' % (
                q, len(reps), nd, nc, len(reps) - nd - nc, time.time() - t)
        print(line, flush=True)
