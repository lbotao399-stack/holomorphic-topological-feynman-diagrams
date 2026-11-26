"""D=3, n=5: exact kernels of the series-parallel graphs of twistor BF Q_5 and bulk3(3,5) (enum_d3.json)
by htk.kernel_exact (series-parallel reduction, self-check against htk.brute at one point).
Sources: random integers in [-3,3]^3 at the source vertices (twistor: bivalent and O, root O;
bulk3: bivalent, root = last bivalent), seed 1000*idx + 7 (atlas/e3 convention).
usage: python3 run_d3_sp.py SET (twistor | bulk3)"""
import json
import os
import random
import sys
import time

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'code')))
import htk                                    # noqa: E402
from htk.graph import from_networkx           # noqa: E402

import signal

d = json.load(open(os.path.join(HERE, 'enum_d3.json')))
which = sys.argv[1]
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 0          # seconds per graph (0: none)
SKIP = set(int(x) for x in sys.argv[3].split(',')) if len(sys.argv) > 3 and sys.argv[3] else set()


def _alarm(sig, frame):
    raise TimeoutError('time limit')


signal.signal(signal.SIGALRM, _alarm)
items = []
if which == 'twistor':
    for key in ('tree', 'unicyclic'):
        for r in d['twistor'][key]:
            items.append((key, r, r['O']))
else:
    for r in d['bulk3']:
        items.append(('bulk3', r, r['S0'][-1]))
for idx, (key, r, root) in enumerate(items):
    if r['cls'] != 'SP' or idx in SKIP:
        continue
    G = nx.Graph([tuple(e) for e in r['edges']])
    Gh, lab = from_networkx(G, root=root)
    rng = random.Random(1000 * idx + 7)
    lam = {lab[v]: [rng.randint(-3, 3) for _ in range(3)] for v in r['S0'] if v != root}
    full = [lam.get(v, [0, 0, 0]) for v in range(Gh.nv - 1)]
    t = time.time()
    try:
        signal.alarm(LIMIT)
        K = htk.kernel_exact(Gh, full, check=1)
        signal.alarm(0)
        msg = str(K)
    except Exception as e:
        signal.alarm(0)
        msg = 'FAILED %r' % (e,)
    print(json.dumps({'set': key, 'idx': idx, 'nv': Gh.nv, 'edges': [list(e) for e in Gh.edges], 'sources': full,
                      'K': msg, 'seconds': round(time.time() - t, 1)}), flush=True)
