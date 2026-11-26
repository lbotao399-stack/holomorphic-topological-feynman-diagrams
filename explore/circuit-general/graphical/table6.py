"""Six-loop K5-cored skeletons (32 census class-I graphs, D=2, n=6): all convergent orbits at level q
that admit a sequential plan (always via completion here); writes tables/K5core_<first census index>_q<q>.json.
A second, structurally different plan is run for every orbit as an internal cross-check.
usage: python3 table6.py Q"""
import json
import os
import sys
import time
from fractions import Fraction as F

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import seq
import mzv
from validate import orbits, to_labels
from coverage6 import skeleton

if __name__ == '__main__':
    q = int(sys.argv[1])
    data = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))['data']['2,6']['graphs']
    graphs = [r for r in data if r['k_gen'] is None]
    sk = {}
    for r in graphs:
        s = skeleton([tuple(e) for e in r['edges']])
        G = nx.Graph(list(s))
        key = next((k for k in sk if nx.is_isomorphic(G, nx.Graph(list(k)))), None)
        if key is None:
            key = s
            sk[key] = []
        sk[key].append(r['index'])
    tot = agree2 = 0
    for s, idx in sk.items():
        if len(s) > int(os.environ.get('MAXE', '99')):
            continue
        reps = orbits(s, q)
        out = {}
        t0 = time.time()
        for a in reps:
            cands = seq.plans_for(list(s), list(a), q)
            if not cands:
                continue
            pref = F(1, seq.gamma_int(q))
            for x in a:
                pref *= seq.gamma_int(q - x)
            t = time.time()
            try:
                v = mzv.scal(seq.run(cands[0]), pref)
            except Exception as e:
                print('  [%d] %s  FAILED %s: %s' % (idx[0], a, type(e).__name__, e), flush=True)
                continue
            sig0 = (cands[0]['first'], tuple((o[0], str(o[1])) for o in cands[0]['ops']))
            other = next((p for p in cands[1:] if (p['first'], tuple((o[0], str(o[1])) for o in p['ops'])) != sig0
                          and p['route'] != cands[0]['route']), None)
            chk = ''
            if other is not None:
                try:
                    v2 = mzv.scal(seq.run(other), pref)
                except Exception as e:
                    v2 = None
                    chk = ' second plan failed: %s' % type(e).__name__
                if v2 is not None:
                    chk = ' second plan (%s): %s' % (other['route'], 'equal' if v2 == v else 'DIFFERENT %s' % mzv.tostr(v2))
                    agree2 += (v2 == v)
            tot += 1
            out[a] = to_labels(v)
            print('  [%d] %s  %s  (%.1fs, %s, %d appends)%s' % (idx[0], a, mzv.tostr(v), time.time() - t,
                                                              cands[0]['route'], cands[0]['appends'], chk), flush=True)
        os.makedirs(os.path.join(HERE, 'tables'), exist_ok=True)
        json.dump({'core': 'K5core_%d' % idx[0], 'census_indices': idx, 'edges': [list(e) for e in s], 'q': q,
                   'orbits_total': len(reps), 'seconds': round(time.time() - t0, 1),
                   'values': [[list(r), {k: str(x) for k, x in val.items()}] for r, val in sorted(out.items())]},
                  open(os.path.join(HERE, 'tables', 'K5core_%d_q%d.json' % (idx[0], q)), 'w'))
        print('skeleton %s (census %s): %d/%d orbits (%.0fs)' % (s, idx, len(out), len(reps), time.time() - t0),
              flush=True)
    print('total %d orbits; second-plan agreements %d' % (tot, agree2))
