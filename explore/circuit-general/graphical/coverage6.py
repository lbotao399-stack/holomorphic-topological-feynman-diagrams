"""Coverage of the 32 six-loop K5-cored graphs (D=2, n=6, census class I): series-reduced skeleton
(chains -> one edge, alpha adds), sequential plans directly or after completion, at levels q.
usage: python3 coverage6.py Q1 [Q2 ...]"""
import json
import os
import sys
import time
from collections import Counter

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, '..', 'cores'))
import gutil
import seq
from validate import orbits


def skeleton(edges):
    G = nx.MultiGraph()
    G.add_edges_from(edges)
    while True:
        biv = [v for v in G.nodes if G.degree(v) == 2 and len(set(G.neighbors(v))) == 2]
        if not biv:
            break
        v = biv[0]
        a, b = list(G.neighbors(v))
        G.remove_node(v)
        G.add_edge(a, b)
    G = nx.convert_node_labels_to_integers(G)
    return tuple(sorted((min(u, v), max(u, v)) for u, v in G.edges()))


if __name__ == '__main__':
    qs = [int(a) for a in sys.argv[1:]] or [3]
    data = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))['data']['2,6']['graphs']
    graphs = [r for r in data if r['k_gen'] is None]
    print('%d class-I graphs' % len(graphs))
    sk = {}
    for r in graphs:
        s = skeleton([tuple(e) for e in r['edges']])
        G = nx.Graph(list(s))
        key = None
        for k in sk:
            if nx.is_isomorphic(G, nx.Graph(list(k))):
                key = k
                break
        if key is None:
            key = s
            sk[key] = []
        sk[key].append(r['index'])
    print('%d distinct skeletons' % len(sk))
    for s, idx in sk.items():
        G = nx.Graph(list(s))
        multi = len(s) != G.number_of_edges()
        degs = sorted((d for _, d in G.degree()), reverse=True)
        seqG = bool(seq.decompositions(seq.merge([(u, v, 1) for u, v in s])))
        line = 'skeleton V=%d E=%d degs=%s multi=%s seq(G)=%s graphs=%s' % (
            G.number_of_nodes(), len(s), degs, multi, seqG, idx)
        print(line, flush=True)
        if multi or len(s) > int(os.environ.get('MAXE', '99')):
            continue
        for q in qs:
            t = time.time()
            reps = orbits(s, q)
            cnt = Counter()
            for a in reps:
                c = seq.plans_for(list(s), list(a), q)
                cnt['direct' if any(p['route'] == 'direct' for p in c) else ('completion' if c else 'none')] += 1
            print('    q=%d: %d orbits, %s (%.0fs)' % (q, len(reps), dict(cnt), time.time() - t), flush=True)
