"""Which simple five-loop skeletons (skel5.json) are skeletons of admissible graphs at D = 2, 3, 4
(n = 5, |I| = D n + 1): count of subdivisions (length compositions, up to skeleton automorphisms not
quotiented) passing graphs_enum.local_ok, and of isomorphism classes."""
import json
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'atlas'))
import graphs_enum      # noqa: E402

sk = json.load(open(os.path.join(HERE, 'skel5.json')))
for D in (2, 3, 4):
    N = D * 5 + 1
    row = []
    for r in sk:
        E = [tuple(e) for e in r['edges']]
        s = r['V']
        classes = []
        nlen = 0
        for lens in graphs_enum.compositions(N, len(E)):
            if not graphs_enum.local_ok(D, s, E, lens):
                continue
            nlen += 1
            nv, Es = graphs_enum.subdivide(s, E, lens)
            G = nx.Graph(Es)
            if not any(nx.is_isomorphic(G, H) for H in classes):
                classes.append(G)
        row.append('%s %d' % (r['name'], len(classes)))
    print('D=%d n=5 (|I|=%d): graphs with simple skeleton =' % (D, N), ', '.join(row), flush=True)
