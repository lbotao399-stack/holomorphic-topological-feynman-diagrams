"""SP skeleton (series-parallel core) of every census graph, (D, n) = (2, 5) and others."""
import json
import os
import sys
from collections import Counter

import networkx as nx

import skel5

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.join(HERE, '..', 'atlas')
cen = json.load(open(os.path.join(ATLAS, 'census.json')))
sk = json.load(open(os.path.join(HERE, 'skel5.json')))
lib = [(r['name'], nx.Graph([tuple(e) for e in r['edges']])) for r in sk]
LOW = {'K4': nx.complete_graph(4), 'W4': nx.wheel_graph(5),
       'prism': nx.Graph([(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (0, 3), (1, 4), (2, 5)]),
       'K33': nx.complete_bipartite_graph(3, 3)}
for k, G in LOW.items():
    lib.append((k, G))


def core_name(edges):
    E = skel5.sp_skeleton(edges)
    if len(E) <= 1:
        return 'SP', E
    G = nx.Graph(E)
    for nm, H in lib:
        if G.number_of_nodes() == H.number_of_nodes() and G.number_of_edges() == H.number_of_edges() \
                and nx.is_isomorphic(G, H):
            return nm, E
    return 'other(V=%d,E=%d)' % (G.number_of_nodes(), G.number_of_edges()), E


if __name__ == '__main__':
    for key in sys.argv[1:] or ['2,5']:
        gl = cen['data'][key]['graphs']
        cnt = Counter()
        bycls = {}
        for g in gl:
            nm, E = core_name(g['edges'])
            cls = 'SP' if g['sp'] else ('B' if g['B'] else ('I' if g['k_gen'] is None else 'R%d' % g['k_gen']))
            cnt[(nm, cls)] += 1
            bycls.setdefault(nm, []).append(g['index'])
        print('(D,n) =', key, len(gl), 'graphs')
        for (nm, cls), c in sorted(cnt.items()):
            print('  core %-8s class %-3s : %d' % (nm, cls, c))
        for nm, idx in sorted(bycls.items()):
            print('  %-8s indices %s' % (nm, idx))
