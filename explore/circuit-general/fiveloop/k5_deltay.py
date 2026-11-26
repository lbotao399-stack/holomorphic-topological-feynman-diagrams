"""D=2, n=6: the 32 census graphs without circuit chart (I-core K5).  Graph level: Delta-Y on every
triangle of Gamma (three unsubdivided skeleton edges), new vertex s; series-parallel core of the result.
Period level: (Y) relation I_q(K5; c) = I_q(K5 - triangle + star; a) at a unique triangle (sum c = 2q),
K5 - triangle + star = K33+e; K5 periods at q = 3 from the K33+e q=3 table; consistency over all
unique triangles and against (C) completion relations."""
import json
import os
from collections import Counter
from fractions import Fraction as F

import networkx as nx

import census_cores
import per5
import gutil
import rel5
import skel5

HERE = os.path.dirname(os.path.abspath(__file__))
cen = json.load(open(os.path.join(HERE, '..', 'atlas', 'census.json')))
gs = cen['data']['2,6']['graphs']
irr = [g for g in gs if not g['sp'] and not g['B'] and g['k_gen'] is None]
print('D=2 n=6 graphs without reduction:', len(irr))
K5 = nx.complete_graph(5)
stat = Counter()
for g in irr:
    E = [tuple(e) for e in g['edges']]
    core, sk = census_cores.core_name(E)
    G = nx.Graph(E)
    tris = [t for t in nx.enumerate_all_cliques(G) if len(t) == 3]
    res = Counter()
    for a, b, c in tris:
        H = G.copy()
        H.remove_edges_from([(a, b), (a, c), (b, c)])
        s = max(G.nodes()) + 1
        H.add_edges_from([(s, a), (s, b), (s, c)])
        cn, _ = census_cores.core_name(list(H.edges()))
        res[cn] += 1
    stat[(core, len(tris), tuple(sorted(res.items())))] += 1
for k, v in sorted(stat.items()):
    print('  SP core %s, triangles %d, Delta-Y cores %s : %d graphs' % (k[0], k[1], dict(k[2]), v))

# period level, q = 3
q = 3
E = sorted(tuple(sorted(e)) for e in K5.edges())
reps = set()
autos = gutil.edge_automorphisms(tuple(E))
for c in gutil.compositions(6 * q, 10, 1, q - 1):
    if gutil.convergent(tuple(E), list(c), q):
        reps.add(gutil.orbit_rep(c, autos))
print('K5 q=3: %d convergent orbits' % len(reps))
for r in sorted(reps):
    r = list(r)
    vals = []
    for tri in rel5.triangles(E):
        if r[tri[3]] + r[tri[4]] + r[tri[5]] == 2 * q:
            ne, npw = rel5.tri_to_star(E, r, tri, q)
            vals.append(per5.value(ne, npw, q))
    vals = [v for v in vals if v is not None]
    same = all(v == vals[0] for v in vals)
    print('  I_3(K5; %s) = %s   [%d unique triangles, all equal: %s]' % (r, per5.fmt(vals[0]) if vals else '?', len(vals), same))
    # completion relations
    nc = bad = 0
    for items in rel5.completion_items(E, r, q):
        tgt = per5.value(items[1][0], items[1][1], q)
        if tgt is None or not vals:
            continue
        nc += 1
        bad += (per5.vadd(vals[0], tgt, items[1][2]) != {})
    print('     (C) relations to evaluable graphs: %d checked, %d bad' % (nc, bad))
