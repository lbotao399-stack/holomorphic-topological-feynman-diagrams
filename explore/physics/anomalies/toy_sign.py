"""Relative signs between graphs: Wick sign of the 0-dimensional Grassmann model
(odd vertex measures mu_v, vertex mu_v f theta theta theta, contraction theta_{e0} theta_{e1} -> kappa p_{e0 e1})
  sign of  mu_0 B_0 mu_1 B_1 ...  ->  (pairs in edge order)(mu_0 ... mu_{m-1})(externals in ext order)
compared with the orientation-rule sign, random labellings, all bulk graphs n = 1..4 (m = 2n+2 vertices).
Claim (explore/anomaly-structure c_orientation_check (2), m = 6): toy = (-1)^{3m(m-1)/2} * rule."""
import json, random, sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/circuit-general/atlas')
import colour as C
from atlas import set_bulk3
from htk_free_relabel import relabel_any
rng = random.Random(23)

def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s

def toy_sign(nv, edges, ext):
    words, s_rule = C.rule_words(nv, edges, ext)
    toy = []
    for v in range(nv):
        toy += [('m', v)] + words[v]
    target = []
    for e in range(len(edges)):
        target += [('e', e, 0), ('e', e, 1)]
    target += [('m', v) for v in range(nv)] + [('x', v) for v in ext]
    pos = {h: i for i, h in enumerate(target)}
    return perm_sign([pos[h] for h in toy]), s_rule

import networkx as nx
for n in (1, 2, 3, 4):
    if n == 1:
        gl = [nx.cycle_graph(4)]
    else:
        gl = [it[0] for it in set_bulk3(3, n)]
    m = 2 * n + 2
    ok = []
    for G in gl:
        G = nx.convert_node_labels_to_integers(G)
        edges = list(G.edges())
        for _ in range(6):
            E, ext = relabel_any(G.number_of_nodes(), edges, rng)
            st, sr = toy_sign(G.number_of_nodes(), E, ext)
            ok.append(st == (-1) ** (3 * m * (m - 1) // 2) * sr)
    print('n=%d (m=%d): %d graphs x 6 random labellings: toy == (-1)^{3m(m-1)/2} rule: %s' % (n, m, len(gl), all(ok)), flush=True)
