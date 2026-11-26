"""D = 3, n = 5 graph sets: bulk3(3,5) (max degree 3, 4 legs) and twistor BF Q_5 (Gamma - O connected,
cycle rank 0 or 1, others degree <= 3).  Skeletons by degree sequence (bulk3: cubic on 8 vertices;
twistor: deg O = 6 (tree) with 4 cubic, deg O = 5 (unicyclic) with 5 cubic), subdivisions with
|I| = 16 checked by graphs_enum.local_ok (admissibility f(S) >= 0), isomorphism classes; for each
graph: series-parallel core (census_cores.core_name), classes SP / B / R_k (atlas reduce_search)."""
import itertools
import json
import os
import sys
import time
from collections import Counter

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.join(HERE, '..', 'atlas')
sys.path.insert(0, ATLAS)
import graphs_enum                                   # noqa: E402
from reduce_search import Network, search            # noqa: E402
from htk import catalog                              # noqa: E402

import census_cores                                  # noqa: E402

D, n = 3, 5
N = D * n + 1


def multigraphs(degs):
    """2-connected loopless multigraphs with the given degree sequence (vertices 0..s-1), up to iso"""
    s = len(degs)
    out = []
    found = []

    def rec(rem, edges):
        v = next((i for i in range(s) if rem[i] > 0), None)
        if v is None:
            G = nx.MultiGraph()
            G.add_nodes_from(range(s))
            G.add_edges_from(edges)
            Gs = nx.Graph(G)
            if not nx.is_connected(Gs) or (s >= 3 and not nx.is_biconnected(Gs)):
                return
            for H in found:
                if nx.is_isomorphic(G, H):
                    return
            found.append(G)
            out.append(list(edges))
            return
        for w in range(v + 1, s):
            if rem[w] > 0 and (not edges or (v, w) >= edges[-1]):
                rem[v] -= 1
                rem[w] -= 1
                rec(rem, edges + [(v, w)])
                rem[v] += 1
                rem[w] += 1
    rec(list(degs), [])
    return out


def graphs_from(skels, keep=None):
    buckets, res = {}, []
    for sk in skels:
        s = max(max(e) for e in sk) + 1
        for lens in graphs_enum.compositions(N, len(sk)):
            seen1 = set()
            bad = False
            for (a, b), l in zip(sk, lens):
                if l == 1:
                    k = (min(a, b), max(a, b))
                    if k in seen1:
                        bad = True
                        break
                    seen1.add(k)
            if bad or not graphs_enum.local_ok(D, s, sk, lens):
                continue
            nv, E = graphs_enum.subdivide(s, sk, lens)
            G = nx.Graph()
            G.add_nodes_from(range(nv))
            G.add_edges_from(E)
            if G.number_of_edges() != N:
                continue
            h = (tuple(sorted(d for _, d in G.degree())), nx.weisfeiler_lehman_graph_hash(G, iterations=4))
            lst = buckets.setdefault(h, [])
            if any(nx.is_isomorphic(G, H) for H in lst):
                continue
            lst.append(G)
            res.append(G)
    return res


def rooted_orbits(G, roots):
    autos = list(nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
    seen, reps = set(), []
    for O in roots:
        if O in seen:
            continue
        seen |= {g[O] for g in autos}
        reps.append(O)
    return reps


def classify(G, S0):
    net = Network(G, S0)
    rp = search(net, ("S", "Y0"))
    if rp["reducible"]:
        return 0
    rr = search(net, ("S", "Y0", "Y1"))
    return rr["min_Y1"]


if __name__ == '__main__':
    t0 = time.time()
    out = {}
    # bulk3: cubic skeletons on 8 vertices
    cub = multigraphs([3] * 8)
    print('cubic loopless 2-connected multigraph skeletons, 8 vertices:', len(cub), flush=True)
    bulk = graphs_from(cub)
    print('bulk3(3,5): %d graphs  %.1fs' % (len(bulk), time.time() - t0), flush=True)
    rows = []
    cnt = Counter()
    for i, G in enumerate(bulk):
        G = nx.convert_node_labels_to_integers(G)
        S0 = sorted(v for v in G if G.degree(v) == 2)
        core, _ = census_cores.core_name([tuple(e) for e in G.edges()])
        sp = bool(catalog.is_sp(G))
        k = classify(G, S0)
        cls = 'SP' if sp else ('B' if k == 0 else ('I' if k is None else 'R%d' % k))
        cnt[(core, cls)] += 1
        rows.append(dict(index=i, nv=G.number_of_nodes(), edges=sorted(tuple(sorted(e)) for e in G.edges()),
                         core=core, cls=cls, S0=S0))
    print('bulk3(3,5) core x class:', flush=True)
    for (c, k), m in sorted(cnt.items()):
        print('   %-8s %-3s %d' % (c, k, m))
    out['bulk3'] = rows
    # twistor
    tw = {}
    for key, degs in (('tree', [6, 3, 3, 3, 3]), ('unicyclic', [5, 3, 3, 3, 3, 3])):
        sks = multigraphs(degs)
        gs = graphs_from(sks)
        items = []
        cnt = Counter()
        for G in gs:
            G = nx.convert_node_labels_to_integers(G)
            cand = []
            for O in G:
                Hh = G.copy()
                Hh.remove_node(O)
                if not nx.is_connected(Hh) or any(G.degree(v) > 3 for v in Hh):
                    continue
                r = Hh.number_of_edges() - Hh.number_of_nodes() + 1
                if r == (0 if key == 'tree' else 1):
                    cand.append(O)
            for O in rooted_orbits(G, cand):
                S0 = sorted(set(v for v in G if G.degree(v) == 2) | {O})
                core, _ = census_cores.core_name([tuple(e) for e in G.edges()])
                sp = bool(catalog.is_sp(G))
                k = classify(G, S0)
                cls = 'SP' if sp else ('B' if k == 0 else ('I' if k is None else 'R%d' % k))
                cnt[(core, cls)] += 1
                items.append(dict(nv=G.number_of_nodes(), edges=sorted(tuple(sorted(e)) for e in G.edges()),
                                  O=O, core=core, cls=cls, S0=S0))
        print('twistor %s: %d skeletons, %d graphs, %d rooted  %.1fs' % (key, len(sks), len(gs), len(items),
                                                                        time.time() - t0), flush=True)
        for (c, k), m in sorted(cnt.items()):
            print('   %-8s %-3s %d' % (c, k, m))
        tw[key] = items
    out['twistor'] = tw
    json.dump(out, open(os.path.join(HERE, 'enum_d3.json'), 'w'))
    print('%.1fs' % (time.time() - t0))
