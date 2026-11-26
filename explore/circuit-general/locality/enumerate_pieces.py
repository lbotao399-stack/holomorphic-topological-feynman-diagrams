"""Two-terminal decompositions Gamma = S cup_{p,q} O of the D = 3, n = 4 graph sets of the atlas
(bulk3(3,4): sources at bivalent vertices; twistor4 non-SP: sources at bivalent vertices and O).
For each 2-vertex cut {p,q} and each union of components C of Gamma - {p,q} (S = edges meeting C):
  N_S, n_S, delta_S = N_S - 1 - 3 n_S (tight iff 0), S series-parallel (S + pq), number of internal
  sources, Gamma' = Gamma - C + e* series-parallel.  Output: per graph the decompositions with
  S non-SP or Gamma' non-SP, delta_S = 0 (tight S) or delta_O = 0."""
import json
import os
import sys
from itertools import combinations

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '../atlas'))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '../../../code')))
from atlas import set_bulk3, set_twistor          # noqa: E402
from htk import catalog                           # noqa: E402

H = 3


def pieces(G, S0):
    out = []
    V = list(G.nodes())
    for p, q in combinations(V, 2):
        R = G.copy()
        R.remove_nodes_from([p, q])
        comps = [set(c) for c in nx.connected_components(R)]
        if len(comps) < 2:
            continue
        for r in range(1, len(comps)):
            for sel in combinations(range(len(comps)), r):
                C = set().union(*[comps[i] for i in sel])
                Sed = [(a, b, k) for a, b, k in G.edges(keys=True) if a in C or b in C] if G.is_multigraph() \
                    else [(a, b) for a, b in G.edges() if a in C or b in C]
                NS = len(Sed)
                VS = len(C) + 2
                nS = NS - VS + 1
                Sg = nx.MultiGraph()
                Sg.add_edges_from([e[:2] for e in Sed])
                Sg.add_edge(p, q)
                sp_S = bool(catalog.is_sp(nx.Graph(Sg)) if not Sg.is_multigraph() else catalog.is_sp(Sg))
                Og = nx.MultiGraph()
                Og.add_edges_from([e[:2] for e in G.edges() if not (e[0] in C or e[1] in C)])
                Og.add_edge(p, q)
                sp_O = bool(catalog.is_sp(Og))
                srcS = len([v for v in C if v in S0])
                deltaS = NS - 1 - H * nS
                deltaO = (G.number_of_edges() - NS) - 1 - H * (G.number_of_edges() - NS - (G.number_of_nodes() - len(C)) + 1 - 1)
                out.append(dict(p=p, q=q, C=sorted(C), NS=NS, nS=nS, dS=deltaS, dO=H - 1 - deltaS,
                                spS=sp_S, spO=sp_O, srcS=srcS, pq_src=[p in S0, q in S0]))
    return out


def main():
    rows = []
    b3 = set_bulk3(3, 4)
    for idx, (G, O, S0) in enumerate(b3):
        S0 = [v for v in G if G.degree(v) == 2]
        rows.append(('bulk3', idx, G, S0[-1], S0))
    tw = set_twistor(4)
    items = [it for k in ('tree', 'unicyclic') for it in tw[k]]
    for idx, (G, O, S0) in enumerate(items):
        rows.append(('twistor4', idx, G, O, S0))
    summary = []
    for setname, idx, G, root, S0 in rows:
        if catalog.is_sp(G):
            continue
        pcs = pieces(G, S0)
        tightS_nonSP = [x for x in pcs if x['dS'] == 0 and not x['spS']]
        tightS_SP_coreO = [x for x in pcs if x['dS'] == 0 and x['spS'] and not x['spO']]
        summary.append(dict(set=setname, idx=idx, edges=sorted(map(list, G.edges())), root=root, S0=S0,
                            n_pieces=len(pcs), tightS_nonSP=tightS_nonSP, tightS_SP_coreO=tightS_SP_coreO))
        print('%-8s #%-3d cuts=%-3d tight nonSP S: %s | tight SP S with nonSP Gamma\': %s' % (
            setname, idx, len(pcs),
            [(x['p'], x['q'], x['NS'], x['nS'], x['srcS'], 'O-SP' if x['spO'] else 'O-core') for x in tightS_nonSP],
            [(x['p'], x['q'], x['NS'], x['srcS']) for x in tightS_SP_coreO]), flush=True)
    json.dump(summary, open(os.path.join(HERE, 'pieces_d3n4.json'), 'w'), default=str)
    nb = sum(1 for s in summary if s['tightS_nonSP'])
    nc = sum(1 for s in summary if s['tightS_SP_coreO'])
    print('non-SP graphs: %d; with a tight non-SP two-terminal S: %d; with a tight SP S and non-SP Gamma\': %d'
          % (len(summary), nb, nc))


if __name__ == '__main__':
    main()
