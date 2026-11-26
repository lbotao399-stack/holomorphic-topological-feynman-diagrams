"""Claim a4 (flows) and the counts used in b1 (40 flows per graph at m=1, 12288 letter maps, |Aut|).

Flow on an edge: from its mu end to its gamma end (all edges mu-gamma, claim a2).  At vertex v:
in = #gamma letters contracted, out = #mu letters contracted.  Letters: b_0 = (ga, mu, mu), c = (ga, ga, ga).
"""
import sys, os
from collections import Counter
from itertools import product
from math import comb
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcommon import theta_graph, PQR

LET = {'b': {'ga': 1, 'mu': 2}, 'c': {'ga': 3, 'mu': 0}}


def falling(n, k):
    r = 1
    for i in range(k):
        r *= n - i
    return r


def classify(nv, edges, o):
    deg = Counter()
    ind = Counter(); outd = Counter()
    for e, (a, b) in enumerate(edges):
        t, h = (a, b) if o[e] == 0 else (b, a)
        outd[t] += 1; ind[h] += 1; deg[a] += 1; deg[b] += 1
    types, nmaps, ext = {}, 1, Counter()
    for v in range(nv):
        i, u = ind[v], outd[v]
        cand = [T for T in ('b', 'c') if LET[T]['ga'] >= i and LET[T]['mu'] >= u]
        # external letters left: interior vertices need >= 1 (Dtheta), at K = 0
        ok = []
        for T in cand:
            left = LET[T]['ga'] - i + LET[T]['mu'] - u
            if deg[v] == 2 and left == 0:
                continue
            ok.append(T)
        if len(ok) != 1:
            return None
        T = ok[0]
        types[v] = T
        nmaps *= falling(LET[T]['ga'], i) * falling(LET[T]['mu'], u)
        ext['mu'] += LET[T]['mu'] - u
        ext['ga'] += LET[T]['ga'] - i
    return types, nmaps, ext, ind, outd


if __name__ == '__main__':
    for pqr in PQR:
        nv, edges, paths = theta_graph(*pqr)
        byM = Counter(); maps_m1 = Counter(); roles = Counter(); struct_ok = True
        nonunique = 0
        for o in product((0, 1), repeat=len(edges)):
            r = classify(nv, edges, o)
            if r is None:
                continue
            types, nmaps, ext, ind, outd = r
            m = sum(1 for v in types if types[v] == 'c')
            byM[(m, ext['mu'], ext['ga'])] += 1
            if m == 1:
                maps_m1[nmaps] += 1
                cv = [v for v in types if types[v] == 'c'][0]
                if cv == 0:
                    roles['c@N'] += 1
                elif cv == nv - 1:
                    roles["c@N'"] += 1
                else:
                    roles['c@interior'] += 1
                # structural statements of a4 at m = 1
                srcs = [v for v in range(1, nv - 1) if (ind[v], outd[v]) == (0, 2)]
                sinks = [v for v in range(1, nv - 1) if (ind[v], outd[v]) == (2, 0)]
                if cv in (0, nv - 1):
                    struct_ok &= (len(srcs) == 1 and len(sinks) == 0)
                else:
                    struct_ok &= (len(srcs) == 0 and sinks == [cv])
                    # path through the sink directed into it; the other two paths: one N->N', one N'->N
                    dirs = []
                    for pv, pe in paths:
                        if cv in pv:
                            k = pv.index(cv)
                            # edges pe[0..k] towards N' (o=0), pe[k+1..] towards N (o=1)
                            struct_ok &= all(o[e] == 0 for e in pe[:k + 1]) and all(o[e] == 1 for e in pe[k + 1:])
                        else:
                            s = set(o[e] for e in pe)
                            struct_ok &= (len(s) == 1)
                            dirs.append(s.pop())
                    struct_ok &= (sorted(dirs) == [0, 1])
        tot = sum(byM.values())
        pred = {(m, 11 - 2 * m, 2 * m - 1): 4 * comb(10, 2 * m - 1) for m in range(1, 6)}
        print(pqr, 'admissible flows', tot, '; by (m, ext mu, ext ga):', dict(sorted(byM.items())),
              '; equals 4*binom(10,2m-1):', dict(byM) == pred)
        print('    m=1: roles', dict(roles), '; letter maps per flow', dict(maps_m1), '; a4 structure at m=1:', struct_ok)
    # automorphism groups
    for pqr in PQR:
        nv, edges, paths = theta_graph(*pqr)
        G = nx.Graph(); G.add_nodes_from(range(nv)); G.add_edges_from(edges)
        naut = sum(1 for _ in GraphMatcher(G, G).isomorphisms_iter())
        print(pqr, '|Aut| =', naut)
    from math import factorial
    print('number of labelled copies 12!/|Aut| :', {pqr: factorial(12) // n for pqr, n in zip(PQR, [4, 2, 2, 4, 4])})
