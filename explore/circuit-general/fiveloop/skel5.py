"""Five-loop SP-irreducible skeletons: simple 2-connected graphs, minimum degree >= 3,
cycle rank 5 (|V| = 5..8, |E| = |V| + 4).  3-connectivity, planarity, 2-vertex cuts,
planar duals.  Also: series-parallel reduction of an arbitrary multigraph to its skeleton."""
import itertools
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))


def sp_skeleton(edges):
    """series-parallel reduction (remove degree-2 vertices, merge parallel edges, drop
    resulting self-loops never occur for 2-connected input); returns simple edge list."""
    E = [tuple(e) for e in edges]
    while True:
        # parallel
        seen = {}
        newE = []
        changed = False
        for u, v in E:
            k = (min(u, v), max(u, v))
            if k in seen:
                changed = True
                continue
            seen[k] = 1
            newE.append(k)
        E = newE
        deg = {}
        for u, v in E:
            deg[u] = deg.get(u, 0) + 1
            deg[v] = deg.get(v, 0) + 1
        d2 = [v for v in deg if deg[v] == 2]
        if d2:
            v = d2[0]
            inc = [e for e in E if v in e]
            a = inc[0][0] if inc[0][1] == v else inc[0][1]
            b = inc[1][0] if inc[1][1] == v else inc[1][1]
            E = [e for e in E if v not in e]
            if a != b:
                E.append((min(a, b), max(a, b)))
            changed = True
        if not changed:
            return E


def canon(edges):
    G = nx.Graph(list(edges))
    return G


def skeletons5():
    out = []
    for s in range(5, 9):
        m = s + 4
        pairs = list(itertools.combinations(range(s), 2))
        found = []
        for es in itertools.combinations(pairs, m):
            deg = [0] * s
            for a, b in es:
                deg[a] += 1
                deg[b] += 1
            if min(deg) < 3:
                continue
            G = nx.Graph()
            G.add_nodes_from(range(s))
            G.add_edges_from(es)
            if not nx.is_biconnected(G):
                continue
            if any(nx.is_isomorphic(G, H) for H in found):
                continue
            found.append(G)
        out.extend(found)
    return out


def describe(G):
    k = nx.node_connectivity(G)
    planar, emb = nx.check_planarity(G)
    degs = tuple(sorted((d for _, d in G.degree()), reverse=True))
    return dict(V=G.number_of_nodes(), E=G.number_of_edges(), degs=degs, conn=k, planar=planar,
                girth=nx.girth(G) if hasattr(nx, 'girth') else None,
                tri=sum(nx.triangles(G).values()) // 3,
                aut=len(list(nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())))


def planar_dual(G):
    """dual multigraph of a 2-connected planar graph: list of (face_a, face_b) per edge of G
    (edge order = list(G.edges()))."""
    ok, emb = nx.check_planarity(G)
    assert ok
    faces = {}
    fid = 0
    visited = set()
    for u, v in emb.edges():
        if (u, v) in visited:
            continue
        face = emb.traverse_face(u, v, mark_half_edges=visited)
        # half edges of this face
        hs = [(face[i], face[(i + 1) % len(face)]) for i in range(len(face))]
        for h in hs:
            faces[h] = fid
        fid += 1
    dual = []
    for u, v in G.edges():
        dual.append((faces[(u, v)], faces[(v, u)]))
    return list(G.edges()), dual


NAMES = {}


def name_of(G):
    d = describe(G)
    V, degs = d['V'], d['degs']
    if V == 5:
        return 'K5-e'
    if V == 6:
        if degs[0] == 5:
            return 'W5'
        if d['conn'] < 3:
            return 'K4+K4'           # 2-sum of two K4
        if not d['planar']:
            return 'K33+e'
        return 'Z5'
    return None


if __name__ == '__main__':
    import json
    sk = skeletons5()
    rows = []
    for i, G in enumerate(sk):
        d = describe(G)
        nm = name_of(G)
        rows.append(dict(index=i, name=nm, edges=sorted(tuple(sorted(e)) for e in G.edges()), **d))
    # names for V = 7, 8
    c7 = [r for r in rows if r['V'] == 7]
    c8 = [r for r in rows if r['V'] == 8]
    for j, r in enumerate(c7):
        r['name'] = r['name'] or 'V7_%d' % j
    for j, r in enumerate(c8):
        r['name'] = r['name'] or 'V8_%d' % j
    for r in rows:
        print('%-7s V=%d E=%d degs=%s conn=%d planar=%s tri=%d |Aut|=%d  %s' % (
            r['name'], r['V'], r['E'], ''.join(map(str, r['degs'])), r['conn'], r['planar'], r['tri'], r['aut'],
            r['edges']))
    json.dump(rows, open(os.path.join(HERE, 'skel5.json'), 'w'), indent=0)
