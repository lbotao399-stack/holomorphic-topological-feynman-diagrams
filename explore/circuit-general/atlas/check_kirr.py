"""Independent check of the I class (D = 2, n = 6): explicit Kron elimination on adjacency sets
(eliminate v: delete v, join its neighbours pairwise), DFS over eliminated sets, valid iff deg in {2,3}.
Reducible iff two vertices remain.  Compared with census.json (k_gen = None)."""
import json
import networkx as nx
from graphs_enum import admissible


def reducible(G):
    adj0 = {v: set(G.neighbors(v)) for v in G}
    seen = set()

    def rec(adj):
        if len(adj) <= 2:
            return True
        key = frozenset(adj)
        if key in seen:
            return False
        seen.add(key)
        for v in list(adj):
            if len(adj[v]) in (2, 3):
                nb = adj[v]
                new = {u: set(w for w in adj[u] if w != v) for u in adj if u != v}
                for a in nb:
                    for b in nb:
                        if a != b:
                            new[a].add(b)
                if rec(new):
                    return True
        return False
    return rec(adj0)


data = json.load(open("census.json"))["data"]["2,6"]["graphs"]
gs = admissible(2, 6)
irr = [i for i, G in enumerate(gs) if not reducible(G)]
census_irr = [r["index"] for r in data if r["k_gen"] is None]
print("D=2 n=6: irreducible by explicit Kron DFS: %d graphs; census I: %d; same set: %s" % (
    len(irr), len(census_irr), sorted(irr) == sorted(census_irr)))
k5 = 0
for i in irr:
    G = gs[i]
    print("  graph %d: degrees %s, planar %s, contains K5 minor-free check skipped" % (
        i, sorted((d for _, d in G.degree()), reverse=True), nx.check_planarity(G)[0]))
# A-core composition (gen) at D=2 n=6
cls = json.load(open("census.json"))["data"]["2,6"]["classes"]["A"]
from collections import Counter
c = Counter()
for k in cls:
    name = k["desc"].split("|V|")[0].strip() or ("|V|" + k["desc"].split("|V|")[1].split(" src")[0])
    c[name] += k["count"]
print("A-cores D=2 n=6 by topology:", dict(c))
npl = [i for i, G in enumerate(gs) if not nx.check_planarity(G)[0]]
print("D=2 n=6: non-planar graphs %d; irreducible ones all non-planar: %s; non-planar and reducible: %d" % (
    len(npl), set(irr) <= set(npl), len(set(npl) - set(irr))))
for n in (4, 5):
    g2 = admissible(2, n)
    print("D=2 n=%d: non-planar %d, irreducible %d" % (n, sum(1 for G in g2 if not nx.check_planarity(G)[0]),
                                                     sum(1 for G in g2 if not reducible(G))))
