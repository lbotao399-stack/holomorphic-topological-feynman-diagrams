"""Admissible graphs on R^T x C^H, D = H + T.

|I| = D n + 1, |V| = (D-1) n + 2, f(S) = D|V_S| - D - 1 - (D-1)|I_S| >= 0 for every vertex set S, |S| >= 2.
f(whole) = 0, no bridge, no cut vertex (f < 0 otherwise), so Gamma = subdivision of a 2-connected
loopless multigraph skeleton, minimum degree >= 3 (n >= 2), sum_v (deg v - 2) = 2n - 2.
Local check on skeleton sets Q: a chain of length l with both ends in Q contributes 1 - D (l = 1) or
min(0, l - D) (l >= 2); partial chains raise f.

D >= 3 (and D = 2, n <= 4): skeleton enumeration.  D = 2: htk.catalog.enumerate_tight (networkx atlas).
"""
import itertools
import sys
import os

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(P, "code"))


def skeletons(n):
    """2-connected loopless multigraphs, min degree >= 3, cycle rank n >= 2, up to isomorphism.
    Returns list of (s, edge list)."""
    out = []
    for s in range(2, 2 * (n - 1) + 1):
        m = s + n - 1
        if 2 * m < 3 * s:
            continue
        pairs = list(itertools.combinations(range(s), 2))
        found = {}
        for ms in itertools.combinations_with_replacement(pairs, m):
            deg = [0] * s
            for a, b in ms:
                deg[a] += 1
                deg[b] += 1
            if min(deg) < 3:
                continue
            Gs = nx.Graph()
            Gs.add_nodes_from(range(s))
            Gs.add_edges_from(ms)
            if not nx.is_connected(Gs):
                continue
            if s >= 3 and not nx.is_biconnected(Gs):
                continue
            MG = nx.MultiGraph()
            MG.add_nodes_from(range(s))
            MG.add_edges_from(ms)
            key = tuple(sorted(deg))
            lst = found.setdefault(key, [])
            if any(nx.is_isomorphic(MG, H) for H in lst):
                continue
            lst.append(MG)
            out.append((s, list(ms)))
    return out


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def local_ok(D, s, sk, lens):
    for q in range(2, s + 1):
        for Q in itertools.combinations(range(s), q):
            Qs = set(Q)
            f = D * q - D - 1
            for (a, b), l in zip(sk, lens):
                if a in Qs and b in Qs:
                    f += (1 - D) if l == 1 else min(0, l - D)
            if f < 0:
                return False
    return True


def subdivide(s, sk, lens):
    nv = s
    E = []
    for (a, b), l in zip(sk, lens):
        prev = a
        for j in range(l):
            if j == l - 1:
                nxt = b
            else:
                nxt = nv
                nv += 1
            E.append((prev, nxt))
            prev = nxt
    return nv, E


def sparse_ok_bruteforce(G, D):
    """f(S) >= 0 for every vertex subset with an edge (exponential; |V| <= 18)."""
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * len(nodes)
    for a, b in G.edges():
        adj[idx[a]] |= 1 << idx[b]
        adj[idx[b]] |= 1 << idx[a]
    for S in range(3, 1 << len(nodes)):
        k = bin(S).count("1")
        if k < 2:
            continue
        e2, T, i = 0, S, 0
        while T:
            if T & 1:
                e2 += bin(adj[i] & S).count("1")
            T >>= 1
            i += 1
        e = e2 // 2
        if e and (D - 1) * e > D * k - D - 1:
            return False
    return True


def admissible_skeleton(D, n):
    """networkx graphs (simple), up to isomorphism, by the skeleton method."""
    N = D * n + 1
    if n == 1:
        return [nx.cycle_graph(N)]
    buckets = {}
    out = []
    for s, sk in skeletons(n):
        m = len(sk)
        for lens in compositions(N, m):
            # simple graph: among parallel skeleton edges at most one of length 1
            seen1 = set()
            bad = False
            for (a, b), l in zip(sk, lens):
                if l == 1:
                    key = (min(a, b), max(a, b))
                    if key in seen1:
                        bad = True
                        break
                    seen1.add(key)
            if bad or not local_ok(D, s, sk, lens):
                continue
            nv, E = subdivide(s, sk, lens)
            G = nx.Graph()
            G.add_nodes_from(range(nv))
            G.add_edges_from(E)
            assert G.number_of_edges() == N and nv == (D - 1) * n + 2
            h = (tuple(sorted(d for _, d in G.degree())), nx.weisfeiler_lehman_graph_hash(G, iterations=4))
            lst = buckets.setdefault(h, [])
            if any(nx.is_isomorphic(G, H) for H in lst):
                continue
            lst.append(G)
            out.append(G)
    return out


def admissible(D, n, method="auto"):
    if method == "auto":
        method = "catalog" if (D == 2 and n >= 5) else "skeleton"
    if method == "catalog":
        from htk import catalog
        return [nx.convert_node_labels_to_integers(G) for G in catalog.enumerate_tight(D, n)]
    return admissible_skeleton(D, n)


if __name__ == "__main__":
    import time
    for D, n in [(2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2)]:
        t0 = time.time()
        a = admissible_skeleton(D, n)
        from htk import catalog
        b = catalog.enumerate_tight(D, n)
        print(D, n, len(a), len(b), "%.1fs" % (time.time() - t0))
