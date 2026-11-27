"""Own graph layer for C4: skeleton multigraphs, admissible subdivisions, cycle bases, period matrices."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
import random
import networkx as nx


def skeletons(n):
    """2-connected loopless multigraphs with min degree >= 3 and cyclomatic number n >= 2, up to isomorphism"""
    out = []
    for s in range(2, 2 * (n - 1) + 1):
        m = s + n - 1
        if 2 * m < 3 * s:
            continue
        pairs = list(combinations(range(s), 2))
        found = []
        for ms in combinations_with_replacement(pairs, m):
            deg = [0] * s
            for a, b in ms:
                deg[a] += 1
                deg[b] += 1
            if min(deg) < 3:
                continue
            G = nx.Graph()
            G.add_nodes_from(range(s))
            G.add_edges_from(ms)
            if not nx.is_connected(G) or (s >= 3 and not nx.is_biconnected(G)):
                continue
            MG = nx.MultiGraph()
            MG.add_nodes_from(range(s))
            MG.add_edges_from(ms)
            if any(nx.is_isomorphic(MG, H) for H in found):
                continue
            found.append(MG)
            out.append((s, list(ms)))
    return out


def random_cubic_skeleton(n, rng, tries=1000):
    """random 2-connected loopless cubic multigraph with cyclomatic number n (s = 2(n-1) vertices)"""
    s = 2 * (n - 1)
    for _ in range(tries):
        half = [v for v in range(s) for _ in range(3)]
        rng.shuffle(half)
        E = [(min(half[2 * i], half[2 * i + 1]), max(half[2 * i], half[2 * i + 1])) for i in range(len(half) // 2)]
        if any(a == b for a, b in E):
            continue
        G = nx.Graph()
        G.add_nodes_from(range(s))
        G.add_edges_from(E)
        if nx.is_connected(G) and (s < 3 or nx.is_biconnected(G)):
            return s, sorted(E)
    raise RuntimeError


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def local_ok(D, s, sk, lens):
    """min over skeleton vertex sets Q of f = D|Q| - D - 1 + sum_{chains inside Q} (l==1 ? 1-D : min(0, l-D)) >= 0"""
    for q in range(2, s + 1):
        for Q in combinations(range(s), q):
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


def sparse_ok(nv, E, D):
    """direct check (D-1)|I_S| <= D|V_S| - D - 1 for every vertex set S with >= 2 vertices (small graphs)"""
    for k in range(2, nv + 1):
        for S in combinations(range(nv), k):
            Ss = set(S)
            e = sum(1 for a, b in E if a in Ss and b in Ss)
            if e and (D - 1) * e > D * k - D - 1:
                return False
    return True


def admissible(D, n, dedupe=True):
    """connected simple graphs with |I| = Dn+1, |V| = (D-1)n+2, (D-1)|I_S| <= D|V_S|-D-1 (D >= 2)"""
    N = D * n + 1
    if n == 1:
        return [(N, [(i, (i + 1) % N) for i in range(N)])]
    cands = []
    for s, sk in skeletons(n):
        for lens in compositions(N, len(sk)):
            if not local_ok(D, s, sk, lens):
                continue
            nv, E = subdivide(s, sk, lens)
            if len(set(tuple(sorted(e)) for e in E)) < len(E):
                continue
            cands.append((nv, E))
    if not dedupe:
        return cands
    buckets, out = {}, []
    for nv, E in cands:
        G = nx.Graph()
        G.add_nodes_from(range(nv))
        G.add_edges_from(E)
        h = nx.weisfeiler_lehman_graph_hash(G)
        lst = buckets.setdefault(h, [])
        if any(nx.is_isomorphic(G, H) for H in lst):
            continue
        lst.append(G)
        out.append((nv, sorted(tuple(sorted(e)) for e in G.edges())))
    return out


def cycle_matrix(nv, E):
    """fundamental cycles of a Kruskal spanning tree (edge order); rows over cotree edges; works for multigraphs"""
    parent = list(range(nv))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    tree, cot = [], []
    for e, (a, b) in enumerate(E):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
            tree.append(e)
        else:
            cot.append(e)
    adj = {v: [] for v in range(nv)}
    for e in tree:
        a, b = E[e]
        adj[a].append((b, e, +1))
        adj[b].append((a, e, -1))

    def path(u, w):
        """edges with signs on the tree path u -> w"""
        prev = {u: None}
        stack = [u]
        while stack:
            x = stack.pop()
            for y, e, sgn in adj[x]:
                if y not in prev:
                    prev[y] = (x, e, sgn)
                    stack.append(y)
        out = []
        x = w
        while prev[x] is not None:
            px, e, sgn = prev[x]
            out.append((e, sgn))
            x = px
        return out
    C = []
    for c in cot:
        a, b = E[c]
        row = [0] * len(E)
        row[c] = 1
        for e, sgn in path(b, a):     # back from b to a through the tree
            row[e] += sgn
        C.append(row)
    return C


def incidence(nv, E):
    """B[e][v] for v < root = nv-1"""
    return [[(1 if v == a else 0) - (1 if v == b else 0) for v in range(nv - 1)] for (a, b) in E]


def check_cycles(nv, E, C):
    B = incidence(nv, E)
    return all(sum(C[r][e] * B[e][v] for e in range(len(E))) == 0 for r in range(len(C)) for v in range(nv - 1))
