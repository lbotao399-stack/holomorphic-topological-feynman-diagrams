"""Circuit reduction of an admissible graph with a source pattern: move search and cores.

Every move eliminates one vertex v (Kron reduction):
    S   deg v = 2                      series (+ parallel merge)      polynomial transfer (any source)
    Y0  deg v = 3, v source-free       star-triangle (+ merges)        polynomial
    Y1  deg v = 3, v sourced           star-triangle with source       rational (lambda_u += R_wx/P lambda_v)
After eliminating a vertex set X the network on V \\ X is the elimination graph G_X:
    original edges among V \\ X  +  a clique on the boundary of every component C of G[X];
a vertex u of V \\ X is sourced iff u in S0 or some component C with C n S0 != {} has u on its boundary
(harmonic measure > 0).  Hence the state of a reduction is the set X alone (order independence of Kron
reduction); validity (deg 2 or 3 at elimination) and the source status of v at its elimination depend on X.

Tight moves (definition A): tight chain moves (k >= 2 chains p..q, interior vertices bivalent,
sum a_i = D (k - 1) + 1), and for D = 2 the LLC moves (b) single source-free star (ux, wx present, uw
absent) and (c) paired source-free stars (v1, v2 joined to u, A, B; AB present; uA, uB absent).
"""
import itertools
from collections import deque

import networkx as nx


# ---------------------------------------------------------------------------- elimination graphs

class Network:
    """Original simple graph (vertices 0..nv-1) with an initial source set S0."""

    def __init__(self, G, S0):
        self.G = nx.convert_node_labels_to_integers(G) if set(G.nodes()) != set(range(G.number_of_nodes())) else G
        self.nv = self.G.number_of_nodes()
        self.adj = [0] * self.nv
        for a, b in self.G.edges():
            self.adj[a] |= 1 << b
            self.adj[b] |= 1 << a
        self.S0 = 0
        for v in S0:
            self.S0 |= 1 << v
        self.full = (1 << self.nv) - 1

    def state(self, X):
        """(neighbour bitmasks of the remaining vertices, sourced bitmask) of G_X."""
        nb = {}
        rem = self.full & ~X
        r = rem
        while r:
            v = (r & -r).bit_length() - 1
            r &= r - 1
            nb[v] = self.adj[v] & rem
        src = self.S0 & rem
        seen = 0
        x = X
        while x:
            v = (x & -x).bit_length() - 1
            x &= x - 1
            if seen >> v & 1:
                continue
            comp, stack, has_src, bdry = 1 << v, [v], bool(self.S0 >> v & 1), 0
            while stack:
                u = stack.pop()
                for w in _bits(self.adj[u]):
                    if X >> w & 1:
                        if not comp >> w & 1:
                            comp |= 1 << w
                            stack.append(w)
                            if self.S0 >> w & 1:
                                has_src = True
                    else:
                        bdry |= 1 << w
            seen |= comp
            for u in _bits(bdry):
                nb[u] |= bdry & ~(1 << u)
            if has_src:
                src |= bdry
        return nb, src


def _bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b


def _pc(m):
    return bin(m).count("1")


def moves(net, X, nb=None, src=None):
    """[(v, kind)] with kind in 'S', 'Y0', 'Y1'."""
    if nb is None:
        nb, src = net.state(X)
    out = []
    for v, m in nb.items():
        d = _pc(m)
        if d == 2:
            out.append((v, "S"))
        elif d == 3:
            out.append((v, "Y1" if src >> v & 1 else "Y0"))
    return out


CLOSURE = True


def safe_closure(net, X):
    """Eliminate bivalent vertices whose elimination changes no source status:
    deg 2 and (v source-free or both neighbours sourced).  Returns (X, nb, src, eliminated list).
    CLOSURE = False: no elimination (exhaustive search over all eliminated sets)."""
    elim = []
    if not CLOSURE:
        nb, src = net.state(X)
        return X, nb, src, elim
    while True:
        nb, src = net.state(X)
        hit = None
        for v, m in nb.items():
            if len(nb) > 2 and _pc(m) == 2 and (not src >> v & 1 or (m & ~src) == 0):
                hit = v
                break
        if hit is None:
            return X, nb, src, elim
        elim.append((hit, "S"))
        X |= 1 << hit


def search(net, allow=("S", "Y0"), max_states=2_000_000):
    """0-1 BFS over eliminated sets X.  Cost = number of Y1 moves.
    Returns dict: reducible, min_Y1 (None if irreducible), stuck (list of X of minimal remaining size),
    nstates, sequence (list of (v, kind)) of an optimal reduction."""
    allow = set(allow)
    start, _, _, el0 = safe_closure(net, 0)
    dist = {start: 0}
    parent = {start: None}
    dq = deque([start])
    done = set()
    stuck = []
    best_stuck = None
    target = None
    while dq:
        X = dq.popleft()
        if X in done:
            continue
        done.add(X)
        if len(done) > max_states:
            raise RuntimeError("state limit")
        nb, src = net.state(X)
        if len(nb) <= 2:
            target = X
            break
        mv = [(v, k) for v, k in moves(net, X, nb, src) if k in allow]
        if not mv:
            size = (len(nb), sum(_pc(m) for m in nb.values()) // 2)
            if best_stuck is None or size < best_stuck:
                best_stuck, stuck = size, [X]
            elif size == best_stuck:
                stuck.append(X)
            continue
        for v, k in mv:
            Y, _, _, el1 = safe_closure(net, X | (1 << v))
            c = dist[X] + (1 if k == "Y1" else 0)
            if Y not in dist or c < dist[Y]:
                dist[Y] = c
                parent[Y] = (X, [(v, k)] + el1)
                if k == "Y1":
                    dq.append(Y)
                else:
                    dq.appendleft(Y)
    res = {"reducible": target is not None, "nstates": len(done), "stuck": stuck,
           "min_Y1": dist[target] if target is not None else None, "sequence": None}
    if target is not None:
        seq, X = [], target
        while parent[X] is not None:
            X0, steps = parent[X]
            seq = steps + seq
            X = X0
        res["sequence"] = el0 + seq
    return res


def core_graph(net, X):
    """The elimination graph G_X as a networkx graph with node attribute src (bool)."""
    nb, src = net.state(X)
    H = nx.Graph()
    for v in nb:
        H.add_node(v, src=bool(src >> v & 1))
    for v, m in nb.items():
        for w in _bits(m):
            if v < w:
                H.add_edge(v, w)
    return H


# ---------------------------------------------------------------------------- tight moves (definition A)

def _chains_between(adjm, p, q):
    """chains from p to q (interior vertices of degree 2): list of interior tuples."""
    out = []
    for w0 in adjm[p]:
        prev, cur, interior = p, w0, []
        while True:
            if cur == q:
                out.append(tuple(interior))
                break
            if cur == p or len(adjm[cur]) != 2:
                break
            interior.append(cur)
            nxt = [u for u in adjm[cur] if u != prev]
            prev, cur = cur, nxt[0]
    return out


def tight_search(G, S0, D, sourced_star=False, limit=200000, paired=True):
    """Definition A.  State: (frozenset of edges (a<b), frozenset sources).  Moves: tight chain moves;
    D = 2: LLC (b) single source-free star, (c) paired source-free stars; sourced_star: also (b) with a
    sourced centre (LLC App. E, rational).  Returns dict reducible, n_sourced (minimal), stuck (minimal
    states), nstates."""
    E0 = frozenset(tuple(sorted(e)) for e in G.edges())
    s0 = frozenset(S0)
    start = (E0, s0)
    dist = {start: 0}
    parent = {start: None}
    dq = deque([start])
    done = set()
    best, stuck = None, []
    target = None
    while dq:
        st = dq.popleft()
        if st in done:
            continue
        done.add(st)
        if len(done) > limit:
            raise RuntimeError("state limit")
        E, src = st
        adjm = {}
        for a, b in E:
            adjm.setdefault(a, set()).add(b)
            adjm.setdefault(b, set()).add(a)
        if len(adjm) == 2:
            target = st
            break
        nxt = []
        verts = sorted(adjm)
        for p, q in itertools.combinations(verts, 2):
            ch = _chains_between(adjm, p, q)
            if len(ch) < 2:
                continue
            for k in range(2, len(ch) + 1):
                for sub in itertools.combinations(ch, k):
                    if sum(len(c) + 1 for c in sub) != D * (k - 1) + 1:
                        continue
                    interior = [v for c in sub for v in c]
                    E2 = set(E)
                    for c in sub:
                        path = (p,) + c + (q,)
                        for a, b in zip(path, path[1:]):
                            E2.discard(tuple(sorted((a, b))))
                    assert tuple(sorted((p, q))) not in E2
                    E2.add(tuple(sorted((p, q))))
                    s2 = set(src) - set(interior)
                    if any(v in src for v in interior):
                        s2 |= {p, q}
                    nxt.append(((frozenset(E2), frozenset(s2)), 0, [("S", v) for v in interior]))
        if D == 2:
            for v in verts:
                nbv = adjm[v]
                if len(nbv) != 3:
                    continue
                sv = v in src
                if sv and not sourced_star:
                    continue
                for u, w in itertools.combinations(sorted(nbv), 2):
                    x = (nbv - {u, w}).pop()
                    if w in adjm[u] or x not in adjm[u] or x not in adjm[w]:
                        continue
                    E2 = set(E) - {tuple(sorted((v, a))) for a in nbv}
                    E2.add(tuple(sorted((u, w))))
                    s2 = set(src) - {v}
                    if sv:
                        s2 |= nbv
                    nxt.append(((frozenset(E2), frozenset(s2)), 1 if sv else 0, [("Y", v)]))
            for v1, v2 in (itertools.combinations(verts, 2) if paired else []):
                if v1 in src or v2 in src:
                    continue
                if len(adjm[v1]) != 3 or adjm[v1] != adjm[v2]:
                    continue
                N = sorted(adjm[v1])
                for u in N:
                    A, B = [a for a in N if a != u]
                    if B not in adjm[A] or A in adjm[u] or B in adjm[u]:
                        continue
                    E2 = set(E) - {tuple(sorted((a, b))) for a in (v1, v2) for b in N}
                    E2 |= {tuple(sorted((u, A))), tuple(sorted((u, B)))}
                    nxt.append(((frozenset(E2), frozenset(set(src) - {v1, v2})), 0, [("Y", v1), ("Y", v2)]))
        if not nxt:
            size = (len(adjm), len(E))
            if best is None or size < best:
                best, stuck = size, [st]
            elif size == best:
                stuck.append(st)
            continue
        for st2, c, steps in nxt:
            d2 = dist[st] + c
            if st2 not in dist or d2 < dist[st2]:
                dist[st2] = d2
                parent[st2] = (st, steps)
                if c:
                    dq.append(st2)
                else:
                    dq.appendleft(st2)
    def path(st):
        seq = []
        while parent[st] is not None:
            st0, steps = parent[st]
            seq = steps + seq
            st = st0
        return seq
    return {"reducible": target is not None, "n_sourced": dist[target] if target else None,
            "stuck": stuck, "nstates": len(done), "path": path(target) if target else None,
            "stuck_path": path(stuck[0]) if stuck else None}


def tight_core_graph(st):
    E, src = st
    H = nx.Graph()
    for a, b in E:
        H.add_edge(a, b)
    for v in H.nodes():
        H.nodes[v]["src"] = v in src
    return H


# ---------------------------------------------------------------------------- classes up to isomorphism

class Classes:
    """Isomorphism classes of graphs with a boolean node attribute src."""

    def __init__(self):
        self.reps = []      # (key, graph, count, members)
        self._buckets = {}

    @staticmethod
    def _key(H):
        degs = tuple(sorted((d, bool(H.nodes[v].get("src", False))) for v, d in H.degree()))
        Hl = H.copy()
        for v in Hl.nodes():
            Hl.nodes[v]["lab"] = "s" if Hl.nodes[v].get("src", False) else "o"
        return (degs, nx.weisfeiler_lehman_graph_hash(Hl, node_attr="lab", iterations=4))

    def add(self, H, tag=None):
        k = self._key(H)
        lst = self._buckets.setdefault(k, [])
        nm = nx.algorithms.isomorphism.categorical_node_match("src", False)
        for idx in lst:
            if nx.is_isomorphic(self.reps[idx][1], H, node_match=nm):
                self.reps[idx][2] += 1
                self.reps[idx][3].append(tag)
                return idx
        self.reps.append([k, H, 1, [tag]])
        lst.append(len(self.reps) - 1)
        return len(self.reps) - 1


def describe_core(H):
    """short name: |V|, |E|, degree/src pattern, recognised names."""
    nvh, neh = H.number_of_nodes(), H.number_of_edges()
    U = nx.Graph(H)
    name = None
    if nvh == 2 and neh == 1:
        name = "edge"
    elif nvh == 4 and neh == 6:
        name = "K4"
    elif nvh == 5 and neh == 10:
        name = "K5"
    elif nvh == 6 and neh == 9 and nx.is_isomorphic(U, nx.complete_bipartite_graph(3, 3)):
        name = "K33"
    elif nvh == 6 and neh == 9 and nx.is_isomorphic(U, nx.circular_ladder_graph(3)):
        name = "prism"
    elif nvh == 6 and neh == 12 and nx.is_isomorphic(U, nx.octahedral_graph()):
        name = "octahedron"
    elif nvh == 5 and neh == 8 and nx.is_isomorphic(U, nx.wheel_graph(5)):
        name = "W4"
    ns = sum(1 for v in H if H.nodes[v].get("src", False))
    degs = sorted((d for _, d in H.degree()), reverse=True)
    return "%s|V|=%d |E|=%d src=%d deg=%s" % ((name + " ") if name else "", nvh, neh, ns,
                                               "".join(map(str, degs)))
