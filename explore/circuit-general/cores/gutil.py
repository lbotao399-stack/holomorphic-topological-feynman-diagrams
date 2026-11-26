"""Small graph utilities: cycle rank, projective convergence, spanning trees, automorphisms."""
from functools import lru_cache
from itertools import combinations, permutations


def cycle_rank(edges):
    """loop number of the subgraph given by a list of (u,v) edges (multigraph allowed)"""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    loops = 0
    for u, v in edges:
        a, b = find(u), find(v)
        if a == b:
            loops += 1
        else:
            parent[a] = b
    return loops


@lru_cache(maxsize=None)
def subset_ranks(edges):
    edges = tuple(edges)
    N = len(edges)
    out = []
    for mask in range(1, (1 << N) - 1):
        sub = [edges[i] for i in range(N) if mask >> i & 1]
        out.append((mask, cycle_rank(sub)))
    return tuple(out)


@lru_cache(maxsize=None)
def subset_table(edges):
    """for each proper nonempty mask: (mask, h); plus lowbit decomposition for fast sums"""
    edges = tuple(edges)
    N = len(edges)
    hs = [0] * (1 << N)
    for mask in range(1, 1 << N):
        sub = [edges[i] for i in range(N) if mask >> i & 1]
        hs[mask] = cycle_rank(sub)
    return N, tuple(hs)


_conv_cache = {}


def convergent(edges, powers, q):
    """sum_{e in S} a_e > q h(S) for every proper nonempty S, and total = n q"""
    edges = tuple(edges)
    key = (edges, tuple(powers), q)
    r = _conv_cache.get(key)
    if r is not None:
        return r
    N, hs = subset_table(edges)
    full = (1 << N) - 1
    if sum(powers) != hs[full] * q:
        _conv_cache[key] = False
        return False
    sums = [0] * (1 << N)
    ok = True
    for mask in range(1, full):
        low = mask & -mask
        i = low.bit_length() - 1
        s = sums[mask ^ low] + powers[i]
        sums[mask] = s
        if s <= q * hs[mask]:
            ok = False
            break
    if len(_conv_cache) > 2000000:
        _conv_cache.clear()
    _conv_cache[key] = ok
    return ok


@lru_cache(maxsize=None)
def spanning_trees(edges):
    edges = tuple(edges)
    verts = sorted({x for e in edges for x in e})
    k = len(verts) - 1
    out = []
    for comb in combinations(range(len(edges)), k):
        if cycle_rank([edges[i] for i in comb]) == 0:
            out.append(comb)
    return tuple(out)


@lru_cache(maxsize=None)
def cotrees(edges):
    edges = tuple(edges)
    N = len(edges)
    return tuple(tuple(i for i in range(N) if i not in t) for t in spanning_trees(edges))


@lru_cache(maxsize=None)
def edge_automorphisms(edges):
    """list of edge permutations induced by vertex automorphisms (simple graphs)"""
    edges = tuple(edges)
    verts = sorted({x for e in edges for x in e})
    eset = {frozenset(e): i for i, e in enumerate(edges)}
    out = []
    for perm in permutations(verts):
        mp = dict(zip(verts, perm))
        img = []
        ok = True
        for u, v in edges:
            j = eset.get(frozenset((mp[u], mp[v])))
            if j is None:
                ok = False
                break
            img.append(j)
        if ok:
            out.append(tuple(img))
    return tuple(sorted(set(out)))


def orbit_rep(powers, autos):
    return min(tuple(powers[p[i]] for i in range(len(powers))) for p in autos)


def compositions(total, parts, lo, hi):
    if parts == 0:
        if total == 0:
            yield ()
        return
    for x in range(max(lo, total - hi * (parts - 1)), min(hi, total - lo * (parts - 1)) + 1):
        for rest in compositions(total - x, parts - 1, lo, hi):
            yield (x,) + rest
