"""Tight states reached by eliminations: compare reachable X with G_X tight ((D-1)|E| = D|V| - D - 1,
|V| >= 3) using S moves only versus S and star moves (Y0, Y1).  Claim [drv, REPORT 1.2]: for D >= 3 no
star move leads to a new tight state."""
import sys
import time
from collections import deque
from graphs_enum import admissible
from reduce_search import Network, _pc


def tight_states(net, D, kinds):
    start = 0
    seen = {start}
    dq = deque([start])
    out = set()
    while dq:
        X = dq.popleft()
        nb, src = net.state(X)
        nv = len(nb)
        ne = sum(_pc(m) for m in nb.values()) // 2
        if nv >= 3 and (D - 1) * ne == D * nv - D - 1:
            out.add(X)
        if nv <= 2:
            continue
        for v, m in nb.items():
            d = _pc(m)
            if (d == 2 and "S" in kinds) or (d == 3 and "Y" in kinds):
                Y = X | (1 << v)
                if Y not in seen:
                    seen.add(Y)
                    dq.append(Y)
    return out


for D, n in [tuple(map(int, a.split(","))) for a in sys.argv[1:]]:
    t0 = time.time()
    extra = 0
    graphs_with_extra = 0
    bsizes = {}
    for G in admissible(D, n):
        net = Network(G, [])
        a = tight_states(net, D, ("S",))
        b = tight_states(net, D, ("S", "Y"))
        if b - a:
            graphs_with_extra += 1
            extra += len(b - a)
        for X in b - a:
            # boundary sizes of the components of G[X]
            import networkx as nx
            nodes = [v for v in G if X >> v & 1]
            H = G.subgraph(nodes)
            key = []
            for C in nx.connected_components(H):
                bd = set(w for v in C for w in G.neighbors(v) if not X >> w & 1)
                key.append(len(bd))
            key = tuple(sorted(key))
            bsizes[key] = bsizes.get(key, 0) + 1
    print("D=%d n=%d: graphs with tight states reachable only through stars: %d (states %d); boundary sizes of the components of G[X]: %s | %.1fs" % (
        D, n, graphs_with_extra, extra, bsizes, time.time() - t0), flush=True)
