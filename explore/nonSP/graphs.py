"""The six K4 subdivisions (C^3, n=3) and prism, K33 (C^2, n=4) from results/graph_counts.txt."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'code'))


def __getattr__(name):
    if name == "Graph":
        from htk.graph import Graph
        return Graph
    raise AttributeError(name)

# edge lists copied from results/graph_counts.txt (d=3 n=3 "not reduced: K4" lines)
K4SUB = [
 ("K4[1,1,1,1,3,3]",  [(0, 4), (0, 7), (1, 2), (1, 5), (1, 7), (2, 3), (2, 4), (3, 6), (4, 5), (5, 6)]),
 ("K4[1,1,1,2,2,3]a", [(0, 4), (0, 7), (1, 2), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4), (5, 6), (5, 7)]),
 ("K4[1,1,1,2,2,3]b", [(0, 4), (0, 7), (1, 2), (1, 4), (1, 5), (2, 3), (2, 6), (3, 4), (5, 6), (6, 7)]),
 ("K4[1,1,1,2,2,3]c", [(0, 4), (0, 7), (1, 2), (1, 5), (2, 3), (2, 4), (3, 6), (3, 7), (4, 5), (5, 6)]),
 ("K4[1,1,2,2,2,2]a", [(0, 1), (0, 4), (0, 7), (1, 2), (2, 3), (2, 5), (3, 6), (3, 7), (4, 5), (5, 6)]),
 ("K4[1,1,2,2,2,2]b", [(0, 1), (0, 5), (0, 6), (1, 2), (1, 7), (2, 3), (3, 4), (3, 6), (4, 5), (4, 7)]),
]
# d=2 n=4
PRISM = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]
K33 = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (4, 5)]

def deg(edges, nv):
    d = [0]*nv
    for a, b in edges:
        d[a] += 1; d[b] += 1
    return d

def skeleton_chains(edges, nv):
    """K4 skeleton: cubic vertices, chains (list of vertex paths between cubic vertices)."""
    d = deg(edges, nv)
    adj = {v: [] for v in range(nv)}
    for e, (a, b) in enumerate(edges):
        adj[a].append((b, e)); adj[b].append((a, e))
    cubic = [v for v in range(nv) if d[v] == 3]
    used = set(); chains = []
    for p in cubic:
        for (w, e) in adj[p]:
            if e in used: continue
            path = [p]; es = [e]; used.add(e); prev = p; cur = w
            while d[cur] == 2:
                path.append(cur)
                (n1, e1), (n2, e2) = adj[cur]
                nxt, en = (n2, e2) if e1 == es[-1] else (n1, e1)
                used.add(en); es.append(en); prev, cur = cur, nxt
            path.append(cur)
            chains.append((path, es))
    return cubic, chains

__all__ = ["Graph", "K4SUB", "PRISM", "K33", "deg", "skeleton_chains", "sys", "os"]

if __name__ == "__main__":
    for name, E in K4SUB:
        cub, ch = skeleton_chains(E, 8)
        print(name, "cubic", cub, "deg2", [v for v in range(8) if deg(E, 8)[v] == 2])
        for path, es in ch:
            print("   chain", path, "edges", es)
