"""Graphs of Q_1 and Q_2 for holomorphic Chern-Simons on C^{4|1} (H = 4, T = 0, D = 4).

Labelling (insertion O = root = last vertex, as in S2 Section 5.2):
  pentagon   O -> 1 -> 2 -> 3 -> 4 -> O, outputs at 1..4 (vertices 0..3), O = 4;
             edges e0 = (O,1), e1 = (1,2), e2 = (2,3), e3 = (3,4), e4 = (4,O);
             displacements z_{e0} = z_1, z_{e4} = -z_2 (S2 Fig. fig-box).
  theta(p,q,r)  poles O (root) and N (vertex 0); paths P_1, P_2, P_3 from O to N with p, q, r edges;
             interior vertices numbered 1, 2, ... along P_1 (from O), then P_2, then P_3;
             edges oriented away from O, in the order P_1, P_2, P_3;
             spoke of P_j = its first edge, displacement z_j on it.
"""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from htk.graph import Graph


def pentagon():
    G = Graph(5, [(4, 0), (0, 1), (1, 2), (2, 3), (3, 4)], "pentagon")
    spokes = {0: (1, +1), 4: (2, -1)}          # edge -> (input label j, sign of z_j on it)
    outputs = [0, 1, 2, 3]
    return G, spokes, outputs


def theta_O(p, q, r):
    V = p + q + r - 1
    O, N = V - 1, 0
    edges, paths, spokes, nxt = [], [], {}, 1
    for jpath, L in enumerate((p, q, r)):
        prev, ints, eids = O, [], []
        for j in range(L):
            if j == L - 1:
                cur = N
            else:
                cur = nxt
                nxt += 1
                ints.append(cur)
            if j == 0:
                spokes[len(edges)] = (jpath + 1, +1)
            eids.append(len(edges))
            edges.append((prev, cur))
            prev = cur
        paths.append((eids, ints))
    G = Graph(V, edges, "theta_O(%d,%d,%d)" % (p, q, r))
    outputs = [v for _e, ints in paths for v in ints]
    return G, spokes, outputs, paths


THETAS = [(1, 4, 4), (2, 3, 4), (3, 3, 3)]
AUT = {(1, 4, 4): 2, (2, 3, 4): 1, (3, 3, 3): 6, 'pent': 2}
