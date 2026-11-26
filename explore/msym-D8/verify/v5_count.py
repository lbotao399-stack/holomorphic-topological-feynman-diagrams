"""C4: graphs of Q_1, Q_2 on C^4 (D = 4): |I| = 4n+1, |V| = 3n+2, 3|I_S| <= 4|V_S| - 5 for every edge subset S.

(a) exhaustive over connected simple graphs with |V| = 3n+2, |I| = 4n+1 and minimum degree >= 1, n = 1 (all graphs)
    and n = 2 (all graphs: enumeration of edge sets is replaced by the 2-core argument checked here by brute force on
    every graph whose 2-core is a subdivided theta, dumbbell or figure eight with 9 edges; trees attached are excluded
    because the 2-core C must satisfy 3|C| <= 4(|C|-1) - 5);
(b) the theta rule deg O = n + p (Gamma \\ O a forest with p components), |Aut| fixing O, outputs = 2n+2.
"""
import itertools
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def sparse_ok(G):
    E = list(G.edges())
    for r in range(1, len(E) + 1):
        for S in itertools.combinations(E, r):
            V = set(v for e in S for v in e)
            if 3 * len(S) > 4 * len(V) - 5:
                return False, S
    return True, None


def aut_fixing(G, O):
    nx.set_node_attributes(G, {v: (v == O) for v in G.nodes()}, 'o')
    gm = GraphMatcher(G, G, node_match=lambda a, b: a['o'] == b['o'])
    return sum(1 for _ in gm.isomorphisms_iter())


# n = 1: all connected graphs with 5 vertices and 5 edges
cnt, good = 0, []
nodes = range(5)
alle = list(itertools.combinations(nodes, 2))
for S in itertools.combinations(alle, 5):
    G = nx.Graph(); G.add_nodes_from(nodes); G.add_edges_from(S)
    if not nx.is_connected(G):
        continue
    cnt += 1
    if sparse_ok(G)[0]:
        good.append(G)
iso = []
for G in good:
    if not any(nx.is_isomorphic(G, H) for H in iso):
        iso.append(G)
print('(a) n = 1: %d connected graphs (5 vertices, 5 edges); sparse: %d, isomorphism classes: %s'
      % (cnt, len(good), [sorted(d for _, d in G.degree()) for G in iso]))


def theta_graph(p, q, r):
    G = nx.Graph()
    O, N = 'O', 'N'
    k = 0
    for L in (p, q, r):
        prev = O
        for j in range(L):
            cur = N if j == L - 1 else 'v%d' % k
            if j < L - 1:
                k += 1
            G.add_edge(prev, cur)
            prev = cur
    return G


def dumbbell(a, b, l):
    G = nx.cycle_graph(a)
    G = nx.relabel_nodes(G, {i: 'a%d' % i for i in range(a)})
    G2 = nx.relabel_nodes(nx.cycle_graph(b), {i: 'b%d' % i for i in range(b)})
    G = nx.union(G, G2)
    prev = 'a0'
    for j in range(l):
        cur = 'b0' if j == l - 1 else 'p%d' % j
        G.add_edge(prev, cur)
        prev = cur
    return G


def figure8(a, b):
    G = nx.relabel_nodes(nx.cycle_graph(a), {i: 'a%d' % i for i in range(a)})
    G2 = nx.relabel_nodes(nx.cycle_graph(b), {i: ('a0' if i == 0 else 'b%d' % i) for i in range(b)})
    return nx.compose(G, G2)


cands = []
for p in range(1, 8):
    for q in range(p, 8):
        r = 9 - p - q
        if r < q:
            continue
        if p == 1 and q == 1:
            continue            # double edge
        cands.append(('theta', (p, q, r), theta_graph(p, q, r)))
for a in range(3, 7):
    for b in range(a, 7):
        l = 9 - a - b
        if l >= 1:
            cands.append(('dumbbell', (a, b, l), dumbbell(a, b, l)))
for a in range(3, 7):
    b = 9 - a
    if b >= a:
        cands.append(('figure8', (a, b), figure8(a, b)))
print('(a) n = 2 cores with 9 edges, 8 vertices:')
for kind, par, G in cands:
    ok, S = sparse_ok(G)
    print('    %-8s %-10s |V| = %d |I| = %d sparse: %s%s' % (kind, par, G.number_of_nodes(), G.number_of_edges(), ok,
                                                           '' if ok else '  (violated by %d-edge subset)' % len(S)))
# (b) placement of O in the surviving thetas
for kind, par, G in cands:
    if kind != 'theta' or not sparse_ok(G)[0]:
        continue
    res = []
    for O in G.nodes():
        H = G.copy(); H.remove_node(O)
        forest = nx.is_forest(H)
        pcomp = nx.number_connected_components(H)
        if forest and G.degree(O) == 2 + pcomp:
            nout = 3 * 7 - 2 * (9 - G.degree(O)) - G.degree(O)
            res.append((O, G.degree(O), pcomp, aut_fixing(G.copy(), O), nout))
    print('(b) theta%s: admissible O (deg O = n + p, Gamma\\O forest): %s  [(O, deg O, p, |Aut fixing O|, #outputs)]' % (par, res))
G = nx.cycle_graph(5)
res = []
for O in G.nodes():
    H = G.copy(); H.remove_node(O)
    res.append((O, G.degree(O), nx.number_connected_components(H), nx.is_forest(H), aut_fixing(G.copy(), O), 3 * 4 - 2 * (5 - 2) - 2))
print('(b) pentagon: (O, deg O, p, forest, |Aut fixing O|, #outputs) = %s' % (res[0],))
