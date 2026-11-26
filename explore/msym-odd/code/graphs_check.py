"""Two-loop graph lists for D = H + T = 2, 3, 4 (htk.catalog, counting of CONVENTIONS Section 7) and
the odd-coordinate rule (S2 eq-theta-rule): O at a vertex with Gamma \\ O a forest, deg O = n + p."""
import networkx as nx
import htk_path  # noqa
from htk import catalog
for D in (2, 3, 4):
    gs = catalog.tight_graphs(D, 2)
    out = []
    for G in gs:
        g = G.to_networkx()
        degs = sorted(d for _, d in g.degree())
        allowed = []
        for v in g.nodes():
            h = g.copy(); h.remove_node(v)
            p = nx.number_connected_components(h)
            forest = nx.is_forest(nx.Graph(h)) and h.number_of_edges() == nx.Graph(h).number_of_edges()
            if forest and g.degree(v) == 2 + p:
                allowed.append((v, g.degree(v)))
        out.append((G.ne, G.nv, degs, len(allowed), sorted(set(d for _, d in allowed))))
    print("D=%d: %d tight two-loop graphs; per graph (|I|, |V|, degrees, #admissible O positions, deg O): %s" % (D, len(gs), out))
    print("      theta_graphs:", [t[0] for t in catalog.theta_graphs(D)])
