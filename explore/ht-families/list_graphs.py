"""List the series-parallel 2-Laman graphs with 6 vertices (d = 2, n = 4) with code/enumerate_graphs.py."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import enumerate_graphs as eg
import networkx as nx

gs = eg.enumerate_tight(2, 4)
print("d=2 n=4: total", len(gs))
sp = [G for G in gs if eg.is_sp(G)]
print("SP:", len(sp))
for k, G in enumerate(sp):
    degs = sorted(dict(G.degree()).values())
    print(k, "degrees", degs, "edges", sorted(tuple(sorted(e)) for e in G.edges()), "| chain_reducible", eg.chain_reducible(G, 2))
