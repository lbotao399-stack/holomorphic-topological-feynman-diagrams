"""Exact kernels for atlas rows left open, with the cores evaluator (explore/circuit-general/cores).
modes:
  twistor3          twistor BF Q_3, D=3: the 14 non-SP graphs (fastker), plus SP cross-checks
                    (fastker vs htk.kernel_exact on the same graph and sources)
  tadpole50         (H,T)=(1,1) tadpole(5,0) of ht-d2/forests.py, T=1 (cores.kernel)
  bulk3 K S         bulk3(D=3, n=4) (hCS / Kodaira-Spencer on C^3), shard S of K (index mod K)
  twistor4 K S      twistor BF Q_4 non-SP graphs, shard S of K
Sources: random integers (seeded per graph) at the source vertices, root = last vertex.
usage: python3 e3_open_kernels.py MODE [K S]"""
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../../../code")))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "../cores")))

import networkx as nx
import htk
from htk.graph import from_networkx
from atlas import set_twistor, set_bulk3
import fastker
import cores

H3 = 3


def sources(Gh, lab, S0, root, seed, H):
    rng = random.Random(seed)
    lam = {lab[v]: [rng.randint(-3, 3) for _ in range(H)] for v in S0 if v != root}
    full = [lam.get(v, [0] * H) for v in range(Gh.nv - 1)]
    full.append([-sum(r[i] for r in full) for i in range(H)])
    return lam, full


def fmt(val):
    return {k: str(v) for k, v in val.items()} if isinstance(val, dict) else str(val)


def run_items(items, label, only_nonsp, shard=None):
    for idx, (G, O, S0) in enumerate(items):
        if shard and idx % shard[0] != shard[1]:
            continue
        sp = htk.catalog.is_sp(nx.Graph(G))
        if only_nonsp and sp:
            continue
        Gh, lab = from_networkx(G, root=O)
        lam, full = sources(Gh, lab, S0, O, 1000 * idx + 7, H3)
        t0 = time.time()
        try:
            val = fmt(fastker.kernel(Gh.nv, list(Gh.edges), full, verbose=False))
        except Exception as e:
            val = "FAILED %r" % (e,)
        msg = "%s #%d SP=%s fastker K = %s (%.0fs)" % (label, idx, sp, val, time.time() - t0)
        if sp:
            msg += " | spred K = %s" % (htk.kernel_exact(Gh, lam, T=0, check=1),)
        print(msg, flush=True)


mode = sys.argv[1]
shard = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else None
if mode == "twistor3":
    s = set_twistor(3)
    items = [it for k in ("tree", "unicyclic") for it in s[k]]
    run_items(items, "twistor3", only_nonsp=False)
elif mode == "tadpole50":
    edges = sorted([(0, 1), (1, 2), (2, 3), (3, 4), (0, 4), (1, 5), (2, 5), (3, 5), (4, 5)])
    for lamq in (1, 2):
        lam = [[lamq], [0], [0], [0], [0], [-lamq]]
        t0 = time.time()
        val = cores.kernel(6, edges, lam, (), verbose=False, T=1)
        print("tadpole(5,0) (H,T)=(1,1), source %d at q=0, -%d at O: K = %s (%.0fs)" % (lamq, lamq, fmt(val), time.time() - t0), flush=True)
elif mode == "bulk3":
    items = []
    for it in set_bulk3(3, 4):
        G = it[0]
        S0 = [v for v in G if G.degree(v) == 2]
        items.append((G, S0[-1], S0))
    run_items(items, "bulk3(3,4)", only_nonsp=False, shard=shard)
elif mode == "twistor4":
    s = set_twistor(4)
    items = [it for k in ("tree", "unicyclic") for it in s[k]]
    run_items(items, "twistor4", only_nonsp=True, shard=shard)
