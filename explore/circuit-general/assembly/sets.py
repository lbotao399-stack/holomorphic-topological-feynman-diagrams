"""Graph sets and sources of atlas/e3_open_kernels.py (same seeds), as htk edge lists."""
import os
import random
import sys
import warnings

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS = os.path.abspath(os.path.join(HERE, '..', 'atlas'))
sys.path.insert(0, ATLAS)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..', 'code')))

import networkx as nx            # noqa: E402
import htk                       # noqa: E402
from htk.graph import from_networkx   # noqa: E402
from atlas import set_twistor, set_bulk3, set_cone   # noqa: E402


def e3_sources(Gh, lab, S0, root, seed, H):
    """identical to atlas/e3_open_kernels.py sources()"""
    rng = random.Random(seed)
    lam = {lab[v]: [rng.randint(-3, 3) for _ in range(H)] for v in S0 if v != root}
    full = [lam.get(v, [0] * H) for v in range(Gh.nv - 1)]
    full.append([-sum(r[i] for r in full) for i in range(H)])
    return lam, full


def items(name):
    """list of (idx, G networkx, root O, source vertices S0)"""
    if name == 'bulk3':
        out = []
        for it in set_bulk3(3, 4):
            G = it[0]
            S0 = [v for v in G if G.degree(v) == 2]
            out.append((G, S0[-1], S0))
        return list(enumerate(out))
    if name == 'twistor4':
        s = set_twistor(4)
        return list(enumerate([it for k in ('tree', 'unicyclic') for it in s[k]]))
    if name == 'twistor3':
        s = set_twistor(3)
        return list(enumerate([it for k in ('tree', 'unicyclic') for it in s[k]]))
    if name in ('bulk3_5', 'bulk3_4', 'bulk3_6'):
        D = int(name[-1])
        out = []
        for it in set_bulk3(D, 3):
            G = it[0]
            S0 = [v for v in G if G.degree(v) == 2]
            out.append((G, S0[-1], S0))
        return list(enumerate(out))
    if name in ('cone5', 'cone4'):
        D = int(name[-1])
        return list(enumerate(set_cone(D, 3)))
    raise KeyError(name)


def instance(name, idx, G, O, S0, H):
    """(htk Graph, full source list) with the e3 seed 1000*idx + 7"""
    Gh, lab = from_networkx(G, root=O)
    lam, full = e3_sources(Gh, lab, S0, O, 1000 * idx + 7, H)
    return Gh, lab, lam, full


def is_sp(G):
    return htk.catalog.is_sp(nx.Graph(G))
