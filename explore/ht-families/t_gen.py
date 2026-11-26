"""check density_batch_general == density_batch (n=4) and == exact brute (n=4, n=6) at points"""
import sys, random
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
import htk
import enumerate_graphs as eg
import hf_brute as hb, hf_numeric as hn
rng = random.Random(3)
for (D, n, k) in [(2, 4, 3), (2, 6, 6), (2, 6, 0)]:
    gs = [Gx for Gx in eg.enumerate_tight(D, n) if eg.is_sp(Gx)]
    G, _ = htk.from_networkx(gs[k], None)
    Gm = hb.G(G.nv, G.edges); C = Gm.cycles()
    lam = [[Fr(rng.randint(-5, 5))] for _ in range(G.nv - 1)]
    for _ in range(2):
        t = [Fr(rng.randint(1, 20), rng.randint(1, 9)) for _ in G.edges]
        c, U = hb.density(Gm, t, lam, T=1, C=C)
        ex = float(c) / float(U) ** 0.5
        tf = np.array([[float(x) for x in t]])
        lf = np.array([[float(l[0])] for l in lam])
        v = hn.density_batch_general(G.nv, G.edges, C, tf, lf)[0]
        print("D=%d n=%d graph %d: exact %.12g  float-general %.12g  rel %.1e" % (D, n, k, ex, v, abs(v - ex) / abs(ex)))
