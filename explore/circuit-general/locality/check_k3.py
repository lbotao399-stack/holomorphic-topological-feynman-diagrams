"""Task 3: three-terminal elements (stars and composite stars) at D = 2: (H,T) = (2,0) and (1,1).
For every census graph (D = 2, n = 4) and every interior set W of the listed types, S = edges at W:
  single star        W = {v}, deg v = 3
  star + leg         W = {v, w}, deg v = 3, w a bivalent neighbour of v (k = 3, m = 2)
  paired stars       W = {v1, v2}, deg 3, same three neighbours (K_{2,3}, n_S = 2)
modes: R (response differentials dZ, dh; T = 1 with the Dirichlet source xi = L_00 phi^D),
       P (dh terms and xi terms dropped).  Sources: interior sourced / unsourced.
Expected [drv]: R exact always; P exact at T = 0 iff dh.lambda_S0 = 0 on the relevant degree;
at T = 1, P fails also for unsourced interior (contact term d.Z.dh.xi)."""
import random
from itertools import combinations

import locality as L
from htk import catalog
from check_locality import rand_lam


def pieces(G):
    out = []
    deg = [G.degree(v) for v in range(G.nv)]
    nb = {v: set() for v in range(G.nv)}
    for a, b in G.edges:
        nb[a].add(b)
        nb[b].add(a)
    for v in range(G.nv - 1):
        if deg[v] == 3:
            out.append(('single star', {v}))
            for w in nb[v]:
                if deg[w] == 2 and w != G.nv - 1:
                    out.append(('star + subdivided leg', {v, w}))
    for v1, v2 in combinations(range(G.nv - 1), 2):
        if deg[v1] == 3 and deg[v2] == 3 and nb[v1] == nb[v2]:
            out.append(('paired stars K23', {v1, v2}))
    return out


rng = random.Random(23)
stats = {}
for (H, T) in ((2, 0), (1, 1)):
    for G in catalog.tight_graphs(2, 4):
        for name, W in pieces(G):
            S = [e for e, (a, b) in enumerate(G.edges) if a in W or b in W]
            try:
                sp = L.Split(G, S)
            except ValueError:
                continue
            if sp.k != 3:
                continue
            for sourced in (True, False):
                for mode in ('R', 'P'):
                    t = L.rand_times(G.ne, rng)
                    lam = rand_lam(G, H, rng)
                    if not sourced:
                        for v in W:
                            lam[v] = [0] * H
                    if T == 0:
                        Db, Df, _ = L.factor_T0(sp, t, lam, mode)
                        ok = 'zero' if Db == 0 and Df == 0 else ('+-1' if Db in (Df, -Df) else 'FAIL')
                    else:
                        Db, c, U = L.factor_T1(sp, t, lam, mode)
                        r = L.same_density(Db, c, U)
                        ok = 'zero' if Db.coeff == 0 and c == 0 else ('+-1' if r else 'FAIL')
                    key = ((H, T), name, 'interior sourced' if sourced else 'interior unsourced', mode)
                    stats.setdefault(key, {}).setdefault(ok, 0)
                    stats[key][ok] += 1
for key in sorted(stats):
    print('  (H,T)=%s %-24s %-20s mode %s : %s' % (key[0], key[1], key[2], key[3], stats[key]))
