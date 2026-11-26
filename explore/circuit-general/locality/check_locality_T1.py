"""T = 1 checks with even loop number.
(a) form level (R, P) on census graphs (H,T) = (1,1) n = 2, 4 and (2,1) n = 2, 4 (random S)
(b) chart level, k = 2 tight bridge at (1,1) n = 4:  D_Gamma J_Phi = kappa sum_j rho_S^{(j)} D'^{[w_e*^j]}
    with rho_S^{(j)} computed on S alone (Dirichlet GFF, harmonic part omega a) and D'^{[w^j]} on Gamma' alone."""
import random
from fractions import Fraction as Fr

import locality as L
from htk.graph import Graph
from htk import catalog
from check_locality import rand_lam

rng = random.Random(7)
print('== (a) form level, T = 1, even n')
stats = {}
for (H, n, ngr, nsamp) in ((1, 2, 99, 3), (2, 2, 99, 3), (1, 4, 25, 2), (2, 4, 6, 1)):
    if H + 1 <= 2 or n <= 3:
        gl = catalog.tight_graphs(H + 1, n)
    else:
        import sys, os
        sys.path.insert(0, os.path.join(L.HERE, '../atlas'))
        from graphs_enum import admissible
        from htk.graph import from_networkx
        gl = [from_networkx(g)[0] for g in admissible(H + 1, n)]
    rng.shuffle(gl)
    for G in gl[:ngr]:
        for _ in range(nsamp):
            cand = list(range(G.nv - 1))
            W = set(rng.sample(cand, rng.randint(1, min(3, len(cand)))))
            S = [e for e, (a, b) in enumerate(G.edges) if a in W or b in W]
            if len(S) == G.ne:
                continue
            try:
                sp = L.Split(G, S)
            except ValueError:
                continue
            if sp.k < 2:
                continue
            tight = (len(S) - 1 == sp.nS * (H + 1)) if sp.k == 2 else None
            for mode in ['R'] + (['P'] if sp.k == 2 else []):
                t = L.rand_times(G.ne, rng)
                lam = rand_lam(G, H, rng)
                Db, c, U = L.factor_T1(sp, t, lam, mode)
                r = L.same_density(Db, c, U)
                ok = 'zero' if Db.coeff == 0 and c == 0 else {1: '+1', -1: '-1', None: 'FAIL'}[r]
                key = ((H, 1), n, 'k=%d' % sp.k, mode, 'tight' if tight else ('nontight' if tight is False else '-'))
                stats.setdefault(key, {}).setdefault(ok, 0)
                stats[key][ok] += 1
for key in sorted(stats):
    print('  (H,T)=%s n=%d %s mode %s %-8s : %s' % (key[0], key[1], key[2], key[3], key[4], stats[key]), flush=True)

print('== (b) chart level, (1,1) n = 4: Wheatstone bridge S (p=0, q=1, internal 2, 3) + O = {0-4, 4-5, 5-1, 4-1}')
# vertices: p=0, q=1, a=2, b=3, u=4, w=5 (root)
edges = [(0, 2), (0, 3), (2, 1), (3, 1), (2, 3), (0, 4), (4, 5), (5, 1), (4, 1)]
G = Graph(6, edges)
S = [0, 1, 2, 3, 4]
tt = L.TwoTerminal(G, S)
for H in (1,):
    vals = set()
    for _ in range(4):
        tS = L.rand_times(5, rng)
        tO = L.rand_times(4, rng)
        lam = rand_lam(G, H, rng)
        DG, J, rho, US, Dps, Up = L.check_point_T1(tt, tS, tO, lam)
        assert DG.radicand == US * Up, (DG.radicand, US, Up)
        rhs = sum((r * d for r, d in zip(rho, Dps)), Fr(0))
        vals.add(DG.coeff * J / rhs if rhs else ('0=0' if DG.coeff == 0 else 'rhs=0'))
        print('   rho_S^(j) =', [str(x) for x in rho], ' D\'^[w^j] =', [str(x) for x in Dps])
    print('  H=%d: D_G J_Phi / sum_j rho_j D\'_j = %s;  radicand U_Gamma = U_S U\'' % (H, sorted(map(str, vals))))
