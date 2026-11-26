"""tight two-terminal S without internal sources, H >= 2, n_S >= 1: D_S = 0 (currents of S all
proportional to the terminal current J: psi^i_r(S) = J^i alpha_r, wedge over i vanishes), hence D_Gamma = 0
pointwise; with internal sources D_Gamma != 0."""
import random
import locality as L
from htk import brute
from check_locality import k4sub, relabel_root, rand_lam
rng = random.Random(31)
for H, chains in ((3, [3, 2, 1, 1, 2, 1]), (3, [3, 1, 2, 1, 1, 2]), (2, [2, 1, 1, 1, 1, 1])):
    nv, ed, ch = k4sub(chains, H)
    G, new = relabel_root(nv, ed, 4)
    S = [e for c in ch[1:] for e in c]
    sp = L.Split(G, S)
    out = []
    for internal_sources in (False, True):
        vals = []
        for _ in range(3):
            lam = rand_lam(G, H, rng)
            if not internal_sources:
                for v in sp.S0:
                    lam[v] = [0] * H
            vals.append(brute.density(G, L.rand_times(G.ne, rng), lam).coeff != 0)
        out.append((internal_sources, vals))
    print('H=%d K4%s: S tight=%s; (internal sources, D_Gamma != 0 at 3 points): %s' % (
        H, chains, len(S) - 1 == sp.nS * H, out))
