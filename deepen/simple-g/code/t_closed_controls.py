"""Controls for closedness.dF_exact: beta (exact, so closed) must give 0; the invariant cochain
K = sum_sigma sgn [p_s0 p_s1][p_s0 p_s2] kappa(X_s0,X_s1) kappa(X_s2,X_s3) is a test input (closedness unknown a priori)."""
import sys, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np
import algebras as al, core, closedness as cl
from algebras import perm_sign
for key in sys.argv[1:]:
    L = al.build(key); T = core.Tensors(L); d = L.dim
    TT = core.tree_full(L)
    Bmu = np.zeros((6,) + (d,) * 4, dtype=np.int64)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = core.quad(core.br(k, l), core.br(k, l))
        V = np.transpose(TT, list(np.argsort([i, j, k, l])))
        for mu in range(6):
            if c[mu]:
                Bmu[mu] += (-1) ** (i + j) * int(c[mu]) * V
    kk = np.einsum('ab,cd->abcd', L.kappa, L.kappa)
    Kmu = np.zeros_like(Bmu)
    for sig in itertools.permutations(range(4)):
        c = core.quad(core.br(sig[0], sig[1]), core.br(sig[0], sig[2]))
        V = np.transpose(kk, list(np.argsort(sig)))
        for mu in range(6):
            if c[mu]:
                Kmu[mu] += perm_sign(sig) * int(c[mu]) * V
    print("%-12s beta: d beta = 0 %s ; test cochain K: nonzero %s, d K = 0 %s"
          % (L.name, not cl.dF_exact(L, T, Bmu).any(), Kmu.any(), not cl.dF_exact(L, T, Kmu).any()))
    # sensitivity control: F_133 coefficient tensors with the polynomial labels mu cyclically shifted
    Tmu, W = cl.T_mu_full(T, '133')
    bad = np.roll(Tmu, 1, axis=0)
    print("%-12s control: d(F_133 with shifted polynomial labels) = 0: %s ; d F_133 = 0: %s"
          % (L.name, not cl.dF_exact(L, T, bad).any(), not cl.dF_exact(L, T, Tmu).any()))
