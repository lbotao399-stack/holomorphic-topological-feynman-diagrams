"""Is d nonzero on the G x SL2-invariant alternating 4-cochains?  Test inputs C = alt(P W) with P a product of two
brackets [p_a p_b][p_c p_d] (all choices) and W in {W_133, W_223, kappa kappa, tree}; exact dF_exact."""
import sys, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np
import algebras as al, core, closedness as cl
from algebras import perm_sign
for key in sys.argv[1:]:
    L = al.build(key); T = core.Tensors(L); d = L.dim
    tens = {'W133': T.W('133', [None] * 4), 'W223': T.W('223', [None] * 4),
            'kk': np.einsum('ab,cd->abcd', L.kappa, L.kappa), 'tree': core.tree_full(L)}
    res = []
    for name, W in tens.items():
        for (a, b, c, e) in itertools.product(range(4), repeat=4):
            if a >= b or c >= e:
                continue
            Tmu = np.zeros((6,) + (d,) * 4, dtype=np.int64)
            for sig in itertools.permutations(range(4)):
                q = core.quad(core.br(sig[a], sig[b]), core.br(sig[c], sig[e]))
                V = np.transpose(W, list(np.argsort(sig)))
                for mu in range(6):
                    if q[mu]:
                        Tmu[mu] += perm_sign(sig) * int(q[mu]) * V
            if not Tmu.any():
                continue
            nz = cl.dF_exact(L, T, Tmu).any()
            res.append((name, (a, b, c, e), nz))
    print(L.name, "nonzero alternated inputs:", len(res), "; inputs with d C != 0:", [r[:2] for r in res if r[2]])
