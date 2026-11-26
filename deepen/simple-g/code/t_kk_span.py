"""(g^{*(x)4})^g versus span{kappa_12 kappa_34, kappa_13 kappa_24, kappa_14 kappa_23}: exact rank of the three
products (all entries) against dim (g^{(x)4})^g of the Weyl formula; for the algebras where they agree, the rank of
the alternated cochains alt([p_a p_b][p_c p_d] kappa_01 kappa_23) (all a<b, c<d) together with beta
(evaluation at b_theta, exact) = dim of the G x SL2-invariant alternating 4-cochains of momentum degree 4."""
import sys, itertools
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np
import algebras as al, core, invariants as inv
from algebras import perm_sign
for key in sys.argv[1:]:
    L = al.build(key); d = L.dim; th = L.theta
    k = L.kappa
    kk = [np.einsum('ab,cd->abcd', k, k), np.einsum('ac,bd->abcd', k, k), np.einsum('ad,bc->abcd', k, k)]
    r_kk = core.rank_exact([[int(x) for x in row] for row in np.stack([t.reshape(-1) for t in kk], 1) if any(row)])
    dims, _ = inv.invariant_dims(L, [(4, 'T')])
    line = "%-14s dim (g^(x)4)^g = %d ; rank of the 3 kappa-kappa products = %d" % (L.name, dims[(4, 'T')], r_kk)
    if r_kk == dims[(4, 'T')]:
        cols = [core.beta_eval(L, th)]
        for (a, b, c, e) in itertools.product(range(4), repeat=4):
            if a >= b or c >= e:
                continue
            out = np.zeros((6, d, d, d), dtype=np.int64)
            for sig in itertools.permutations(range(4)):
                q = core.quad(core.br(sig[a], sig[b]), core.br(sig[c], sig[e]))
                V = np.transpose(kk[0], list(np.argsort(sig)))[th]
                for mu in range(6):
                    if q[mu]:
                        out[mu] += perm_sign(sig) * int(q[mu]) * V
            cols.append(out)
        rows = core.to_rows(cols)
        line += " ; rank[beta, alt([p_a p_b][p_c p_d] kappa_01 kappa_23) (36 inputs)] = %d, without beta = %d" % (
            core.rank_exact(rows), core.rank_exact([r[1:] for r in rows]))
    print(line, flush=True)
