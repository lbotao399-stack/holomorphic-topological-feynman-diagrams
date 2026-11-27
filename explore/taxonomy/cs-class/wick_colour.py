"""(4,1) and (5,0) thetas: standard Wick signs sigma_Gamma (flagship rule), invariance of sigma K C under
vertex relabelling (root anywhere), edge reversal and edge order (with a control), and the gl_N ribbon
expansion of C_std against explicit sl_3, sl_4 contractions (exact)."""
import random
from fractions import Fraction as Fr
from itertools import permutations
import numpy as np
from common import Log, parity, perm_from_lists
import colour as col
import twoloop
import lie

log = Log(__file__)
SHAPES = [(1, 5, 5), (2, 4, 5), (3, 3, 5), (3, 4, 4)]
rng = random.Random(7)

log('## standard Wick signs (local order L\', target (e,0),(e,1) by edge, then externals)')
SIG = {}
for shape in SHAPES:
    SIG[shape] = col.standard_sigma(shape)
    log(shape, 'sigma =', SIG[shape])


def physical_roles(shape):
    """for each old vertex: list of role keys in standard local order."""
    nv, edges = col.theta_edges(shape)
    A, B, C = col.theta_paths(shape)
    roles = {}
    # path edge lists
    paths, e = [], 0
    for n_ in shape:
        paths.append(list(range(e, e + n_)))
        e += n_
    roles[0] = [('E', paths[k][0]) for k in range(3)]
    roles[nv - 1] = [('E', paths[k][-1]) for k in range(3)]
    for k, ints in enumerate((A, B, C)):
        for j, v in enumerate(ints):
            roles[v] = [('x', v), ('E', paths[k][j]), ('E', paths[k][j + 1])]
    return nv, edges, roles


def check_invariance(shape, H, T, ntrial=6, control=False):
    nv, edges, roles = physical_roles(shape)
    ext_old = list(range(1, nv - 1))
    lam_phys = {v: [Fr(rng.randint(-5, 5)) for _ in range(H)] for v in range(nv)}
    tot = [-sum(lam_phys[v][i] for v in range(nv - 1)) for i in range(H)]
    lam_phys[nv - 1] = tot

    def kernel(nv_, edges_, lam_by_new):
        th = twoloop.Theta(nv_, edges_, root=nv_ - 1)
        lamfull = [lam_by_new[v] for v in range(nv_)]
        return twoloop.kernel_htk_closed(th, lamfull, H, T)
    K_old = kernel(nv, edges, lam_phys)
    sig_old = SIG[shape] if (H, T) == (4, 1) or True else None
    sig_old = col.wick_sign(nv, edges, ext_old, col.local_orders_Lprime(nv, edges, ext_old))
    ok_all = True
    for trial in range(ntrial):
        # new label of old vertex u: perm
        perm = list(range(nv))
        rng.shuffle(perm)
        new_of = {u: perm[u] for u in range(nv)}
        old_of = {perm[u]: u for u in range(nv)}
        # edges: random order, random reversals
        order = list(range(len(edges)))
        rng.shuffle(order)
        rev = [rng.random() < 0.5 for _ in edges]
        new_edges, new_id_of_old, rev_of_new = [], {}, {}
        for k, eo in enumerate(order):
            a, b = edges[eo]
            if rev[eo]:
                a, b = b, a
            new_edges.append((new_of[a], new_of[b]))
            new_id_of_old[eo] = k
            rev_of_new[k] = rev[eo]
        ext_new = sorted(new_of[u] for u in ext_old)
        lo = col.local_orders_Lprime(nv, new_edges, ext_new)
        sig_new = col.wick_sign(nv, new_edges, ext_new, lo)
        # epsilon: local order L' (new) vs physical role order, per vertex
        eps = 1
        for vnew in range(nv):
            u = old_of[vnew]
            keys = []
            for item in lo[vnew]:
                if item[0] == 'x':
                    keys.append(('x', u))
                else:
                    e_new, end = item
                    keys.append(('E', order[e_new]))
            eps *= perm_from_lists(keys, roles[u])
        if control and trial == 0:
            eps = -eps
        # external reordering: new vertex order (as old labels) -> old order
        ext_new_as_old = [old_of[v] for v in ext_new]
        s_ext = perm_from_lists(ext_new_as_old, sorted(ext_new_as_old))
        lam_new = {new_of[u]: lam_phys[u] for u in range(nv)}
        # htk grounding: root = last label; sources sum to zero already
        K_new = kernel(nv, new_edges, lam_new)
        lhs = sig_new * eps * s_ext * K_new
        rhs = sig_old * K_old
        ok = lhs == rhs
        ok_all &= ok
    return ok_all, K_old


log('## invariance of sigma K C (6 random relabellings each, sources random integers)')
for HT in ((4, 1), (5, 0)):
    for shape in SHAPES:
        ok, K0 = check_invariance(shape, *HT)
        okc, _ = check_invariance(shape, *HT, ntrial=1, control=True)
        log(HT, shape, 'invariant:', ok, ' control (one local-order sign flipped) invariant:', okc)
        assert ok and not okc

log('## ribbon expansion of C_std (gl_N, trace basis) vs explicit sl_3, sl_4 contraction (exact)')
for shape in SHAPES:
    rib = col.ribbon_theta_gl(shape)
    npow = sorted({n for v in rib.values() for n in v})
    nfpf = sum(1 for k in rib if not any(len(w) == 1 for w in k))
    log(shape, 'trace structures:', len(rib), ' fixed-point-free:', nfpf, ' N powers:', npow)
for N in (3, 4):
    L = lie.sl_tr(N)
    d = L.dim
    fl = np.zeros((d, d, d), dtype=object)
    a, b, c, v = L.f_low
    for x, y, z, w in zip(a, b, c, v):
        fl[x, y, z] = Fr(int(w))
    G = L.G.toarray()
    kinv = np.array([[Fr(int(G[i, j]), L.s) for j in range(d)] for i in range(d)], dtype=object)
    adm = []
    for k in range(d):
        M = L.ad(k).toarray()
        adm.append(np.array([[Fr(int(M[i, j])) for j in range(d)] for i in range(d)], dtype=object))

    def ad_of(cv):
        M = np.zeros((d, d), dtype=object)
        for k in range(d):
            if cv[k]:
                M = M + cv[k] * adm[k]
        return M
    basis_mats = []
    for m in L.basis:
        Mx = [[Fr(0)] * N for _ in range(N)]
        for (i, j), w in m.items():
            Mx[i][j] = Fr(w)
        basis_mats.append(np.array(Mx, dtype=object))
    for shape in SHAPES:
        cd = col.colour_perm_dict(shape, N)
        Xs = {vv: [Fr(rng.randint(-3, 3)) for _ in range(d)] for vv in range(1, 9)}
        mats = {vv: sum((Xs[vv][k] * basis_mats[k] for k in range(d)), np.zeros((N, N), dtype=object)) for vv in Xs}
        explicit = col.eval_theta_dense(shape, fl, kinv, ad_of, Xs)
        val = Fr(0)
        for perm, cf in cd.items():
            # successor perm on positions 0..7 (letter = vertex pos+1)
            seen, term = set(), Fr(1)
            for s in range(8):
                if s in seen:
                    continue
                cyc, x = [], s
                while x not in seen:
                    seen.add(x); cyc.append(x); x = perm[x]
                M = np.eye(N, dtype=object) * Fr(1)
                for x in cyc:
                    M = M.dot(mats[x + 1])
                term *= sum(M[i, i] for i in range(N))
            val += cf * term
        log('sl_%d' % N, shape, 'explicit =', explicit, ' trace basis =', val, ' equal:', explicit == val)
        assert explicit == val
log.close()
