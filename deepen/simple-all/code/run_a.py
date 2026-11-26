"""Part (a): classes of F_133, F_223, O_2 in the G x SL2-invariant part of H^4_(2,2), full coefficient tensors.

For each algebra (lie.build key):
  W_133, W_223 (t4.Engine.weights, integer s^7 W), F_G^mu = alt([p0p1][p2p3] W_G) and beta^mu (mu over QB), all
  entries (no fixed slot).  g-invariance of W_133, W_223 under the generators E_{+-alpha_i} checked on all entries.
  X in C beta  <=>  X[i] beta[i0] - X[i0] beta[i] = 0 for all entries i (i0 an entry with beta[i0] != 0):
  decided exactly (Python integers); otherwise a nonzero 2x2 minor is printed.
  dim (Lambda^3 g*)^g, (S^3 g*)^g, (g^{*(x)3})^g: inv3.invariant_dim (nullity mod two primes).
"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import scipy.sparse as sps
from fractions import Fraction as Fr
import lie, t4, inv3


def act(L, T, y):
    """(y.T) = - sum_i T(.., [y, .]_i, ..) on all entries; returns number of nonzero entries."""
    d = L.dim
    M = L.ad(y).tocsr().astype(np.int64)            # M[c, b]: coefficient of e_c in [y, e_b]
    keys, vals = [], []
    for i in range(4):
        others = [k for k in range(4) if k != i]
        r = (T.idx[:, others[0]] * d + T.idx[:, others[1]]) * d + T.idx[:, others[2]]
        Ti = sps.csr_matrix((T.val, (r, T.idx[:, i])), shape=(d ** 3, d))
        R = (Ti @ M).tocoo()
        idx = np.zeros((len(R.data), 4), dtype=np.int64)
        rr = R.row
        idx[:, others[2]] = rr % d; rr = rr // d
        idx[:, others[1]] = rr % d; rr = rr // d
        idx[:, others[0]] = rr
        idx[:, i] = R.col
        keys.append(((idx[:, 0] * d + idx[:, 1]) * d + idx[:, 2]) * d + idx[:, 3])
        vals.append(-R.data)
    keys, vals = np.concatenate(keys), np.concatenate(vals)
    o = np.argsort(keys)
    keys, vals = keys[o], vals[o]
    u, st = np.unique(keys, return_index=True)
    s = np.add.reduceat(vals, st) if len(vals) else vals
    return int(np.count_nonzero(s))


def in_C_beta(Xv, Bv):
    """Xv, Bv: aligned lists (6 mu) of object arrays.  Returns (True, c) with X = c beta, or (False, minor)."""
    for mu in range(6):
        nz = np.nonzero(Bv[mu])[0]
        if len(nz):
            i0 = (mu, nz[0]); break
    b0, x0 = int(Bv[i0[0]][i0[1]]), int(Xv[i0[0]][i0[1]])
    for mu in range(6):
        D = Xv[mu] * b0 - Bv[mu] * x0
        nz = np.nonzero(D)[0]
        if len(nz):
            j = nz[0]
            return False, (b0, int(Xv[mu][j]), x0, int(Bv[mu][j]), int(D[j]))
    return True, Fr(x0, b0)


def run(key, invk='LST'):
    t0 = time.time()
    L = lie.build(key)
    r, okr = L.killing_ratio()
    out = ["== %s  dim %d ; form: %s ; Killing = %s * form (all entries: %s) ; s = %d"
           % (L.name, L.dim, L.form, r, okr, L.s)]
    gens = inv3.generators(L)
    gd = inv3.generated_dim(L, gens)
    line = "   generators E_{+-alpha_i}: %d, generated dim %d" % (len(gens), gd)
    for k in invk:
        nul, nu, ne, ng = inv3.invariant_dim(L, k, gens=gens)
        line += " ; dim %s nullity mod p = %s (unknowns %d, eqs %d)" % ({'L': '(Lambda^3)^g', 'S': '(S^3)^g', 'T': '(T^3)^g'}[k], nul, nu, ne)
    out.append(line)
    t1 = time.time()
    E = t4.Engine(L)
    W133, W223 = E.weights()
    inv_ok = all(act(L, W, y) == 0 for W in (W133, W223) for y in gens)
    out.append("   W_133, W_223: nnz %d, %d ; max |s^7 W| %d ; y.W = 0 for all generators, all entries: %s ; Om nnz %d"
               % (W133.nnz(), W223.nnz(), max(int(abs(W133.val).max()), int(abs(W223.val).max())), inv_ok, E.Om.nnz()))
    F133, F223 = t4.alt_bb(W133), t4.alt_bb(W223)
    B = t4.beta_t(E.ff)
    keys, vals = t4.merge(F133 + F223 + B)
    F133v, F223v, Bv = vals[0:6], vals[6:12], vals[12:18]
    s7 = L.s ** 7
    cols = {'F133': F133v, 'F223': F223v, 'F133+F223': [a + b for a, b in zip(F133v, F223v)],
            'O2: -3F133+2F223': [-3 * a + 2 * b for a, b in zip(F133v, F223v)],
            'O2 flip: -3F133-2F223': [-3 * a - 2 * b for a, b in zip(F133v, F223v)]}
    out.append("   entries (union over mu of nonzero keys) %d ; nonzero mu components: F133 %s, F223 %s, beta %s"
               % (len(keys), [int(np.count_nonzero(x)) for x in F133v], [int(np.count_nonzero(x)) for x in F223v],
                  [int(np.count_nonzero(x)) for x in Bv]))
    res = {}
    for name, Xv in cols.items():
        ok, info = in_C_beta(Xv, Bv)
        res[name] = (ok, info)
        if ok:
            c = info / s7
            out.append("   %-22s = c beta on all entries: c = %s (form used), c_Killing = c / r^2 = %s" % (name, c, c / r ** 2))
        else:
            b0, xj, x0, bj, D = info
            out.append("   %-22s NOT in C beta: minor beta[i0] X[j] - X[i0] beta[j] = %d*%d - %d*%d = %d (s^7 units)"
                       % (name, b0, xj, x0, bj, D))
    out.append("   time %.1fs (invariants %.1fs)" % (time.time() - t0, t1 - t0))
    return "\n".join(out), res


if __name__ == "__main__":
    invk = 'LST'
    args = sys.argv[1:]
    if args and args[0].startswith('--inv='):
        invk = args[0][6:]; args = args[1:]
    for key in args:
        txt, _ = run(key, invk)
        print(txt, flush=True)
