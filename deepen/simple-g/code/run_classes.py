"""Classes of F_133, F_223, omega_133 = -3F_133, omega_223 = 2F_223, F_133 + F_223 and
O_2 ~ omega_133 + omega_223 (relative sign of claim c3) / omega_133 - omega_223 (opposite sign)
in the G x SL2-invariant part of H^4 at weight (2,2): B^4_inv = C beta (see REPORT), so a closed
invariant cochain X is exact iff X in C beta.  All ranks exact over Q (Python integers)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-g/code')
import numpy as np, sympy as sp
import algebras as al, core


def killing_ratio(L):
    adK = np.einsum('acd,bdc->ab', L.fu, L.fu)        # tr(ad a ad b)
    k = L.kappa
    i, j = np.nonzero(k)
    r = sp.Rational(int(adK[i[0], j[0]]), int(k[i[0], j[0]]))
    ok = all(sp.Rational(int(adK[a, b])) == r * int(k[a, b]) for a in range(L.dim) for b in range(L.dim))
    return r, ok


def run(key):
    t0 = time.time()
    L = al.build(key)
    r, okr = killing_ratio(L)
    T = core.Tensors(L)
    th = L.theta
    F133 = core.F_eval(T, '133', th)
    F223 = core.F_eval(T, '223', th)
    be = core.beta_eval(L, th)
    t1 = time.time()
    out = ["== %s  dim %d  form: %s ; Killing = %s * form (checked on all entries: %s) ; s = %d ; theta = basis %d"
           % (L.name, L.dim, L.form, r, okr, L.s, th)]
    if L.wt is not None:
        d = L.dim
        W = np.array(L.wt)
        tot = W[:, None, None, :] + W[None, :, None, :] + W[None, None, :, :] + W[th][None, None, None, :]
        bad = np.any(tot != 0, axis=3)
        wok = all(not X[:, bad].any() for X in (F133, F223, be))
        out.append("   h-weight selection rule (entries with wt(b)+wt(c)+wt(d) != -wt(theta) vanish): %s ; "
                   "compatible triples %d of %d" % (wok, int((~bad).sum()), d ** 3))
    cols = {'beta': be, 'F133': F133, 'F223': F223,
            'F133+F223': F133 + F223,
            'omega133+omega223': -3 * F133 + 2 * F223,
            'omega133-omega223': -3 * F133 - 2 * F223}
    nz = {k: bool(v.any()) for k, v in cols.items()}
    out.append("   nonzero: " + ", ".join("%s %s" % (k, nz[k]) for k in ('beta', 'F133', 'F223')))
    rows3 = core.to_rows([be, F133, F223])
    out.append("   rows used (nonzero) %d ; rank[beta] = %d, rank[beta,F133] = %d, rank[beta,F223] = %d, rank[beta,F133,F223] = %d"
               % (len(rows3), core.rank_exact([[x[0]] for x in rows3]), core.rank_exact([[x[0], x[1]] for x in rows3]),
                  core.rank_exact([[x[0], x[2]] for x in rows3]), core.rank_exact(rows3)))
    s7 = T.s7
    res = {}
    for name in ('F133', 'F223', 'F133+F223', 'omega133+omega223', 'omega133-omega223'):
        rows = core.to_rows([be, cols[name]])
        ok, sol = core.in_span(rows, 1)
        res[name] = ok
        if ok:
            c = sol[0] / s7
            out.append("   %-18s exact: %s = c beta, c = %s (form used) ; c_Killing = c / r^2 = %s" % (name, name, c, c / r ** 2))
        else:
            out.append("   %-18s NOT exact (rank[beta, X] = 2 > rank[beta] = 1)" % name)
    out.append("   time %.1fs (cochains %.1fs)" % (time.time() - t0, t1 - t0))
    return "\n".join(out), res


if __name__ == "__main__":
    for key in sys.argv[1:]:
        txt, _ = run(key)
        print(txt, flush=True)
