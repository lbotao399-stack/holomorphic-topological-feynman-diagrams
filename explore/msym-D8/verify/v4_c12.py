"""C12: hatQ1^2 + sum_G k_G {hatQ0, hatQ_G} = 0 for su(2), own operators (v4_ops.py), own theta weights
(theta_weights.pkl from v3_thetas.py), test operators not used by the explorer.
Usage: python3 v4_c12.py NAME [NAME ...]  (names below); results appended to v4_c12.log, rows to v4_c12_rows.pkl."""
import sys, time, pickle
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
import v4_ops as V
from v4_ops import Theory, letter as Lt, op, e, padd, Fr, log

V.LOG = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v4_c12.log', 'a')
THETAS = [(1, 4, 4), (2, 3, 4), (3, 3, 3)]
AUT = {(1, 4, 4): 2, (2, 3, 4): 1, (3, 3, 3): 6}
g = lambda A, k=(0, 0, 0, 0): Lt('g', A, k)
c = lambda A, k=(0, 0, 0, 0): Lt('c', A, k)
TESTS = {
    'g0 g0 c1': op(g(0), g(0), c(1)),
    'g2 g1 c1': op(g(2), g(1), c(1)),
    'g0 g1 g1': op(g(0), g(1), g(1)),
    'd2c0 g1 g2': op(c(0, e(1)), g(1), g(2)),
    'd3g2 g2 c0': op(g(2, e(2)), g(2), c(0)),
    'g1 d4g0 c2': op(g(1), g(0, e(3)), c(2)),
    'g0 g1 c2 [explorer op]': op(g(0), g(1), c(2)),
    'g0 g1 g2 [explorer op]': op(g(0), g(1), g(2)),
}
ROWS = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v4_c12_rows.pkl'


def solve(eqs, n):
    basis = []
    for r in eqs:
        r = list(r)
        for (cc, br) in basis:
            if r[cc] != 0:
                fct = r[cc] / br[cc]
                r = [x - fct * y for x, y in zip(r, br)]
        piv = next((i for i in range(n) if r[i] != 0), None)
        if piv is None:
            if r[n] != 0:
                return len(basis), None, False
            continue
        basis.append((piv, r))
    rank = len(basis)
    if rank < n:
        return rank, None, True
    sol = [None] * n
    for cc, br in sorted(basis, key=lambda x: -x[0]):
        sol[cc] = (br[n] - sum(br[j] * sol[j] for j in range(n) if j != cc and sol[j] is not None)) / br[cc]
    cons = all(sum(r[i] * sol[i] for i in range(n)) == r[n] for r in eqs)
    return rank, sol, cons


def main(names):
    tw = pickle.load(open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/theta_weights.pkl', 'rb'))
    n, f = V.su2()
    T = Theory(n, f, tw)
    try:
        rows = pickle.load(open(ROWS, 'rb'))
    except Exception:
        rows = {}
    for nm in names:
        if nm in rows:
            log('%s: rows already computed' % nm)
            continue
        p = TESTS[nm]
        t0 = time.time()
        q1 = T.Q1(p)
        lhs = T.Q1(q1)
        log('%s: |hatQ1 p| = %d, |hatQ1^2 p| = %d (%.1fs)' % (nm, len(q1), len(lhs), time.time() - t0))
        q0p = T.Q0(p)
        anti = {}
        for G in THETAS:
            ta = time.time()
            qg = T.Qtheta(G, p)
            a = padd(T.Q0(qg), T.Qtheta(G, q0p))
            anti[G] = a
            log('   theta%s: |hatQ_G p| = %d, |{hatQ0, hatQ_G} p| = %d (%.1fs)' % (G, len(qg), len(a), time.time() - ta))
        rows[nm] = (lhs, anti)
        pickle.dump(rows, open(ROWS, 'wb'))
    # equations over all computed operators
    eqs = []
    for nm, (lhs, anti) in rows.items():
        mons = set(lhs)
        for a in anti.values():
            mons |= set(a)
        sub = []
        for m in mons:
            sub.append([Fr(anti[G].get(m, 0)) for G in THETAS] + [-Fr(lhs.get(m, 0))])
        r_, s_, c_ = solve(sub, 3)
        log('   operator %s: %d equations, rank %d, consistent %s, solution %s' % (nm, len(sub), r_, c_, None if s_ is None else [str(x) for x in s_]))
        if r_ == 1:
            ref = anti[THETAS[0]]
            m0 = next(iter(ref))
            av = [Fr(anti[G].get(m0, 0)) / Fr(ref[m0]) for G in THETAS]
            bv = Fr(lhs.get(m0, 0)) / Fr(ref[m0])
            prop = all(Fr(anti[G].get(m, 0)) == av[i] * Fr(ref.get(m, 0)) for i, G in enumerate(THETAS) for m in mons) \
                and all(Fr(lhs.get(m, 0)) == bv * Fr(ref.get(m, 0)) for m in mons)
            log('      rank 1: {hatQ0,hatQ_G} = a_G X, hatQ1^2 = b X with a = %s, b = %s (exact on all monomials: %s); '
                'a.(2,-4,-2/3) + b = %s' % ([str(x) for x in av], bv, prop, sum(x * y for x, y in zip(av, (2, -4, Fr(-2, 3)))) + bv))
        eqs += sub
    rank, sol, cons = solve(eqs, 3)
    log('ALL operators %s: %d equations, rank %d, consistent %s, solution (k144, k234, k333) = %s'
        % (sorted(rows), len(eqs), rank, cons, None if sol is None else [str(x) for x in sol]))
    if sol is not None:
        log('   |k| * |Aut| / 4 = %s' % [str(abs(sol[i]) * AUT[G] / 4) for i, G in enumerate(THETAS)])
        res = all(sum(r[i] * sol[i] for i in range(3)) == r[3] for r in eqs)
        log('   residual zero on all equations: %s' % res)
    test = [Fr(2), Fr(-4), Fr(-2, 3)]
    log('   claimed k = (2, -4, -2/3) satisfies all equations: %s' % all(sum(r[i] * test[i] for i in range(3)) == r[3] for r in eqs))


if __name__ == '__main__':
    main(sys.argv[1:])
