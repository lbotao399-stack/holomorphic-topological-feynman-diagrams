"""C10 closed forms of hat Q_1, C11 {hatQ0, hatQ1} = 0 on operators not used by the explorer; su(2) and so(4).
Control: moments of the reversed order 0 < r_1 < .. < r_4 < 1 in place of 1 > r_1 > .. > r_4 > 0."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
import v4_ops as V
from v4_ops import Theory, letter as Lt, op, e, padd, pmul, poly_mono, mk, EPS4, Fr, H, log

V.LOG = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v4_c11.log', 'w')


def closed_forms(T, pairs):
    ok = True
    for (A, B) in pairs:
        col = T.pent_colour(A, B)
        gg, gc = {}, {}
        for D, cD in col.items():
            for perm, s in EPS4:
                cs = [mk(0, D[q], e(perm[q])) for q in range(4)]
                gc = padd(gc, poly_mono(cs, [], 1), Fr(-1, 12) * cD * s)
                for q in range(4):
                    odd = cs[:q] + cs[q + 1:]
                    gg = padd(gg, poly_mono(odd, [mk(1, D[q], e(perm[q]))], 1), Fr(1, 12) * cD * s * (-1) ** q)
        r1 = T.Q1(op(Lt('g', A), Lt('g', B)))
        r2 = T.Q1(op(Lt('g', A), Lt('c', B)))
        r3 = T.Q1(op(Lt('c', A), Lt('c', B)))
        log('   C10 closed forms, A=%d B=%d: hatQ1(g g) %s (%d terms), hatQ1(g c) %s (%d terms), hatQ1(c c) = 0 %s'
            % (A, B, r1 == gg, len(r1), r2 == gc, len(r2), r3 == {}))
        ok &= (r1 == gg) and (r2 == gc) and (r3 == {})
    return ok


def run(alg, n, f, tests, control=False):
    T = Theory(n, f)
    if control:
        V.pent_moment_saved = V.pent_moment
        V.pent_moment = lambda a, b, c, d: Fr(1, (a + 1) * (a + b + 2) * (a + b + c + 3) * (a + b + c + d + 4))
    log('algebra %s, dim %d%s' % (alg, n, ' [CONTROL: reversed-order moments]' if control else ''))
    if not control:
        q0ok = all(T.Q0(T.Q0(p)) == {} for _nm, p in tests[:6])
        log('   hatQ0^2 = 0 on the first 6 test operators: %s' % q0ok)
        log('   closed forms: %s' % closed_forms(T, [(0, 1), (1, 1), (2, 0)] if n == 3 else [(0, 5), (2, 2), (1, 3)]))
    allok = True
    for nm, p in tests:
        t0 = time.time()
        q1 = T.Q1(p)
        a = T.Q0(q1)
        b = T.Q1(T.Q0(p))
        s = padd(a, b)
        log('   {hatQ0,hatQ1}(%-22s) = %s | |hatQ1 p| = %d, |hatQ0 hatQ1 p| = %d, |hatQ1 hatQ0 p| = %d | %.1fs'
            % (nm, 'ZERO' if s == {} else 'NONZERO (%d terms)' % len(s), len(q1), len(a), len(b), time.time() - t0))
        allok &= (s == {})
    log('   all zero: %s' % allok)
    if control:
        V.pent_moment = V.pent_moment_saved
    return allok


def tests_su2():
    g, c = (lambda A, k=(0, 0, 0, 0): Lt('g', A, k)), (lambda A, k=(0, 0, 0, 0): Lt('c', A, k))
    return [
        ('g0', op(g(0))),
        ('d3g1', op(g(1, e(2)))),
        ('d2d4g2', op(g(2, (0, 1, 0, 1)))),
        ('g1 g2', op(g(1), g(2))),
        ('g2 g2', op(g(2), g(2))),
        ('c1 g2', op(c(1), g(2))),
        ('d3g0 d3g1', op(g(0, e(2)), g(1, e(2)))),
        ('d4g2 c0', op(g(2, e(3)), c(0))),
        ('g1 d2d3c2', op(g(1), c(2, (0, 1, 1, 0)))),
        ('d1d4g0 g2', op(g(0, (1, 0, 0, 1)), g(2))),
        ('d2g1 d4c1', op(g(1, e(1)), c(1, e(3)))),
        ('d3d3g2 g0', op(g(2, e(2, 2)), g(0))),
        ('g1 g2 d4c0', op(g(1), g(2), c(0, e(3)))),
        ('g0 g0 g2', op(g(0), g(0), g(2))),
        ('d2g0 g1 c1', op(g(0, e(1)), g(1), c(1))),
    ]


def tests_so4():
    g, c = (lambda A, k=(0, 0, 0, 0): Lt('g', A, k)), (lambda A, k=(0, 0, 0, 0): Lt('c', A, k))
    return [
        ('g0', op(g(0))),
        ('d2g3', op(g(3, e(1)))),
        ('g0 g5', op(g(0), g(5))),
        ('g1 g1', op(g(1), g(1))),
        ('g2 c3', op(g(2), c(3))),
        ('d3g0 g4', op(g(0, e(2)), g(4))),
        ('d2g1 d4c4', op(g(1, e(1)), c(4, e(3)))),
        ('g0 g3 g5', op(g(0), g(3), g(5))),
    ]


if __name__ == '__main__':
    n, f = V.su2()
    ok1 = run('su2', n, f, tests_su2())
    ok_ctrl = run('su2', n, f, [t for t in tests_su2() if t[0] in ('d3g0 d3g1', 'd1d4g0 g2', 'd2g1 d4c1')], control=True)
    n, f = V.so(4)
    ok2 = run('so4', n, f, tests_so4())
    log('SUMMARY: su2 all zero %s; so4 all zero %s' % (ok1, ok2))
