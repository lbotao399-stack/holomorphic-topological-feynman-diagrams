"""2d: hatQ1 is a second-order operator with vanishing first-order part:
hatQ1(abc) = hatQ1(ab) c + (-1)^|a| a hatQ1(bc) + (-1)^{(|a|+1)|b|} b hatQ1(ac)  (letters a, b, c with jets).
Also hatQ1(single letter) = 0.  Random tests.  Log: c_q1second.log"""
import random
import polyops as po
import theories
from wick import Algebra
from explicit import GraphOp
import rgraphs


def mul(p, q, pi):
    out = {}
    for m1, v1 in p.items():
        for m2, v2 in q.items():
            s, m = po.canon(m1 + m2, pi)
            if s:
                po.padd(out, m, s * v1 * v2)
    return out


for th, gname, matter in (('2d02', 'sl2', 'fund*4'), ('2d22', 'sl2', 'fund'), ('2d44', 'sl2', 'fund'), ('2d88', 'sl2', None)):
    H, sl = theories.build(th, gname, matter)
    alg = Algebra(sl, H)
    pi = alg.pi
    g = dict(rgraphs.rooted_trees(1, 1)[0]); g['name'] = 'B1'
    Q1 = GraphOp(alg, g)
    rng = random.Random(3)
    ok = nt = tests = 0
    single = all(not Q1.apply({((A, (j,)),): 1}) for A in range(alg.N) for j in range(3))
    while tests < 60:
        lets = [(rng.randrange(alg.N), (rng.randint(0, 2),)) for _ in range(3)]
        s, m = po.canon(lets, pi)
        if not s:
            continue
        a, b, c = [{(l,): 1} for l in lets]
        pa, pb = pi[lets[0][0]], pi[lets[1][0]]
        lhs = Q1.apply(mul(mul(a, b, pi), c, pi))
        rhs = po.psum(mul(Q1.apply(mul(a, b, pi)), c, pi),
                      po.pscale(mul(a, Q1.apply(mul(b, c, pi)), pi), (-1) ** pa),
                      po.pscale(mul(b, Q1.apply(mul(a, c, pi)), pi), (-1) ** ((pa + 1) * pb)))
        tests += 1
        ok += (po.psum(lhs, po.pscale(rhs, -1)) == {})
        nt += bool(lhs)
    print('%s %s %s: hatQ1(letter) = 0: %s; second-order rule on %d/%d random triples (hatQ1 != 0 on %d)'
          % (th, gname, matter, single, ok, tests, nt), flush=True)
