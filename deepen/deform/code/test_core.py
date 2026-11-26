"""Tests of ce.py: d^2 = 0, comparison of d with anomaly-structure formal_gl.d_formal and with
cochain.d_CE (explicit tensors), {s3, .} = -d, {s3, s3} = 0, symmetry and Leibniz rule of the
bracket, bracket against explicit dual-basis contraction (gl_N and sl_N), Moyal identities."""
import sys, random, time
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/deform/code')
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/anomaly-structure/code')
import ce, evalx
import formal_gl as fg, cochain as cc

rng = random.Random(7)


def rand_cochain(k, deg, nterms=4, sl=False, Npow=True):
    raw = {}
    taus = ce.all_taus(k)
    if sl:
        taus = [t for t in taus if not ce.has_singleton(t)]
    nv = k
    for _ in range(nterms):
        tau = rng.choice(taus)
        e1 = list(rng.choice(list(ce.compositions(deg[0], nv))))
        e2 = list(rng.choice(list(ce.compositions(deg[1], nv))))
        fm = []
        for i in range(nv):
            fm += [e1[i], e2[i]]
        npl = {0: Fr(rng.randint(-3, 3))}
        if Npow:
            npl[1] = Fr(rng.randint(-3, 3))
        ce.add_term(raw, (tau, tuple(fm)), npl)
    return ce.alt_raw(raw, k, sl=sl)


def to_fg(c):
    return {k: {p: int(v) for p, v in d.items()} for k, d in c.items()}


def from_fg(c):
    return {k: {p: Fr(v) for p, v in d.items()} for k, d in c.items()}


def eq(a, b):
    return ce.clean(ce.lin((1, a), (-1, b))) == {}


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= bool(ok)
    print("%-78s %s" % (name, ok), flush=True)


t0 = time.time()
# 1. d^2 = 0
for k, deg in [(2, (1, 0)), (3, (1, 1)), (3, (2, 1)), (3, (2, 2)), (4, (1, 1))]:
    for sl in (False, True):
        w = rand_cochain(k, deg, sl=sl)
        dd = ce.d_ce(ce.d_ce(w, k, sl), k + 1, sl)
        report("d^2 = 0, k=%d deg=%s sl=%s (|w| = %d terms)" % (k, deg, sl, len(w)), dd == {} and len(w) > 0)

# 2. d against formal_gl.d_formal (sympy-based, anomaly-structure)
mom = cc.Momenta(6, 2)
for k, deg in [(3, (1, 1)), (3, (2, 2)), (4, (2, 2))]:
    w = rand_cochain(k, deg, nterms=3)
    mine = ce.d_ce(w, k)
    theirs = from_fg(fg.d_formal(mom, to_fg(w), k))
    report("d_ce == formal_gl.d_formal, k=%d deg=%s" % (k, deg), eq(mine, theirs))

# 3. alternation: alt_raw output is alternating
for k, deg in [(3, (1, 1)), (4, (2, 2))]:
    w = rand_cochain(k, deg)
    report("alt_raw output alternating, k=%d deg=%s" % (k, deg), ce.is_alternating(w, k))

# 4. {s3, .} = -d, {s3, s3} = 0
s3 = ce.s3_cochain()
report("{s3, s3} = 0 (gl)", ce.bracket(s3, 3, s3, 3) == {})
report("{s3, s3} = 0 (sl)", ce.bracket(s3, 3, s3, 3, sl=True) == {})
for k, deg in [(2, (1, 0)), (3, (1, 1)), (3, (2, 2)), (4, (2, 2))]:
    for sl in (False, True):
        w = rand_cochain(k, deg, sl=sl)
        lhs = ce.bracket(s3, 3, w, k, sl=sl)
        rhs = ce.scale(ce.d_ce(w, k, sl), -1)
        report("{s3, w} = -d w, k=%d deg=%s sl=%s" % (k, deg, sl), eq(lhs, rhs))
        lhs2 = ce.bracket(w, k, s3, 3, sl=sl)
        sgn = 1 if eq(lhs2, lhs) else (-1 if eq(lhs2, ce.scale(lhs, -1)) else 0)
        report("  {w, s3} = %s {s3, w}" % sgn, sgn != 0)

# 5. symmetry {a,b} = {b,a} (arity 3,3) and Leibniz d{a,b} = -{s3,{a,b}}
for sl in (False, True):
    a = rand_cochain(3, (1, 1), sl=sl)
    b = rand_cochain(3, (2, 1), sl=sl)
    ab = ce.bracket(a, 3, b, 3, sl=sl)
    ba = ce.bracket(b, 3, a, 3, sl=sl)
    report("{a,b} = {b,a} for arity 3,3 (sl=%s)" % sl, eq(ab, ba) and ab != {})
    lhs = ce.d_ce(ab, 4, sl)
    da, db = ce.d_ce(a, 3, sl), ce.d_ce(b, 3, sl)
    r1 = ce.bracket(da, 4, b, 3, sl=sl)
    r2 = ce.bracket(a, 3, db, 4, sl=sl)
    found = [(s1, s2) for s1 in (1, -1) for s2 in (1, -1) if eq(lhs, ce.lin((s1, r1), (s2, r2)))]
    report("d{a,b} = s1 {da,b} + s2 {a,db}: (s1,s2) = %s (sl=%s)" % (found, sl), len(found) == 1)

# 6. bracket against explicit dual-basis contraction
for (m, dm), (n, dn), N, sl in [((3, (1, 1)), (3, (1, 1)), 3, False), ((3, (1, 1)), (3, (0, 1)), 2, False),
                                 ((3, (1, 1)), (4, (1, 0)), 2, False), ((3, (1, 1)), (3, (1, 1)), 3, True),
                                 ((3, (1, 0)), (4, (1, 1)), 3, True)]:
    a = rand_cochain(m, dm, sl=sl)
    b = rand_cochain(n, dn, sl=sl)
    ab = ce.bracket(a, m, b, n, sl=sl)
    K = m + n - 2
    good = True
    for trial in range(2):
        mats = [evalx.rand_mat(N, 3, traceless=sl, rng=rng) for _ in range(K)]
        moms = evalx.rand_moms(K, 4, rng=rng)
        v1 = evalx.eval_cochain(ab, mats, moms, N)
        v2 = evalx.bracket_explicit(a, m, b, n, mats, moms, N, sl=sl)
        good &= (v1 == v2)
    report("formal bracket == explicit contraction, (m,n)=(%d,%d), N=%d, sl=%s" % (m, n, N, sl), good)

# 7. d against explicit cochain.d_CE (explicit structure constants, anomaly-structure) on gl2, sl3 values
import lie, numpy as np
for g, N, sl in [(lie.gl(2), 2, False), (lie.sl(3), 3, True)]:
    w = rand_cochain(3, (1, 1), sl=sl, Npow=False)
    dw = ce.d_ce(w, 3, sl)
    good = True
    for trial in range(3):
        mats = [evalx.rand_mat(N, 3, traceless=sl, rng=rng) for _ in range(4)]
        moms = evalx.rand_moms(4, 4, rng=rng)
        v1 = evalx.eval_cochain(dw, mats, moms, N)
        # direct definition
        v2 = Fr(0)
        import itertools
        for i, j in itertools.combinations(range(4), 2):
            rest = [r for r in range(4) if r not in (i, j)]
            C = [[a_ - b_ for a_, b_ in zip(ra, rb)] for ra, rb in zip(evalx.matmul(mats[i], mats[j]), evalx.matmul(mats[j], mats[i]))]
            pm = [moms[i][0] + moms[j][0], moms[i][1] + moms[j][1]]
            v2 += (-1) ** (i + j) * evalx.eval_cochain(w, [C] + [mats[r] for r in rest], [pm] + [moms[r] for r in rest], N)
        good &= (v1 == v2)
    report("d_ce == direct CE formula on explicit %s matrices" % g.name, good)

# 8. Moyal identities (gl_N): s_eps = tr123 e^{eps u/2} - tr132 e^{-eps u/2}, u = [p1 p2]
u = ce.bracket_poly(0, 1, 3)
d_gl = {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(1)}}
f_gl = {((0, 1, 2),): {0: Fr(1)}, ((0, 2, 1),): {0: Fr(-1)}}
jM = ce.scale(ce.canon_raw(ce.raw_from(u, d_gl), 3), Fr(1, 2))
iM = ce.scale(ce.canon_raw(ce.raw_from(ce.poly_mul(u, u), f_gl), 3), Fr(1, 8))
report("j_M, i_M alternating", ce.is_alternating(jM, 3) and ce.is_alternating(iM, 3))
report("Moyal: d j_M = 0", ce.d_ce(jM, 3) == {})
lhs = ce.scale(ce.bracket(jM, 3, jM, 3), Fr(1, 2))
rhs = ce.d_ce(iM, 3)
report("Moyal: (1/2){j_M, j_M} = d i_M   ({j_M,j_M} != 0: %s)" % (lhs != {}), eq(lhs, rhs) and lhs != {})
print("all:", ok_all, " (%.1fs)" % (time.time() - t0))
