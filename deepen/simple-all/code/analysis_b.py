"""Part (b) analysis of logs/b_data.json (exact rationals).

1. Identities checked on every algebra (Killing normalization):
   e_k = (-1)^k d (k = 1..4);  m1 = 0, m2 = d, m3 = -d/4;  n1 = d, n2 = d/2, n3 = 0;  m4 - n4 = d/8
   W133: kk = (d, d/4, 0), ff = (-d/8, -m4, -n4);   W223: kk = (-d/2, d/4, d/4), ff = (0, n4, n4)
   Gram: <kk_i,kk_j> = d^2 delta_ij + d (1 - delta_ij);  <kk_i,ff_j>, <ff_i,ff_j> as printed (universal multiples of d).
2. Combinatorial transfer: v(F_G)[mu, Y_pi] = sum_s sgn(s) c_mu(s) <W_G, Y_pi o s^{-1}>, with Y_pi o t = sign * Y_pi'
   (kk: sign +1; ff: antisymmetry inside a pair, symmetry under exchange of the pairs); same for beta (Y = ff
   permuted by (i,j,k,l)) and Gkk.  Checked against the directly computed v on every algebra.
3. Solution of v(F223) = x v(beta) + a v(Gkk) in terms of (d, n4); beta-annihilating functional mu.
"""
import json, sys, itertools
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/deepen/simple-all/code')
import numpy as np
import t4

P3 = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]


def canon_pair(a, b, c, e):
    """kk: {a,b}|{c,e} -> index of P3."""
    s1 = {a, b}
    if 0 not in s1:
        s1 = {c, e}
    other = (s1 - {0}).pop()
    return {1: 0, 2: 1, 3: 2}[other]


def canon_ff(a, b, c, e):
    """kappa([x_a,x_b],[x_c,x_e]) = sign * ff_{P3[k]}."""
    sign = 1
    if 0 not in (a, b):
        a, b, c, e = c, e, a, b
    if a != 0:
        a, b = b, a; sign = -sign
    if c > e:
        c, e = e, c; sign = -sign
    k = {1: 0, 2: 1, 3: 2}[b]
    assert (c, e) == tuple(sorted(set(range(4)) - {0, b}))
    return k, sign


def inv(s):
    out = [0] * 4
    for i, x in enumerate(s):
        out[x] = i
    return tuple(out)


def transfer_alt():
    """matrices T_kk[mu][pi][pi'] and T_ff[mu][pi][pi'] for X = alt([p0p1][p2p3] W)."""
    Tk = np.zeros((6, 3, 3), dtype=np.int64)
    Tf = np.zeros((6, 3, 3), dtype=np.int64)
    for sig in itertools.permutations(range(4)):
        c = t4.quad(t4.br(sig[0], sig[1]), t4.br(sig[2], sig[3]))
        sg = t4.perm_sign(sig)
        si = inv(sig)
        for mu in range(6):
            if not c[mu]:
                continue
            for i, p in enumerate(P3):
                # Y_p o s^{-1}: indices s^{-1}(p_k)
                q = tuple(si[x] for x in p)
                Tk[mu, i, canon_pair(*q)] += sg * c[mu]
                k, sgn = canon_ff(*q)
                Tf[mu, i, k] += sg * c[mu] * sgn
    return Tk, Tf


def transfer_beta():
    """beta^mu = sum_{i<j} (-1)^{i+j} quad([p_k p_l]^2)_mu ff o (i,j,k,l):  <ff o t, Y> = <ff, Y o t^{-1}>."""
    Tk = np.zeros((6, 3, 3), dtype=np.int64)        # [mu, probe pi, gram index of <ff_pi', Y>]
    Tf = np.zeros((6, 3, 3), dtype=np.int64)
    for i, j in itertools.combinations(range(4), 2):
        k, l = [r for r in range(4) if r not in (i, j)]
        c = t4.quad(t4.br(k, l), t4.br(k, l))
        t = (i, j, k, l)
        # (ff o t)(x) = kappa([x_i,x_j],[x_k,x_l]) = sign ff_{pi'}(x)
        kk_, sgn = canon_ff(*t)
        for mu in range(6):
            if c[mu]:
                for pi in range(3):
                    Tk[mu, pi, kk_] += (-1) ** (i + j) * c[mu] * sgn
                    Tf[mu, pi, kk_] += (-1) ** (i + j) * c[mu] * sgn
    return Tk, Tf


def main(path):
    recs = json.load(open(path))
    Tk, Tf = transfer_alt()
    print("transfer (alternation), mu = uv block: kk rows (probe pi) x cols (W.kk pi'):", Tk[3].tolist(),
          " ff:", Tf[3].tolist())
    d_, n4_ = sp.symbols('d n4')
    W133s = {'kk': [d_, d_ / 4, 0], 'ff': [-d_ / 8, -(n4_ + d_ / 8), -n4_]}
    W223s = {'kk': [-d_ / 2, d_ / 4, d_ / 4], 'ff': [0, n4_, n4_]}
    rows = []
    allok = True
    for R in recs:
        d = R['dim']
        m = {k: Fr(v) for k, v in R['moments'].items()}
        W = {k: Fr(v) for k, v in R['W'].items()}
        g = {k: Fr(v) for k, v in R['gram'].items()}
        v = {k: [Fr(x) for x in vv] for k, vv in R['v'].items()}
        n4 = m['n4']
        chk = {}
        chk['e_k'] = all(m['e%d' % k] == (-1) ** k * d for k in range(1, 5))
        chk['m1-3,n1-3'] = (m['m1'], m['m2'], m['m3'], m['n1'], m['n2'], m['n3']) == (0, d, Fr(-d, 4), d, Fr(d, 2), 0)
        chk['m4-n4=d/8'] = m['m4'] - m['n4'] == Fr(d, 8)
        chk['W133'] = ([W['W133.kk%d' % i] for i in range(3)] == [d, Fr(d, 4), 0] and
                       [W['W133.ff%d' % i] for i in range(3)] == [Fr(-d, 8), -m['m4'], -n4])
        chk['W223'] = ([W['W223.kk%d' % i] for i in range(3)] == [Fr(-d, 2), Fr(d, 4), Fr(d, 4)] and
                       [W['W223.ff%d' % i] for i in range(3)] == [0, n4, n4])
        chk['gram kk.kk'] = all(g['kk%d.kk%d' % (i, j)] == (d * d if i == j else d) for i in range(3) for j in range(3))
        # transfer check
        ok_tr = True
        for G in ('133', '223'):
            wk = [W['W%s.kk%d' % (G, i)] for i in range(3)]
            wf = [W['W%s.ff%d' % (G, i)] for i in range(3)]
            for mu in range(6):
                for pi in range(3):
                    pk = sum(int(Tk[mu, pi, q]) * wk[q] for q in range(3))
                    pf = sum(int(Tf[mu, pi, q]) * wf[q] for q in range(3))
                    ok_tr &= (pk == v['F' + G][6 * mu + pi]) and (pf == v['F' + G][6 * mu + 3 + pi])
        # Gkk via transfer with W = kk: W.kk_pi' = <kk_0123, kk_pi'>, W.ff_pi' = <kk_0123, ff_pi'>
        wk = [g['kk0.kk%d' % i] for i in range(3)]
        wf = [g['kk0.ff%d' % i] for i in range(3)]
        for mu in range(6):
            for pi in range(3):
                pk = sum(int(Tk[mu, pi, q]) * wk[q] for q in range(3))
                pf = sum(int(Tf[mu, pi, q]) * wf[q] for q in range(3))
                ok_tr &= (pk == v['Gkk'][6 * mu + pi]) and (pf == v['Gkk'][6 * mu + 3 + pi])
        chk['transfer F, Gkk'] = ok_tr
        # x, a
        lk = lambda X: v[X][6 * 3 + 0]          # mu = uv, kk_0123
        lf = lambda X: v[X][6 * 3 + 3]          # mu = uv, ff_0123
        mu_F = lk('F223') * lf('beta') - lf('F223') * lk('beta')
        mu_G = lk('Gkk') * lf('beta') - lf('Gkk') * lk('beta')
        chk['probe values'] = (lk('beta'), lf('beta'), lk('Gkk'), lf('Gkk'), lk('F223'), lf('F223')) == \
            (2 * d, d, 8 * d * (d - 1), 8 * d, -6 * d, -8 * n4)
        chk['mu(F223) = 2d(8n4-3d)'] = mu_F == 2 * d * (8 * n4 - 3 * d)
        sol = None
        M = sp.Matrix([[sp.Rational(b.numerator, b.denominator), sp.Rational(c.numerator, c.denominator)]
                       for b, c in zip(v['beta'], v['Gkk'])])
        y = sp.Matrix([sp.Rational(f.numerator, f.denominator) for f in v['F223']])
        rk = M.rank(); rk2 = M.row_join(y).rank()
        if rk == 2 and rk2 == 2:
            sol = (M.T * M).inv() * M.T * y
            assert (M * sol - y).is_zero_matrix
        line = "%-16s d %3d  n4 %-14s 8n4-3d %-12s mu(F223) %-14s rank[v(beta),v(Gkk)] %d rank[+v(F223)] %d" % (
            R['name'], d, n4, 8 * n4 - 3 * d, mu_F, rk, rk2)
        if sol is not None:
            a_formula = (8 * n4 - 3 * d) / Fr(4 * d * (d - 3)) if d != 3 else None
            line += "  x = %s, a = %s, a = (8n4-3d)/(4d(d-3)): %s" % (sol[0], sol[1], sp.Rational(a_formula.numerator, a_formula.denominator) == sol[1])
        rows.append(line)
        bad = [k for k, x in chk.items() if not x]
        allok &= not bad
        print("%-16s checks: %s" % (R['name'], ", ".join("%s %s" % (k, x) for k, x in chk.items())))
    print()
    for l in rows:
        print(l)
    print("\nall identities hold on all %d algebras: %s" % (len(recs), allok))
    print("Gram (Killing), e.g. %s: " % recs[-1]['name'] + " ".join("%s=%s" % (k, v) for k, v in recs[-1]['gram'].items()))


if __name__ == "__main__":
    main(sys.argv[1])
