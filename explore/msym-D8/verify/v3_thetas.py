"""C6-C9: theta kernels on C^4, independent of htk and of the explorer's code.

Coordinates (three-chain move for ALL three thetas; the explorer used two moves for (1,4,4), (2,3,4)):
  c_j = (1/R_j) / sum_i (1/R_i)  (c_1 = 1 - c_2 - c_3),  x_{j,k} = t_{j,k} / R_j  (x_{j,1} = 1 - sum_{k>=2} x_{j,k}),
  section t_{j,k} = x_{j,k} / c_j;   xi = (c_2, c_3, x_{1,2..p}, x_{2,2..q}, x_{3,2..r}).
Own formula (Kirchhoff on the theta, lam_N = 0): with sigma_v = sum_{k > i(v)} x_{j,k} for v = v_{j,i},
  l_j := y(spoke of P_j) = - c_j sum_v (1 - sigma_v) lam_v - sum_{v in P_j} sigma_v lam_v ;
  {spoke of P2, spoke of P3} is the cotree of the tree P1 u (P2 - spoke) u (P3 - spoke), so
  g(xi) := D(t(xi)) |det[t; d_xi t]| = sgn * det[ d(l_2, l_3) / d xi ],  sgn = sign det[t; d_xi t].
Checks:
  (i)   C6 displacement map, exact, fresh mesh solver, random rational t;
  (ii)  D(t(xi)) det[t; d_xi t] == det[d(l_2,l_3)/d xi] exactly at random rational xi (fresh solver D, BFS tree);
  (iii) claimed densities of C7, C8 (explorer's two-move sections) and C9 (det L) == D |det[t;d t]| exactly;
  (iv)  K(lam;0) by exact integration of det[d l/d xi] (numeric lam, sympy DomainMatrix det) == claimed closed forms;
  (v)   bracket expansion K = sum C_{VV'} <V><V'> and first-order weights (own), consistency with (iv);
  (vi)  float tensor-Gauss quadrature of the fresh-solver density D |det[t; d t]| (z = 0 and first order in z,
        first order from the fresh-solver spoke currents) vs the exact values;
  (vii) symmetries: P2 <-> P3 for (1,4,4); path swaps for (3,3,3).
Output: v3_thetas.log, theta_weights.pkl (weights for v4_ops.py).
"""
import sys, random, itertools, time, pickle, math
from fractions import Fraction as Fr
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify')
import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from common import Network, theta, bracket, product_rule, perm_sign, H

OUT = open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/v3_thetas.log', 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    OUT.write(s + '\n'); OUT.flush()


THETAS = [(1, 4, 4), (2, 3, 4), (3, 3, 3)]


def closed_form(pqr, lam):
    b = lambda *S: bracket(lam, S)
    if pqr == (1, 4, 4):
        return Fr(-1, 11520) * (b(1, 2, 3, 4) + 2 * b(1, 2, 3, 5) + 3 * b(1, 2, 3, 6)) * (b(1, 4, 5, 6) + 2 * b(2, 4, 5, 6) + 3 * b(3, 4, 5, 6))
    if pqr == (2, 3, 4):
        return Fr(-1, 17280) * (b(1, 2, 3, 4) + 2 * b(1, 2, 3, 5) + 3 * b(1, 2, 3, 6)) * (3 * b(1, 4, 5, 6) + 2 * b(2, 4, 5, 6) + 4 * b(3, 4, 5, 6))
    terms = [(-1, (1, 2, 3, 4), (1, 3, 5, 6)), (-2, (1, 2, 3, 4), (1, 4, 5, 6)), (-2, (1, 2, 3, 4), (2, 3, 5, 6)),
             (-4, (1, 2, 3, 4), (2, 4, 5, 6)), (-1, (1, 2, 3, 5), (1, 4, 5, 6)), (-2, (1, 2, 3, 5), (2, 4, 5, 6)),
             (-1, (1, 2, 3, 5), (3, 4, 5, 6)), (-2, (1, 2, 3, 6), (1, 4, 5, 6)), (-4, (1, 2, 3, 6), (2, 4, 5, 6)),
             (-2, (1, 2, 3, 6), (3, 4, 5, 6)), (1, (1, 2, 4, 5), (1, 3, 5, 6)), (2, (1, 2, 4, 5), (2, 3, 5, 6)),
             (-2, (1, 2, 4, 5), (3, 4, 5, 6)), (2, (1, 2, 4, 6), (1, 3, 5, 6)), (4, (1, 2, 4, 6), (2, 3, 5, 6)),
             (-4, (1, 2, 4, 6), (3, 4, 5, 6))]
    return Fr(1, 8640) * sum(c * b(*S) * b(*T) for c, S, T in terms)


class Theta:
    def __init__(self, pqr):
        self.pqr = pqr
        self.nv, self.edges, self.paths = theta(*pqr)
        self.net = Network(self.nv, self.edges)
        self.outputs = [v for _e, ints in self.paths for v in ints]
        # symbols
        self.cs = [sp.Symbol('c2'), sp.Symbol('c3')]
        self.xs = {}
        for j, L in enumerate(pqr):
            for k in range(2, L + 1):
                self.xs[(j, k)] = sp.Symbol('x%d_%d' % (j + 1, k))
        self.xi = self.cs + [self.xs[(j, k)] for j, L in enumerate(pqr) for k in range(2, L + 1)]
        self.blocks = [self.cs] + [[self.xs[(j, k)] for k in range(2, L + 1)] for j, L in enumerate(pqr) if L >= 2]
        c = [1 - self.cs[0] - self.cs[1], self.cs[0], self.cs[1]]
        self.c = c
        xfull = {}
        for j, L in enumerate(pqr):
            rest = [self.xs[(j, k)] for k in range(2, L + 1)]
            xfull[(j, 1)] = 1 - sum(rest)
            for k in range(2, L + 1):
                xfull[(j, k)] = self.xs[(j, k)]
        self.xfull = xfull
        # sigma_v, path of v
        self.sigma, self.jof = {}, {}
        for j, (eids, ints) in enumerate(self.paths):
            L = len(eids)
            for i, v in enumerate(ints, start=1):          # v = v_{j,i}
                self.sigma[v] = sum(xfull[(j, k)] for k in range(i + 1, L + 1))
                self.jof[v] = j
        # response R[j][v] = coefficient of lam_v in l_j
        self.R = {j: {v: sp.expand(-c[j] * (1 - self.sigma[v]) - (self.sigma[v] if self.jof[v] == j else 0)) for v in self.outputs}
                  for j in range(3)}
        # displacement coefficients X_v = sum_j cX[v][j] z_j
        self.cX = {v: [sp.expand((self.sigma[v] if self.jof[v] == j else 0) + (1 - self.sigma[v]) * c[j]) for j in range(3)]
                   for v in self.outputs}
        # section
        self.tex = [None] * len(self.edges)
        for j, (eids, ints) in enumerate(self.paths):
            for k, e in enumerate(eids, start=1):
                self.tex[e] = xfull[(j, k)] / c[j]
        self.Jt = sp.Matrix([self.tex] + [[sp.diff(x, s) for x in self.tex] for s in self.xi])
        self.Jt_f = sp.lambdify(self.xi, list(self.Jt), 'numpy')
        self.t_f = sp.lambdify(self.xi, self.tex, 'numpy')

    def random_xi(self, rng):
        vals = {}
        a = [Fr(rng.randint(1, 60)) for _ in range(3)]
        sa = sum(a)
        vals[self.cs[0]], vals[self.cs[1]] = a[1] / sa, a[2] / sa
        for j, L in enumerate(self.pqr):
            if L >= 2:
                b = [Fr(rng.randint(1, 60)) for _ in range(L)]
                sb = sum(b)
                for k in range(2, L + 1):
                    vals[self.xs[(j, k)]] = b[k - 1] / sb
        return vals

    def subsQ(self, expr, vals):
        r = sp.Rational(expr.subs({s: sp.Rational(v.numerator, v.denominator) for s, v in vals.items()}))
        return Fr(int(r.p), int(r.q))

    def ell(self, lam, j):
        return [sum(self.R[j][v] * lam[v][a] for v in self.outputs) for a in range(H)]

    def jac_ell(self, lam):
        rows = []
        for j in (1, 2):
            lj = self.ell(lam, j)
            for a in range(H):
                rows.append([sp.diff(lj[a], s) for s in self.xi])
        return sp.Matrix(rows)

    def integrate(self, poly_expr):
        P = sp.Poly(sp.expand(poly_expr), *self.xi)
        pos = {s: i for i, s in enumerate(self.xi)}
        tot = Fr(0)
        for mon, cf in P.terms():
            w = Fr(1)
            for blk in self.blocks:
                ex = [mon[pos[s]] for s in blk]
                num = 1
                for x in ex:
                    num *= math.factorial(x)
                w *= Fr(num, math.factorial(len(blk) + sum(ex)))
            tot += Fr(int(cf.p), int(cf.q)) * w
        return tot


def section_two_moves(th, rng, which):
    """Explorer's sections of C7 / C8: (P1||P2) -> a, a || P3, gauge t_b = 1:
    t_a = 1/u, R_3 = 1/(1-u), R_1 = t_a/om_1, R_2 = t_a/(1-om_1), om_1 = t_a/R_1; t_{j,k} = x_{j,k} R_j."""
    om = sp.Symbol('om'); u = sp.Symbol('u')
    names = [om]
    xs = {}
    for j, L in enumerate(th.pqr):
        for k in range(2, L + 1):
            xs[(j, k)] = sp.Symbol('y%d_%d' % (j + 1, k))
    coords = [om] + [xs[(j, k)] for j in (0, 1) for k in range(2, th.pqr[j] + 1)] + [u] + [xs[(2, k)] for k in range(2, th.pqr[2] + 1)]
    ta = 1 / u
    R = [ta / om, ta / (1 - om), 1 / (1 - u)]
    tex = [None] * len(th.edges)
    xfull = {}
    for j, L in enumerate(th.pqr):
        rest = [xs[(j, k)] for k in range(2, L + 1)]
        xfull[(j, 1)] = 1 - sum(rest)
        for k in range(2, L + 1):
            xfull[(j, k)] = xs[(j, k)]
    for j, (eids, ints) in enumerate(th.paths):
        for k, e in enumerate(eids, start=1):
            tex[e] = xfull[(j, k)] * R[j]
    return coords, tex, xfull, om, u


def check_two_move_density(th, rng, npts=4):
    coords, tex, xfull, om, u = section_two_moves(th, rng, th.pqr)
    Jm = sp.Matrix([tex] + [[sp.diff(x, s) for x in tex] for s in coords])
    ok = True
    for it in range(npts):
        vals = {om: Fr(rng.randint(1, 99), 100), u: Fr(rng.randint(1, 99), 100)}
        for j, L in enumerate(th.pqr):
            if L >= 2:
                b = [Fr(rng.randint(1, 60)) for _ in range(L)]
                sb = sum(b)
                for k in range(2, L + 1):
                    vals[sp.Symbol('y%d_%d' % (j + 1, k))] = b[k - 1] / sb
        sub = {s: sp.Rational(v.numerator, v.denominator) for s, v in vals.items()}
        t = [Fr(int(sp.Rational(x.subs(sub)).p), int(sp.Rational(x.subs(sub)).q)) for x in tex]
        J = sp.Rational(Jm.subs(sub).det())
        lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in th.outputs}
        D, y = th.net.density(t, lam)
        g_brute = D * abs(Fr(int(J.p), int(J.q)))
        # sigma_v, mu
        sig = {}
        for j, (eids, ints) in enumerate(th.paths):
            L = len(eids)
            for i, v in enumerate(ints, start=1):
                s_ = sum(xfull[(j, k)] for k in range(i + 1, L + 1))
                sig[v] = Fr(int(sp.Rational(sp.sympify(s_).subs(sub)).p), int(sp.Rational(sp.sympify(s_).subs(sub)).q))
        mu = lambda Ps: [sum((1 - sig[v]) * lam[v][a] for j in Ps for v in th.paths[j][1]) for a in range(H)]
        om1, uu = vals[om], vals[u]
        det4 = lambda cols: Fr(int(sp.Matrix([[sp.Rational(c[a].numerator, c[a].denominator) if isinstance(c[a], Fr) else c[a] for c in cols] for a in range(H)]).det().p),
                               int(sp.Matrix([[sp.Rational(c[a].numerator, c[a].denominator) if isinstance(c[a], Fr) else c[a] for c in cols] for a in range(H)]).det().q))
        L_ = lambda v: lam[v]
        if th.pqr == (1, 4, 4):
            g_claim = -om1 ** 3 * uu ** 4 * det4([L_(1), L_(2), L_(3), mu([2])]) * det4([mu([1]), L_(4), L_(5), L_(6)])
        else:
            om2 = 1 - om1
            g_claim = -om1 ** 2 * om2 * uu ** 4 * det4([L_(1), L_(2), L_(3), mu([2])]) * det4([mu([0, 1]), L_(4), L_(5), L_(6)])
        log('   (iii) %s two-move section point %d: D|det[t;dt]| = %s, claimed g = %s, equal: %s'
            % (th.pqr, it, g_brute, g_claim, g_brute == g_claim))
        ok &= g_brute == g_claim
    return ok


def main():
    rng = random.Random(4242)
    weights = {}
    t_all = time.time()
    for pqr in THETAS:
        t0 = time.time()
        th = Theta(pqr)
        log('=' * 90)
        log('theta%s: edges %s; paths %s; BFS tree %s, cotree %s' % (pqr, th.edges, th.paths, th.net.tree, th.net.cotree))
        # (i) displacement map
        ok = True
        for it in range(5):
            t = [Fr(rng.randint(1, 40), rng.randint(1, 9)) for _ in range(len(th.edges))]
            lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in th.outputs}
            z = [[Fr(rng.randint(-9, 9), rng.randint(1, 4)) for _ in range(H)] for _ in range(3)]
            y, ell, _ = th.net.solve(t, lam)
            assert th.net.check_kcl(y, lam)
            lhs = -sum(y[th.paths[j][0][0]][a] * z[j][a] for j in range(3) for a in range(H))
            R = [sum(t[e] for e in th.paths[j][0]) for j in range(3)]
            cj = [(1 / R[j]) / sum(1 / x for x in R) for j in range(3)]
            ZN = [sum(cj[j] * z[j][a] for j in range(3)) for a in range(H)]
            rhs = Fr(0)
            for j, (eids, ints) in enumerate(th.paths):
                for i, v in enumerate(ints, start=1):
                    sig = sum(t[e] for e in eids[i:]) / R[j]
                    rhs += sum(lam[v][a] * (sig * z[j][a] + (1 - sig) * ZN[a]) for a in range(H))
            ok &= (lhs == rhs)
        log('(i) C6 displacement map: -sum_e y_e.z_e == sum_v lam_v.(sigma_v z_j(v) + (1-sigma_v) Z_N), 5 rational points: %s' % ok)

        # (ii) pointwise: D(t(xi)) det[t; d t] == det[d l / d xi]
        ok = True
        sgns = set()
        lam_sym = None
        for it in range(4):
            vals = th.random_xi(rng)
            t = [th.subsQ(x, vals) for x in th.tex]
            J = th.subsQ(th.Jt.det() if it == 0 else th.Jt.det(), vals) if False else None
            Jm = th.Jt.subs({s: sp.Rational(v.numerator, v.denominator) for s, v in vals.items()})
            J = sp.Rational(Jm.det())
            J = Fr(int(J.p), int(J.q))
            sgns.add(1 if J > 0 else -1)
            lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in th.outputs}
            D, y = th.net.density(t, lam)
            # own formula for the spoke currents vs the solver
            okc = all(th.subsQ(th.ell(lam, j)[a], vals) == y[th.paths[j][0][0]][a] for j in range(3) for a in range(H))
            Mj = th.jac_ell(lam).subs({s: sp.Rational(v.numerator, v.denominator) for s, v in vals.items()})
            dl = sp.Rational(Mj.det())
            dl = Fr(int(dl.p), int(dl.q))
            log('(ii) point %d: spoke currents = own formula: %s; D(t) det[t;dt/dxi] = %s, det[d(l2,l3)/dxi] = %s, equal: %s; sign det[t;dt/dxi] = %+d'
                % (it, okc, D * J, dl, D * J == dl, 1 if J > 0 else -1))
            ok &= okc and (D * J == dl)
        assert len(sgns) == 1
        sgn = sgns.pop()
        log('(ii) all points: %s; orientation sign of the section: %+d' % (ok, sgn))

        # (iii) claimed densities
        if pqr in ((1, 4, 4), (2, 3, 4)):
            okd = check_two_move_density(th, rng)
            log('(iii) claimed two-move density of %s == brute force at all points: %s' % ('C7' if pqr == (1, 4, 4) else 'C8', okd))
        else:
            okL = True
            for it in range(4):
                vals = th.random_xi(rng)
                lam = {v: [rng.randint(-9, 9) for _ in range(H)] for v in th.outputs}
                t = [th.subsQ(x, vals) for x in th.tex]
                Jm = th.Jt.subs({s: sp.Rational(v.numerator, v.denominator) for s, v in vals.items()})
                J = sp.Rational(Jm.det()); J = Fr(int(J.p), int(J.q))
                D, y = th.net.density(t, lam)
                w = [th.subsQ(sp.sympify(cc), vals) for cc in th.c]
                sig = {v: th.subsQ(sp.sympify(th.sigma[v]), vals) for v in th.outputs}
                yred = [-sum((1 - sig[v]) * lam[v][a] for v in th.outputs) for a in range(H)]
                ints = [p_[1] for p_ in th.paths]
                Qv = {}
                for j in range(3):
                    Qv[(j, 2)] = lam[ints[j][0]]
                    Qv[(j, 3)] = [lam[ints[j][0]][a] + lam[ints[j][1]][a] for a in range(H)]
                rows = []
                for i in (1, 2):
                    for a in range(H):
                        row = [yred[a] if i == ip else Fr(0) for ip in (1, 2)]
                        row += [(w[i] - (1 if i == j else 0)) * Qv[(j, l)][a] for j in range(3) for l in (2, 3)]
                        rows.append(row)
                from common import det_exact
                dL = det_exact(rows)
                log('   (iii) det L at point %d: D|det[t;dt]| = %s, det L = %s, equal: %s' % (it, D * abs(J), dL, D * abs(J) == dL))
                okL &= (D * abs(J) == dL)
            log('(iii) C9 density g = det L at all points: %s' % okL)

        # (iv) exact K(lam;0) at numeric lam, direct det, vs closed forms
        dom = sp.QQ[tuple(th.xi)]
        okiv = True
        Kvals = []
        for it in range(3):
            lam = {v: [rng.randint(-7, 7) for _ in range(H)] for v in th.outputs}
            Mj = th.jac_ell(lam)
            DM = DomainMatrix([[dom.from_sympy(sp.expand(x)) for x in Mj.row(r)] for r in range(8)], (8, 8), dom)
            dpoly = dom.to_sympy(DM.det())
            K0 = sgn * th.integrate(dpoly)
            cf = closed_form(pqr, lam)
            log('(iv) lam %d: exact integral of sgn*det[dl/dxi] = %s, claimed closed form = %s, equal: %s' % (it, K0, cf, K0 == cf))
            okiv &= (K0 == cf)
            Kvals.append((lam, K0))
        log('(iv) all: %s' % okiv)

        # (v) bracket expansion + first-order weights (polynomial ring QQ[xi])
        tv = time.time()
        ring = sp.QQ.old_poly_ring(*th.xi) if False else None
        P_ = lambda e: sp.Poly(sp.expand(e), *th.xi, domain='QQ')
        Bm = {j: {v: [P_(sp.diff(th.R[j][v], s)) for s in th.xi] for v in th.outputs} for j in (1, 2)}
        subsets = list(itertools.combinations(th.outputs, 4))
        Ssets = list(itertools.combinations(range(8), 4))
        perms4 = [(p, perm_sign(p)) for p in itertools.permutations(range(4))]

        def det4(rows):
            tot = None
            for p, s in perms4:
                term = rows[0][p[0]] * rows[1][p[1]] * rows[2][p[2]] * rows[3][p[3]]
                if term.is_zero:
                    continue
                term = term if s > 0 else -term
                tot = term if tot is None else tot + term
            return tot
        minors = {j: {} for j in (1, 2)}
        for j in (1, 2):
            for V in subsets:
                for S in Ssets:
                    d = det4([[Bm[j][v][k] for k in S] for v in V])
                    if d is not None and not d.is_zero:
                        minors[j][(V, S)] = d
        cVV = {}
        for S in Ssets:
            Sc = tuple(k for k in range(8) if k not in S)
            s_ = (-1) ** (10 + sum(k + 1 for k in S))
            for V in subsets:
                a = minors[1].get((V, S))
                if a is None:
                    continue
                for Vp in subsets:
                    b = minors[2].get((Vp, Sc))
                    if b is None:
                        continue
                    term = a * b if s_ > 0 else -(a * b)
                    cVV[(V, Vp)] = cVV[(V, Vp)] + term if (V, Vp) in cVV else term
        cVV = {k: v for k, v in cVV.items() if not v.is_zero}
        C0 = {k: sgn * th.integrate(v.as_expr()) for k, v in cVV.items()}
        C0 = {k: v for k, v in C0.items() if v != 0}
        cXP = {w_: [P_(th.cX[w_][j]) for j in range(3)] for w_ in th.outputs}
        C1 = {}
        for k, v in cVV.items():
            for w_ in th.outputs:
                for j in range(3):
                    val = sgn * th.integrate((v * cXP[w_][j]).as_expr())
                    if val != 0:
                        C1[(k, w_, j)] = val
        log('(v) %d nonzero coefficient polynomials c_{VV\'}(xi); %d nonzero integrals C_{VV\'}; %d first-order integrals (%.1fs)'
            % (len(cVV), len(C0), len(C1), time.time() - tv))

        def K_brackets(lam):
            return sum(c * bracket(lam, V) * bracket(lam, Vp) for (V, Vp), c in C0.items())

        def K1_brackets(lam, z):
            return sum(c * bracket(lam, V) * bracket(lam, Vp) * sum(lam[w_][a] * z[j][a] for a in range(H))
                       for ((V, Vp), w_, j), c in C1.items())
        okv = all(K_brackets(lam) == K0 for lam, K0 in Kvals)
        log('(v) bracket expansion reproduces the direct integrals of (iv): %s' % okv)
        # direct first order at one lam, z: int sgn*det * sum_v lam_v . X_v
        lam, _ = Kvals[0]
        z = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(3)]
        Mj = th.jac_ell(lam)
        DM = DomainMatrix([[dom.from_sympy(sp.expand(x)) for x in Mj.row(r)] for r in range(8)], (8, 8), dom)
        dpoly = dom.to_sympy(DM.det())
        Xdot = sum(lam[w_][a] * th.cX[w_][j] * z[j][a] for w_ in th.outputs for j in range(3) for a in range(H))
        K1_direct = sgn * th.integrate(sp.expand(dpoly * Xdot))
        log('(v) first order in z, lam 0, z = %s: direct = %s, bracket weights = %s, equal: %s'
            % (z, K1_direct, K1_brackets(lam, z), K1_direct == K1_brackets(lam, z)))

        # (vi) float quadrature of the brute-force density
        tq = time.time()
        dims = [2] + [L - 1 for L in pqr if L >= 2]
        nodes, Wq = product_rule(dims, 5)
        cols = np.concatenate(nodes, axis=1)
        args = [cols[:, i] for i in range(cols.shape[1])]
        lam, K0 = Kvals[1]
        lamf = {v: [float(x) for x in lam[v]] for v in th.outputs}
        Pn = cols.shape[0]
        Dl, Jl, yl = [], [], []
        for s0 in range(0, Pn, 40000):
            a_ = [x[s0:s0 + 40000] for x in args]
            m_ = len(a_[0])
            Tm = np.array([np.broadcast_to(np.asarray(x, dtype=float), (m_,)) for x in th.t_f(*a_)]).T
            Jfull = th.Jt_f(*a_)
            Jarr = np.empty((m_, 9, 9))
            for r in range(9):
                for c_ in range(9):
                    Jarr[:, r, c_] = np.broadcast_to(np.asarray(Jfull[9 * r + c_], dtype=float), (m_,))
            Jl.append(np.abs(np.linalg.det(Jarr)))
            Dv_, yv_ = th.net.density_float(Tm, lamf)
            Dl.append(Dv_)
            yl.append(yv_[:, [th.paths[j][0][0] for j in range(3)], :])
        Dv = np.concatenate(Dl); Jd = np.concatenate(Jl); yv = np.concatenate(yl)
        val0 = np.sum(Wq * Dv * Jd)
        rel0 = abs(val0 - float(K0)) / abs(float(K0))
        z = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(3)]
        lin = -sum(yv[:, j, a] * z[j][a] for j in range(3) for a in range(H))
        val1 = np.sum(Wq * Dv * Jd * lin)
        ex1 = K1_brackets(lam, z)
        rel1 = abs(val1 - float(ex1)) / abs(float(ex1))
        log('(vi) quadrature (5 Gauss nodes per direction, %d points, fresh float solver): K(lam;0) = %.14g vs exact %s = %.14g, rel %.1e; '
            'first order (spoke currents of the solver) = %.14g vs exact %s, rel %.1e (%.1fs)'
            % (Pn, val0, K0, float(K0), rel0, val1, ex1, rel1, time.time() - tq))

        # (vii) symmetries
        def apply_perm(lam, z, vperm, jperm):
            return {v: lam[vperm[v]] for v in th.outputs}, [z[jperm[j]] for j in range(3)]
        if pqr == (1, 4, 4):
            perms = [({1: 4, 2: 5, 3: 6, 4: 1, 5: 2, 6: 3}, [0, 2, 1])]
        elif pqr == (3, 3, 3):
            perms = [({1: 3, 2: 4, 3: 1, 4: 2, 5: 5, 6: 6}, [1, 0, 2]), ({1: 1, 2: 2, 3: 5, 4: 6, 5: 3, 6: 4}, [0, 2, 1]),
                     ({1: 5, 2: 6, 3: 3, 4: 4, 5: 1, 6: 2}, [2, 1, 0])]
        else:
            perms = []
        for vp, jp in perms:
            res = []
            for it in range(3):
                lam = {v: [rng.randint(-7, 7) for _ in range(H)] for v in th.outputs}
                z = [[rng.randint(-5, 5) for _ in range(H)] for _ in range(3)]
                l2, z2 = apply_perm(lam, z, vp, jp)
                a0, b0 = K_brackets(l2), K_brackets(lam)
                a1, b1 = K1_brackets(l2, z2), K1_brackets(lam, z)
                res.append((a0 / b0 if b0 else None, a1 / b1 if b1 else None))
            log('(vii) %s: vertex map %s, input map %s: ratios K(lam o pi; z o pi)/K(lam; z) at z^0, z^1 (3 points): %s'
                % (pqr, vp, [j + 1 for j in jp], [(str(x), str(y_)) for x, y_ in res]))

        # weights in monomials for v4_ops: m = tuple over outputs (vertex order) of 4-exponent tuples
        tw = time.time()
        idx = {v: i for i, v in enumerate(th.outputs)}
        bm = {}

        def bmon(V):
            if V not in bm:
                out = {}
                for perm in itertools.permutations(range(4)):
                    key = tuple(sorted((V[b], perm[b]) for b in range(4)))
                    out[key] = out.get(key, 0) + perm_sign(perm)
                bm[V] = out
            return bm[V]

        def to_mon(pairs):
            m = [[0] * H for _ in th.outputs]
            for (v, a) in pairs:
                m[idx[v]][a] += 1
            return tuple(tuple(r) for r in m)
        W0 = {}
        for (V, Vp), cf in C0.items():
            for k1, s1 in bmon(V).items():
                for k2, s2 in bmon(Vp).items():
                    m = to_mon(k1 + k2)
                    W0[m] = W0.get(m, 0) + cf * s1 * s2
        W0 = {m: c for m, c in W0.items() if c}
        W1 = {}
        for j in range(3):
            for b in range(H):
                acc = {}
                for ((V, Vp), w_, jj), cf in C1.items():
                    if jj != j:
                        continue
                    for k1, s1 in bmon(V).items():
                        for k2, s2 in bmon(Vp).items():
                            m = to_mon(k1 + k2 + ((w_, b),))
                            acc[m] = acc.get(m, 0) + cf * s1 * s2
                W1[(j, b)] = {m: c for m, c in acc.items() if c}
        log('(v) monomial weights: K(lam;0) %d monomials; first order (j,b): %s monomials (%.1fs)'
            % (len(W0), [len(W1[(j, b)]) for j in range(3) for b in range(H)], time.time() - tw))
        weights[pqr] = dict(outputs=th.outputs, arms=[list(p_[1]) for p_ in th.paths], W0=W0, W1=W1, C0=C0, C1=C1, sgn=sgn)
        log('theta%s done (%.1fs)' % (pqr, time.time() - t0))
    pickle.dump(weights, open((__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/explore/msym-D8/verify/theta_weights.pkl', 'wb'))
    log('saved theta_weights.pkl; total %.1fs' % (time.time() - t_all))


if __name__ == '__main__':
    main()
