"""3d N=2 HT-twisted gauge theory: letters, exact polynomials, Q_0 and Q_2 from the GKW graph rule.

Letters (name, idx, d): 'c' ghost c^a (odd), 'B' B_a (even), 'p<k>' phi of multiplet k, 's<k>' psi of
multiplet k; d = number of holomorphic derivatives.  Parity of phi^k = par[k] (0 default), psi = 1+par.
Polynomials: dict {sorted letter tuple: Fraction}.

Graph rule (Higher.tex Sec. 4.1, 4.4): for a labelled 2-Laman graph Gamma on vertices
(I_0,...,I_n, O = root), pairs (x_e at e(0), y_e at e(1)) are moved to the front of the ordered product
I_0 I_1 ... I_n O in the edge order of Gamma (Koszul signs), replaced by eta(x_e,y_e), and the survivors
are multiplied in their original order.  The kernel K_Gamma (htk conventions, same edge order and
orientation) acts by: letter with d derivatives at e(0): (-y_e)^d, at e(1): (+y_e)^d; lambda_v -> d/dz
acting on the survivors at v.  Q_n O = (1/(n+1)!) sum over labelled contractions, times N_n.
"""
from fractions import Fraction as Fr
from itertools import product, permutations
from math import factorial, comb
import pickle, os

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- letters and polynomials
class Algebra:
    def __init__(self, par_phi):
        self.par_phi = list(par_phi)

    def parity(self, L):
        n = L[0]
        if n == 'c':
            return 1
        if n == 'B':
            return 0
        k = int(n[1:])
        return (self.par_phi[k] + (1 if n[0] == 's' else 0)) % 2

    def canon(self, letters):
        """sort letters (Koszul sign); returns (sign, tuple) or (0, None)"""
        L = list(letters)
        sign = 1
        # insertion sort with sign
        for i in range(1, len(L)):
            j = i
            while j > 0 and L[j - 1] > L[j]:
                if self.parity(L[j - 1]) and self.parity(L[j]):
                    sign = -sign
                L[j - 1], L[j] = L[j], L[j - 1]
                j -= 1
        for i in range(1, len(L)):
            if L[i] == L[i - 1] and self.parity(L[i]):
                return 0, None
        return sign, tuple(L)

    def mono_parity(self, m):
        return sum(self.parity(L) for L in m) % 2

    def perm_sign(self, letters, order):
        """Koszul sign of rearranging letters into [letters[i] for i in order]"""
        par = [self.parity(letters[i]) for i in order]
        s = 1
        for i in range(len(order)):
            if not par[i]:
                continue
            for j in range(i + 1, len(order)):
                if par[j] and order[i] > order[j]:
                    s = -s
        return s


def padd(P, m, c):
    if c == 0:
        return
    v = P.get(m, 0) + c
    if v == 0:
        P.pop(m, None)
    else:
        P[m] = v


def pclean(P):
    return {m: c for m, c in P.items() if c != 0}


def pscale(P, c):
    return {m: v * c for m, v in P.items()} if c != 0 else {}


def psum(*Ps):
    out = {}
    for P in Ps:
        for m, c in P.items():
            padd(out, m, c)
    return out


def pmul(A, P, Q):
    out = {}
    for m1, c1 in P.items():
        for m2, c2 in Q.items():
            s, m = A.canon(m1 + m2)
            if s:
                padd(out, m, s * c1 * c2)
    return out


def deriv_dist(letters, a):
    """d^a applied to the product of letters (Leibniz): list of (multinomial, new letters)"""
    if a == 0:
        return [(1, tuple(letters))]
    if not letters:
        return []
    out = []
    k = len(letters)

    def rec(i, rem, cur, mult):
        if i == k - 1:
            L = letters[i]
            out.append((mult, tuple(cur + [(L[0], L[1], L[2] + rem)])))
            return
        for t in range(rem + 1):
            L = letters[i]
            rec(i + 1, rem - t, cur + [(L[0], L[1], L[2] + t)], mult * comb(rem, t))

    rec(0, a, [], 1)
    return out


def pstr(P, maxterms=60):
    if not P:
        return "0"
    items = sorted(P.items(), key=lambda kv: (len(kv[0]), kv[0]))
    s = []
    for m, c in items[:maxterms]:
        s.append("%s*%s" % (c, "".join(lname(L) for L in m) or "1"))
    if len(items) > maxterms:
        s.append("... (%d terms)" % len(items))
    return " + ".join(s)


def lname(L):
    n, i, d = L
    return "[%s%s%s]" % ("d" * d if d <= 3 else "d^%d" % d, n, i)


# ---------------------------------------------------------------- theory
class Theory:
    """gauge algebra: dimension G, structure constants f[(a,b,c)] with [t_a,t_b] = f_ab^c t_c;
    multiplets: list of (dimV, rho) with rho[a] = {(i,j): value} = rho(t_a)^i_j, par = parity of phi;
    W: dict {tuple of (k,i) with k multiplet index: coeff}; kappa: CS coefficient,
    I_CS = (kappa/2) sum_a c^a dc^a; eta signs s_g = eta(c,B), s_m = eta(phi,psi)."""

    def __init__(self, G, f, multiplets, W=None, kappa=0, s_g=1, s_m=1, signs=(1, 1, 1), name=""):
        self.G, self.f, self.multiplets = G, dict(f), multiplets
        self.W = W or {}
        self.kappa = Fr(kappa)
        self.s_g, self.s_m = s_g, s_m
        self.signs = signs  # coefficients of (B c c, psi c phi, -W)
        self.name = name
        self.A = Algebra([m[2] if len(m) > 2 else 0 for m in multiplets])
        self.I = self.interaction()
        self.I_list = list(self.I.items())

    def interaction(self):
        A, I = self.A, {}
        sB, sM, sW = self.signs
        for (b, c, a), v in self.f.items():  # (1/2) f_bc^a B_a c^b c^c
            s, m = A.canon([('B', a, 0), ('c', b, 0), ('c', c, 0)])
            if s:
                padd(I, m, s * Fr(v) / 2 * sB)
        for k, mult in enumerate(self.multiplets):
            dimV, rho = mult[0], mult[1]
            for a in range(self.G):
                for (i, j), v in rho[a].items():  # psi_i c^a rho_a^i_j phi^j
                    s, m = A.canon([('s%d' % k, i, 0), ('c', a, 0), ('p%d' % k, j, 0)])
                    if s:
                        padd(I, m, s * Fr(v) * sM)
        for key, v in self.W.items():  # -W
            s, m = A.canon([('p%d' % k, i, 0) for (k, i) in key])
            if s:
                padd(I, m, -s * Fr(v) * sW)
        if self.kappa:
            for a in range(self.G):
                s, m = A.canon([('c', a, 0), ('c', a, 1)])
                padd(I, m, s * self.kappa / 2)
        return pclean(I)

    def eta(self, x, y):
        if x[1] != y[1]:
            return 0
        nx, ny = x[0], y[0]
        if {nx, ny} == {'c', 'B'}:
            return self.s_g
        if nx[0] in 'ps' and ny[0] in 'ps' and nx[0] != ny[0] and nx[1:] == ny[1:]:
            return self.s_m
        return 0


# ---------------------------------------------------------------- kernels
def load_kernels(M=4):
    import sympy as sp
    tab = pickle.load(open(os.path.join(HERE, "kernels_%d.pkl" % M), "rb"))
    l = sp.symbols("l0 l1 l2")
    z = sp.symbols("z0:5")
    out = {}
    for edges, Ks in tab.items():
        mom = {}
        for m, K in enumerate(Ks):
            K = sp.expand(K / sp.pi)
            Pz = sp.Poly(K, *z)
            for zmon, coeff in Pz.terms():
                sgn = (-1) ** sum(zmon)
                fac = 1
                for p in zmon:
                    fac *= factorial(p)
                Pl = sp.Poly(sp.expand(coeff), *l)
                d = {}
                for lm, c in Pl.terms():
                    d[lm] = Fr(int(sp.numer(c)), int(sp.denom(c))) * sgn * fac
                mom[zmon] = d
        out[edges] = (mom, M)
    return out


class Engine:
    def __init__(self, th, kernels=None, M=4):
        self.th, self.A = th, th.A
        self.K = kernels if kernels is not None else load_kernels(M)
        self.M = M
        self.cache0, self.cache2 = {}, {}
        self._ccache = {}

    # ----- one contraction term -> output poly
    def _term(self, mons, coeffs, edges, legpairs, vsurv_lambda, moment):
        """mons: list of letter tuples (vertices in order, O last); legpairs: list over edges of
        ((v,pos),(w,pos)) with v<w; moment: dict {(a_v...): Fraction} lambda polynomial."""
        A = self.A
        letters, offs = [], []
        for m in mons:
            offs.append(len(letters))
            letters += list(m)
        used = set()
        order = []
        eta = 1
        for (v, i), (w, j) in legpairs:
            gi, gj = offs[v] + i, offs[w] + j
            order += [gi, gj]
            used.add(gi), used.add(gj)
            eta *= self.th.eta(letters[gi], letters[gj])
        if eta == 0:
            return {}
        surv_idx = [g for g in range(len(letters)) if g not in used]
        sign = A.perm_sign(letters, order + surv_idx)
        base = sign * eta
        for c in coeffs:
            base *= c
        # survivors grouped by vertex
        nv = len(mons)
        groups = [[letters[offs[v] + i] for i in range(len(mons[v])) if offs[v] + i not in used]
                  for v in range(nv)]
        out = {}
        for lam_exp, kc in moment.items():
            parts = [[(1, ())]]
            ok = True
            for v in range(nv - 1):
                a = lam_exp[v]
                dd = deriv_dist(groups[v], a)
                if not dd:
                    ok = False
                    break
                parts.append(dd)
            if not ok:
                continue
            parts.append([(1, tuple(groups[nv - 1]))])
            for combo in product(*parts[1:]):
                mult = 1
                lets = []
                for mlt, ls in combo:
                    mult *= mlt
                    lets += list(ls)
                s, m = A.canon(lets)
                if s:
                    padd(out, m, base * kc * mult * s)
        return out

    # ----- Q_0
    def Q0_mono(self, O):
        if O in self.cache0:
            return self.cache0[O]
        out = {}
        for mI, cI in self.th.I_list:
            for i, x in enumerate(mI):
                for j, y in enumerate(O):
                    if self.th.eta(x, y) == 0:
                        continue
                    p = x[2] + y[2]
                    sgn = (-1) ** x[2]
                    moment = {(p,): Fr(sgn)}
                    t = self._term([mI, O], [cI], None, [((0, i), (1, j))], None, moment)
                    for m, c in t.items():
                        padd(out, m, c)
        self.cache0[O] = out
        return out

    def Q0(self, P):
        out = {}
        for O, c in P.items():
            for m, v in self.Q0_mono(O).items():
                padd(out, m, c * v)
        return out

    def _compat(self, m1, m2):
        key = (m1, m2)
        c = self._ccache.get(key)
        if c is None:
            c = any(self.th.eta(x, y) != 0 for x in m1 for y in m2)
            self._ccache[key] = c
        return c

    # ----- Q_2
    def Q2_mono(self, O, only=None):
        """only: None, 'fan' (O degree 3) or 'tad' (O degree 2)"""
        key = (O, only)
        if key in self.cache2:
            return self.cache2[key]
        th = self.th
        out = {}
        pairs_all = [(i, j) for i in range(4) for j in range(i + 1, 4)]
        Il = th.I_list
        for missing in pairs_all:
            if only == 'fan' and 3 in missing:
                continue
            if only == 'tad' and 3 not in missing:
                continue
            edges = tuple(p for p in pairs_all if p != missing)
            mom_tab, M = self.K[edges]
            eset = set(edges)
            compat = self._compat
            okO = [compat(m, O) for m, _ in Il]
            triples = []
            for a0, (m0, c0) in enumerate(Il):
                if (0, 3) in eset and not okO[a0]:
                    continue
                for a1, (m1, c1) in enumerate(Il):
                    if (1, 3) in eset and not okO[a1]:
                        continue
                    if (0, 1) in eset and not compat(m0, m1):
                        continue
                    for a2, (m2, c2) in enumerate(Il):
                        if (2, 3) in eset and not okO[a2]:
                            continue
                        if (0, 2) in eset and not compat(m0, m2):
                            continue
                        if (1, 2) in eset and not compat(m1, m2):
                            continue
                        triples.append(((m0, c0), (m1, c1), (m2, c2)))
            for (m0, c0), (m1, c1), (m2, c2) in triples:
                mons = [m0, m1, m2, O]
                # candidate leg pairs for each edge
                cands = []
                okall = True
                for (v, w) in edges:
                    cl = [(i, j) for i, x in enumerate(mons[v]) for j, y in enumerate(mons[w])
                          if th.eta(x, y) != 0]
                    if not cl:
                        okall = False
                        break
                    cands.append(cl)
                if not okall:
                    continue
                for choice in product(*cands):
                    usedl = set()
                    good = True
                    for (v, w), (i, j) in zip(edges, choice):
                        if (v, i) in usedl or (w, j) in usedl:
                            good = False
                            break
                        usedl.add((v, i)), usedl.add((w, j))
                    if not good:
                        continue
                    pexp, sgn = [], 1
                    for (v, w), (i, j) in zip(edges, choice):
                        dx, dy = mons[v][i][2], mons[w][j][2]
                        pexp.append(dx + dy)
                        sgn *= (-1) ** dx
                    pexp = tuple(pexp)
                    if sum(pexp) > M:
                        raise ValueError("kernel order too low: %s" % (pexp,))
                    mom = mom_tab.get(pexp)
                    if not mom:
                        continue
                    moment = {k: v * sgn for k, v in mom.items()}
                    legpairs = [((v, i), (w, j)) for (v, w), (i, j) in zip(edges, choice)]
                    t = self._term(mons, [c0, c1, c2], edges, legpairs, None, moment)
                    for m, c in t.items():
                        padd(out, m, c)
        out = {m: c / 6 for m, c in out.items() if c != 0}
        self.cache2[key] = out
        return out

    def Q2(self, P, only=None):
        out = {}
        for O, c in P.items():
            for m, v in self.Q2_mono(O, only).items():
                padd(out, m, c * v)
        return out


def mono(A, *letters, coeff=1):
    s, m = A.canon(letters)
    return {m: Fr(coeff) * s} if s else {}
