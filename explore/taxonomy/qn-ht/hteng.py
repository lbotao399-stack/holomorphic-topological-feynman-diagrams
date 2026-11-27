"""Unified (H,1) engine: letters, Q_0, Q_2 (all rooted thetas), Q_4 at z = 0 for H = 1 (fan(4),
tadpole(5,0)), for BF-type theories with cubic (and quadratic) underived interaction.

Letters: (name, idx, alpha), alpha an H-tuple (holomorphic derivatives d_1..d_H).
  'c'  c^a (odd); 'B' B_a (parity D); 'p<k>' phi of BF multiplet k (even); 's<k>' psi_k (parity D+1);
  'u<k>' u of symplectic multiplet k (even, D odd only).   D = H + 1.
Pairing on ordered pairs: eta(c^a,B_b) = s_g delta, eta(phi^i,psi_j) = s_m delta, eta(u^al,u^be) = Om^{al be};
  eta(y,x) = (-1)^{D + |x||y|} eta(x,y)   (edge reversal: kernel (-1)^D, Koszul (-1)^{|x||y|}).
Rule (Higher.tex Sec. 4.1, 3dN2 REPORT Sec. 2.1, general D): word I_0 I_1 ... I_n O; contracted pairs
(x_e at e(0), y_e at e(1)) moved to the front in edge order (Koszul), replaced by eta; survivors in order;
kernel K (htk conventions, same labelling, edge order, orientation): letter with derivative alpha at
e(1) -> y_e^alpha, lam_v -> d acting on the survivors at v.  Implemented as sequential contraction:
vertices added in decreasing label (left multiplication), edges (v,w), v<w, contracted in order (-v,w).
Q_n^hat O = sum_{rooted shapes} (1/|Aut_O|) [contractions] K/pi^{n/2};  Q_n = N_n Q_n^hat."""
import sys, os
from fractions import Fraction as Fr
from itertools import combinations, product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thetak


# ------------------------------------------------------------------ polynomial helpers
def padd(P, m, c):
    if c == 0:
        return
    v = P.get(m, 0) + c
    if v == 0:
        P.pop(m, None)
    else:
        P[m] = v


def psum(*Ps):
    out = {}
    for P in Ps:
        for m, c in P.items():
            padd(out, m, c)
    return out


def pscale(P, c):
    return {m: v * c for m, v in P.items()} if c else {}


def lname(L):
    n, i, a = L
    d = "".join("d%d" % (k + 1) * a[k] if a[k] <= 2 else "d%d^%d" % (k + 1, a[k]) for k in range(len(a)))
    return "%s%s_%s" % (d, n, i)


def pstr(P, maxterms=40):
    if not P:
        return "0"
    items = sorted(P.items(), key=lambda kv: (len(kv[0]), kv[0]))
    s = []
    for m, c in items[:maxterms]:
        s.append("%s*%s" % (c, " ".join(lname(L) for L in m) or "1"))
    if len(items) > maxterms:
        s.append("... (%d terms)" % len(items))
    return " + ".join(s)


def multi_add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def leibniz(letters, beta):
    """d^beta (multi-index) of a product of letters: list of (multinomial, new letters)"""
    if not any(beta):
        return [(1, tuple(letters))]
    if not letters:
        return []
    k = len(letters)
    # distribute each component independently
    per = []
    for i, b in enumerate(beta):
        opts = []

        def rec(j, rem, cur, mult):
            if j == k - 1:
                opts.append((mult, cur + [rem]))
                return
            for t in range(rem + 1):
                rec(j + 1, rem - t, cur + [t], mult * comb(rem, t))
        rec(0, b, [], 1)
        per.append(opts)
    out = []
    for choice in product(*per):
        mult = 1
        new = []
        for (m, _) in choice:
            mult *= m
        for j, L in enumerate(letters):
            a = tuple(L[2][i] + choice[i][1][j] for i in range(len(beta)))
            new.append((L[0], L[1], a))
        out.append((mult, tuple(new)))
    return out


# ------------------------------------------------------------------ theory
class Theory:
    """H; gauge: G, f {(a,b,c): v} ([t_a,t_b] = f_ab^c t_c);
    bf: list of (dimV, rho) with rho[a] = {(i,j): v} (rho(t_a)^i_j); W {tuple((k,i)...): v} on phi;
    sym: list of (dimU, Om, rho) with Om {(al,be): v} = Om^{al be} (eta of u), rho[a] = {(al,be): v};
    I = sB (1/2) f_bc^a B_a c^b c^c + sM psi_i c^a rho_a^i_j phi^j - sW W(phi)
        + sU (1/2) c^a om_{al ga} rho_a^ga_be u^al u^be,  om = (Om^{-1}) lower (om_{al be} Om^{be ga} = delta)."""

    def __init__(self, H, G, f, bf=(), W=None, sym=(), s_g=1, s_m=1, signs=None, name="",
                 omlow=None):
        self.H, self.D = H, H + 1
        if signs is None:
            # Q_0^2 = 0 (test_q0.log): sB = 1; sM = +1 (D even), -1 (D odd); sU = -1 (D odd)
            signs = (1, 1 if (H + 1) % 2 == 0 else -1, 1, -1)
        self.G, self.f = G, dict(f)
        self.bf, self.sym = list(bf), list(sym)
        self.W = W or {}
        self.s_g, self.s_m = s_g, s_m
        self.signs = signs
        self.name = name
        self.z = (0,) * H
        self.omlow = omlow
        self.I = self.interaction()
        self.I_list = list(self.I.items())
        self._partner_index()

    # parity
    def parity(self, L):
        n = L[0]
        if n == 'c':
            return 1
        if n == 'B':
            return self.D % 2
        if n[0] == 'p' or n[0] == 'u':
            return 0
        if n[0] == 's':
            return (self.D + 1) % 2
        raise ValueError(n)

    def canon(self, letters):
        L = list(letters)
        sign = 1
        par = self.parity
        for i in range(1, len(L)):
            j = i
            while j > 0 and L[j - 1] > L[j]:
                if par(L[j - 1]) and par(L[j]):
                    sign = -sign
                L[j - 1], L[j] = L[j], L[j - 1]
                j -= 1
        for i in range(1, len(L)):
            if L[i] == L[i - 1] and par(L[i]):
                return 0, None
        return sign, tuple(L)

    def mono(self, *letters, coeff=1):
        s, m = self.canon(letters)
        return {m: Fr(coeff) * s} if s else {}

    def L(self, n, i, *alpha):
        a = tuple(alpha) + (0,) * (self.H - len(alpha))
        return (n, i, a)

    def interaction(self):
        I = {}
        z = self.z
        sB, sM, sW, sU = self.signs
        for (b, c, a), v in self.f.items():
            s, m = self.canon([('B', a, z), ('c', b, z), ('c', c, z)])
            if s:
                padd(I, m, s * Fr(v) / 2 * sB)
        for k, (dimV, rho) in enumerate(self.bf):
            for a in range(self.G):
                for (i, j), v in rho[a].items():
                    s, m = self.canon([('s%d' % k, i, z), ('c', a, z), ('p%d' % k, j, z)])
                    if s:
                        padd(I, m, s * Fr(v) * sM)
        for key, v in self.W.items():
            s, m = self.canon([('p%d' % k, i, z) for (k, i) in key])
            if s:
                padd(I, m, -s * Fr(v) * sW)
        for k, (dimU, Om, rho) in enumerate(self.sym):
            om = self.omlow[k]
            for a in range(self.G):
                for (ga, be), v in rho[a].items():
                    for al in range(dimU):
                        w = om.get((al, ga), 0)
                        if w:
                            s, m = self.canon([('c', a, z), ('u%d' % k, al, z), ('u%d' % k, be, z)])
                            if s:
                                padd(I, m, s * Fr(v) * w / 2 * sU)
        return {m: c for m, c in I.items() if c}

    def eta(self, x, y):
        nx, ny = x[0], y[0]
        if nx == 'c' and ny == 'B':
            return self.s_g if x[1] == y[1] else 0
        if nx[0] == 'p' and ny[0] == 's' and nx[1:] == ny[1:]:
            return self.s_m if x[1] == y[1] else 0
        if nx[0] == 'u' and ny[0] == 'u' and nx == ny:
            return self.sym[int(nx[1:])][1].get((x[1], y[1]), 0)
        # reversed pairs
        if (nx == 'B' and ny == 'c') or (nx[0] == 's' and ny[0] == 'p'):
            v = self.eta(y, x)
            if v == 0:
                return 0
            sg = (-1) ** (self.D + self.parity(x) * self.parity(y))
            return sg * v
        return 0

    def _partner_index(self):
        """for each undifferentiated field key (name, idx): I-monomials containing an eta-partner"""
        self.partners = {}

    def letter_pool(self):
        """undifferentiated letters"""
        z = self.z
        out = [('c', a, z) for a in range(self.G)] + [('B', a, z) for a in range(self.G)]
        for k, (dimV, rho) in enumerate(self.bf):
            out += [('p%d' % k, i, z) for i in range(dimV)] + [('s%d' % k, i, z) for i in range(dimV)]
        for k, sp_ in enumerate(self.sym):
            out += [('u%d' % k, i, z) for i in range(sp_[0])]
        return out


# ------------------------------------------------------------------ shapes with kernels
class Shape:
    """labelled rooted graph: nv, edges (v,w) v<w in contraction order (-v,w), root nv-1, aut,
    moment(key) -> {beta: Fraction} (K/pi^{n/2}), beta over survivor vertices `surv`."""

    def __init__(self, nv, edges, aut, surv, moment, name):
        self.nv, self.edges, self.aut, self.surv, self.moment_fn, self.name = nv, edges, aut, surv, moment, name
        self.root = nv - 1
        assert all(v < w for v, w in edges)
        assert edges == sorted(edges, key=lambda e: (-e[0], e[1]))
        self.root_edges = [e for e, (v, w) in enumerate(edges) if w == self.root]
        self.deg = [0] * nv
        for v, w in edges:
            self.deg[v] += 1
            self.deg[w] += 1
        self.by_v = {}
        for e, (v, w) in enumerate(edges):
            self.by_v.setdefault(v, []).append((e, w))
        self.sidx = {v: i for i, v in enumerate(surv)}
        self.mcache = {}

    def moment(self, key):
        if key not in self.mcache:
            self.mcache[key] = self.moment_fn(key)
        return self.mcache[key]


def theta_shapes(H, seed=1, order_seed=None):
    out = []
    for k, (arms, rs, aut) in enumerate(thetak.rooted_shapes(H)):
        T = thetak.RootedTheta(H, arms, rs, order_seed=None if order_seed is None else order_seed + 17 * k)
        T.fix_sign(seed=seed, check_moments=1)
        sh = Shape(T.nv, T.edges, aut, T.surv, T.moment, "theta%s@%s" % (arms, rs))
        sh.theta = T
        out.append(sh)
    return out


def perm_sign(seq_from, seq_to):
    """sign of the permutation taking list seq_from to seq_to"""
    pos = [seq_from.index(x) for x in seq_to]
    s = 1
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            if pos[i] > pos[j]:
                s = -s
    return s


def h1_q4_shapes():
    """H = 1, z = 0, cubic theories: fan(4) K = (pi^2/180)(l0 l4)^2 and tadpole(5,0) K = (pi^2/630) l0^4
    in the edge orders of explore/physics/3dN2-supercharge/fanengine.py (|Aut| = 2 each), converted to
    the order (-v,w): K -> sgn(pi)^{D+1} K, D = 2."""
    shapes = []
    # fan(4): fanengine order
    R = 5
    fe = []
    for i in range(4):
        fe += [(i, i + 1), (i, R)]
    fe.append((4, R))
    mine = sorted(fe, key=lambda e: (-e[0], e[1]))
    sg = perm_sign(fe, mine) ** 3
    surv = [0, 4]

    def mfan(key, sg=sg):
        if any(any(p) for p in key):
            raise ValueError("Q4 only at z = 0")
        return {((2,), (2,)): Fr(1, 180) * sg}
    shapes.append(Shape(6, mine, 2, surv, mfan, "fan(4)"))
    # tadpole(5,0)
    te = sorted([(min(a, b), max(a, b)) for a, b in [(i, (i + 1) % 5) for i in range(5)]] + [(v, R) for v in range(1, 5)])
    mine = sorted(te, key=lambda e: (-e[0], e[1]))
    sg2 = perm_sign(te, mine) ** 3

    def mtad(key, sg=sg2):
        if any(any(p) for p in key):
            raise ValueError("Q4 only at z = 0")
        return {((4,),): Fr(1, 630) * sg}
    shapes.append(Shape(6, mine, 2, [0], mtad, "tadpole(5,0)"))
    return shapes


# ------------------------------------------------------------------ engine
class Engine:
    def __init__(self, th, shapes=None, seed=1):
        self.th = th
        self.H = th.H
        self.shapes = shapes if shapes is not None else theta_shapes(th.H, seed)
        self.cache0, self.cacheF, self.cacheQ = {}, {}, {}
        par = th.parity
        self.Itag = [(m, c) for m, c in th.I_list]
        # partner index: undiff field (name, idx) -> list of (I-monomial, coeff)
        self.Iby = {}
        for m, c in th.I_list:
            for x in m:
                self.Iby.setdefault((x[0], x[1]), []).append((m, c))
        pool = th.letter_pool()
        self.partners_of = {}
        for y in pool:
            ps = [x for x in pool if th.eta(x, y) != 0]
            self.partners_of[(y[0], y[1])] = ps

    # ---------------- Q_0
    def Q0_mono(self, O):
        if O in self.cache0:
            return self.cache0[O]
        th, par = self.th, self.th.parity
        out = {}
        for mI, cI in th.I_list:
            for i, x in enumerate(mI):
                for j, y in enumerate(O):
                    et = th.eta(x, y)
                    if not et:
                        continue
                    # word = mI + O ; move x (pos i) then y (pos len(mI)+j) to the front
                    s1 = (-1) ** (par(x) * sum(par(t) for t in mI[:i]))
                    rest = list(mI[:i]) + list(mI[i + 1:])
                    s2 = (-1) ** (par(y) * (sum(par(t) for t in rest) + sum(par(t) for t in O[:j])))
                    Irest = rest
                    Orest = list(O[:j]) + list(O[j + 1:])
                    base = cI * et * s1 * s2
                    for mult, lets in leibniz(Irest, y[2]):
                        s, m = th.canon(list(lets) + Orest)
                        if s:
                            padd(out, m, base * mult * s)
        self.cache0[O] = out
        return out

    def Q0(self, P):
        out = {}
        for O, c in P.items():
            for m, v in self.Q0_mono(O).items():
                padd(out, m, c * v)
        return out

    # ---------------- vertex function F_shape(root word)
    def F(self, si, rw):
        key = (si, rw)
        if key in self.cacheF:
            return self.cacheF[key]
        sh = self.shapes[si]
        th, par = self.th, self.th.parity
        root = sh.root
        # state: {(word, mkey): coeff}; word = tuple of tagged letters (v, letter)
        state = {(tuple((root, L) for L in rw), ()): Fr(1)}
        remaining_edges = list(sh.deg)  # uncontracted edges per vertex
        for v in range(sh.nv - 2, -1, -1):
            elist = sh.by_v.get(v, [])
            new = {}
            # candidate I monomials: must supply partners for the first edge's letter
            for (word, mk), cf in state.items():
                # letters available at each neighbour
                e0, w0 = elist[0]
                ys0 = [L for (t, L) in word if t == w0]
                cand = {}
                for y in ys0:
                    for x in self.partners_of.get((y[0], y[1]), []):
                        for (m, c) in self.Iby.get((x[0], x[1]), []):
                            cand[m] = c
                for m, c in cand.items():
                    if len(m) < len(elist):
                        continue
                    tw = tuple((v, L) for L in m) + word
                    terms = {(tw, mk): cf * c}
                    for (e, w) in elist:
                        nt = {}
                        for (wd, k2), c2 in terms.items():
                            n = len(wd)
                            for i in range(n):
                                if wd[i][0] != v:
                                    continue
                                x = wd[i][1]
                                for j in range(n):
                                    if wd[j][0] != w:
                                        continue
                                    y = wd[j][1]
                                    et = th.eta(x, y)
                                    if not et:
                                        continue
                                    s1 = sum(par(wd[t][1]) for t in range(i)) * par(x)
                                    rest = wd[:i] + wd[i + 1:]
                                    jj = j if j < i else j - 1
                                    s2 = sum(par(rest[t][1]) for t in range(jj)) * par(y)
                                    rest2 = rest[:jj] + rest[jj + 1:]
                                    sg = -1 if (s1 + s2) % 2 else 1
                                    k3 = k2 + (y[2],) if w == root else k2
                                    padd(nt, (rest2, k3), c2 * et * sg)
                        terms = nt
                        if not terms:
                            break
                    for (wd, k2), c2 in terms.items():
                        # canonical sort of tagged word (Koszul)
                        s, cw = self._canon_tagged(wd)
                        if s:
                            padd(new, (cw, k2), c2 * s)
            state = new
            if not state:
                break
            # prune: vertex w fully processed when all its edges contracted -> survivors fixed
        out = {}
        for (word, mk), cf in state.items():
            if any(t == root for (t, L) in word):
                continue
            mom = sh.moment(mk)
            groups = {}
            for (t, L) in word:
                groups.setdefault(t, []).append(L)
            # every survivor vertex must be in sh.surv
            if any(t not in sh.sidx for t in groups):
                # survivors at a vertex without lambda (should not happen for thetas)
                raise ValueError("survivor at non-survivor vertex %s in %s" % (groups, sh.name))
            order = [t for (t, L) in word]
            for beta, kc in mom.items():
                parts = []
                ok = True
                for t in sh.surv:
                    b = beta[sh.sidx[t]]
                    g = groups.get(t, [])
                    if not g:
                        if any(b):
                            ok = False
                            break
                        parts.append([(1, ())])
                        continue
                    parts.append(leibniz(g, b))
                if not ok:
                    continue
                # survivors keep word order: word is canonical-tagged (sorted by tag), groups in tag order
                tags_sorted = sorted(groups)
                for combo in product(*parts):
                    mult = 1
                    lets = []
                    bytag = {}
                    for t, (ml, ls) in zip(sh.surv, combo):
                        mult *= ml
                        bytag[t] = ls
                    for t in tags_sorted:
                        lets += list(bytag[t])
                    s, mm = th.canon(lets)
                    if s:
                        padd(out, mm, cf * kc * mult * s)
        self.cacheF[key] = out
        return out

    def _canon_tagged(self, wd):
        par = self.th.parity
        L = list(wd)
        sign = 1
        for i in range(1, len(L)):
            j = i
            while j > 0 and L[j - 1] > L[j]:
                if par(L[j - 1][1]) and par(L[j][1]):
                    sign = -sign
                L[j - 1], L[j] = L[j], L[j - 1]
                j -= 1
        for i in range(1, len(L)):
            if L[i] == L[i - 1] and par(L[i][1]):
                return 0, None
        return sign, tuple(L)

    # ---------------- Q_shape on monomials
    def Qshape_mono(self, si, O):
        sh = self.shapes[si]
        th, par = self.th, self.th.parity
        d = sh.deg[sh.root]
        out = {}
        n = len(O)
        for S in combinations(range(n), d):
            rw = tuple(O[i] for i in S)
            Rw = [O[i] for i in range(n) if i not in S]
            # Koszul sign of moving S to the front preserving order
            sg = 0
            for k, i in enumerate(S):
                sg += par(O[i]) * sum(par(O[t]) for t in range(i) if t not in S)
            sg = -1 if sg % 2 else 1
            Fv = self.F(si, rw)
            for m, c in Fv.items():
                s, mm = th.canon(list(m) + Rw)
                if s:
                    padd(out, mm, c * sg * s)
        return out

    def Q2_mono(self, O, which=None):
        which = tuple(which) if which is not None else None
        key = (O, which)
        if key in self.cacheQ:
            return self.cacheQ[key]
        out = {}
        for si, sh in enumerate(self.shapes):
            if which is not None and si not in which:
                continue
            q = self.Qshape_mono(si, O)
            for m, c in q.items():
                padd(out, m, c / sh.aut)
        self.cacheQ[key] = out
        return out

    def Q2(self, P, which=None):
        out = {}
        for O, c in P.items():
            for m, v in self.Q2_mono(O, which).items():
                padd(out, m, c * v)
        return out
