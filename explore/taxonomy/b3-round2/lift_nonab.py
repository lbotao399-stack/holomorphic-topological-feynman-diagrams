"""d_2 = [Q_2^hat] on H(Q_0) for GL(2) theories with chiral (anti)fundamentals, 3d N=2 HT twist (H=1).
Engine: explore/physics/3dN2-supercharge/ht3d.py (imported, not edited).  Sectors: fixed Jw (spin weight x4),
gauge Cartan weight (0,0) (H(Q_0) is g-invariant), complex graded by Rw.  Optional CS form
I_CS = (1/2) k_ab c^a dc^b with k = s1 tr_F(XY) + s2 trX trY.
Usage: python3 lift_nonab.py <theory> <Jwmax> [s1 s2]    Log: lift_<theory>[_k<s1>_<s2>]_<Jwmax>.log
theories: gl2f (1F), gl2ff (2F), gl2fa (1F+1AF), u1 (U(1)+chiral q=1, control)."""
import sys, os, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "physics", "3dN2-supercharge")
sys.path.insert(0, SRC)
from ht3d import Engine, load_kernels, pstr, padd, pclean
from cohom import Echelon, kernel_basis, poly_to_vec, vec_to_poly
import theories as TH3

name, Jmax = sys.argv[1], int(sys.argv[2])
s1 = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0)
s2 = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
tag = name + ("_k%s_%s" % (str(s1).replace('/', 'o'), str(s2).replace('/', 'o')) if (s1 or s2) else "")
log = open(os.path.join(HERE, "lift_%s_%d.log" % (tag, Jmax)), "w")
def out(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n"); log.flush()

BAS = [(0, 0), (0, 1), (1, 0), (1, 1)]  # gl(2) basis E_ij, index 2i+j
if name == "gl2f":
    th = TH3.gl_theory(2, 1, 0); kinds = ['F']
elif name == "gl2ff":
    th = TH3.gl_theory(2, 2, 0); kinds = ['F', 'F']
elif name == "gl2fa":
    th = TH3.gl_theory(2, 1, 1); kinds = ['F', 'A']
elif name == "u1":
    th = TH3.u1([1]); kinds = ['q1']
else:
    raise SystemExit("unknown theory")
G = th.G
# CS form
kf = {}
if G == 4:
    for a, (i, j) in enumerate(BAS):
        for b, (k, l) in enumerate(BAS):
            v = s1 * (1 if (j == k and l == i) else 0) + s2 * (1 if (i == j and k == l) else 0)
            if v:
                kf[(a, b)] = v
elif s1 or s2:
    kf[(0, 0)] = s1 + s2
if kf:
    A = th.A
    for (a, b), v in kf.items():
        s, m = A.canon([('c', a, 0), ('c', b, 1)])
        if s:
            padd(th.I, m, s * v / 2)
    th.I = pclean(th.I)
    th.I_list = list(th.I.items())
rs = [Fr(1, 2)] * len(th.multiplets)
S = 4
rq = S
def jw(L):
    n, i, d = L
    if n == 'c': return S * d
    if n == 'B': return S + S * d
    return (1 if n[0] == 'p' else 3) + S * d
def rw(L):
    n = L[0]
    if n == 'c': return S
    if n == 'B': return 0
    return 2
nfl = len(kinds)
def charges(L):
    """(Cartan weight h_1, h_2 or U(1) charge, flavour numbers per multiplet)"""
    n, i, _ = L
    if G == 4:
        if n in ('c', 'B'):
            a, b = BAS[i]
            w = [(1 if a == 0 else 0) - (1 if b == 0 else 0), (1 if a == 1 else 0) - (1 if b == 1 else 0)]
            if n == 'B':
                w = [-x for x in w]
            return tuple(w) + (0,) * nfl
        k = int(n[1:])
        w = [1 if i == 0 else 0, 1 if i == 1 else 0]
        if kinds[k] == 'A':
            w = [-x for x in w]
        fl = [1 if j == k else 0 for j in range(nfl)]
        t = tuple(w) + tuple(fl)
        return t if n[0] == 'p' else tuple(-x for x in t)
    if n in ('c', 'B'):
        return (0, 0)
    return (1, 1) if n[0] == 'p' else (-1, -1)
ng = 2 if G == 4 else 1
def gauge_zero(ch):
    return all(x == 0 for x in ch[:ng])

out("theory", th.name, "; CS form k = %s tr_F(XY) + %s trX trY" % (s1, s2) if G == 4 else "; kappa = %s" % (s1 + s2))
out("I =", pstr(th.I, 80))
E = Engine(th, load_kernels(4))
A = E.A
letters = []
for d in range(Jmax + 1):
    for a in range(G):
        letters += [('c', a, d), ('B', a, d)]
    for k, m in enumerate(th.multiplets):
        for i in range(m[0]):
            letters += [('p%d' % k, i, d), ('s%d' % k, i, d)]
letters = sorted(L for L in letters if jw(L) <= Jmax)

def enum(J):
    res = {}
    def rec(start, cur, j):
        if j == J:
            s, m = A.canon(cur)
            if s:
                ch = tuple(sum(x) for x in zip(*[charges(L) for L in m]))
                if gauge_zero(ch):
                    R = sum(rw(L) for L in m)
                    res.setdefault(ch, {}).setdefault(R, set()).add(m)
        for idx in range(start, len(letters)):
            L = letters[idx]
            w = jw(L)
            if j + w > J or (w == 0 and L[0] != 'c'):
                continue
            if A.parity(L) and cur and cur[-1] == L:
                continue
            rec(idx, cur + [L], j + w)
    rec(0, [], 0)
    return {ch: {R: sorted(v) for R, v in d.items()} for ch, d in res.items()}

t0 = time.time()
nq0sq = 0
tot = {}
for J in range(1, Jmax + 1):
    secs = enum(J)
    for ch in sorted(secs):
        bas = secs[ch]
        idx = {R: {m: i for i, m in enumerate(b)} for R, b in bas.items()}
        Z, Bim = {}, {}
        for R, b in bas.items():
            imgs = []
            for m in b:
                P = E.Q0({m: Fr(1)})
                if P:
                    PP = E.Q0(P)
                    if PP:
                        nq0sq += 1
                imgs.append(poly_to_vec(P, idx.get(R + rq, {})) if P else {})
            Z[R] = kernel_basis(None, imgs, None)
            ech = Echelon()
            for v in imgs:
                if v:
                    ech.add(v)
            Bim[R + rq] = ech
        H = {R: len(Z[R]) - (Bim[R].rank() if R in Bim else 0) for R in bas}
        Hn = {R: h for R, h in sorted(H.items()) if h}
        if not Hn:
            continue
        ranks, nz = {}, 0
        for R, h in Hn.items():
            if not H.get(R + rq):
                continue
            ech = Echelon(); base = Bim.get(R + rq, Echelon()); ech.rows = dict(base.rows); r0 = ech.rank()
            reps = []
            for zv in Z[R]:
                P = vec_to_poly(zv, bas[R])
                Q2P = E.Q2(P)
                if Q2P:
                    nz += 1
                    if ech.add(poly_to_vec(Q2P, idx[R + rq])):
                        reps.append((P, Q2P))
            if ech.rank() > r0:
                ranks[R] = ech.rank() - r0
                for P, Q2P in reps[:2]:
                    out("    lifted: O =", pstr(P, 8), "\n       Q2^ O =", pstr(Q2P, 8))
        tot.setdefault(J, []).append((ch, Hn, ranks))
        out("Jw=%d charges=%s: dim C = %s ; dim H(Q0) = %s ; Q2^z!=0 on %d cocycles%s" % (
            J, ch, {R: len(b) for R, b in sorted(bas.items())}, Hn, nz,
            (" ; rank d2 (R -> R+%d) = %s" % (rq, ranks)) if ranks else ""))
out("Q0^2 != 0 on %d monomials" % nq0sq)
out("summary (Jw: H by Rw summed over flavour; rank d2):")
for J, lst in sorted(tot.items()):
    Hs, rk = {}, {}
    for ch, Hn, ranks in lst:
        for R, h in Hn.items():
            Hs[R] = Hs.get(R, 0) + h
        for R, r in ranks.items():
            rk[R] = rk.get(R, 0) + r
    out("  Jw=%d: H = %s ; rank d2 = %s" % (J, dict(sorted(Hs.items())), dict(sorted(rk.items()))))
out("time %.0fs" % (time.time() - t0))
