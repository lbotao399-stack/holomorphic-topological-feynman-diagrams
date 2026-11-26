"""d_2 = [Q_2] on H(Q_0) with general W: sectors of fixed integer spin weight Jw (letter weights below)
and flavour charges; complex graded by Rw (integer R-weight, Q has Rw = rq).  Usage:
  python3 lift_general.py <theory> <Jwmax>      Log: liftg_<theory>.log"""
import sys, time
from fractions import Fraction as Fr
from ht3d import Engine, load_kernels, pstr
from cohom import Echelon, kernel_basis, poly_to_vec, vec_to_poly
from theories import *

name, Jmax = sys.argv[1], int(sys.argv[2])
from math import lcm
# theory, R-charges r_k of the chirals, flavour charge vectors of phi_k (psi_k opposite)
SPEC = {
    "N8": (lambda: su2_adj(3, W=W_N8()), [Fr(2, 3)] * 3, [(1, 0), (-1, 1), (0, -1)]),
    "SQED1": (lambda: u1([1, -1, 0], W={((1, 0), (2, 0), (0, 0)): 1}), [Fr(1, 2), Fr(1, 2), Fr(1)], [(1,), (1,), (-2,)]),
    "U1": (lambda: u1([1]), [Fr(1, 2)], [(1,)]),
    "U1k": (lambda: u1([1], kappa=Fr(1, 2)), [Fr(1, 2)], [(1,)]),
    "U1q": (lambda: u1([1, 2, -3]), [Fr(1, 2)] * 3, [(1, 0, 0), (0, 1, 0), (0, 0, 1)]),
    "su2adj": (lambda: su2_adj(1), [Fr(1, 2)], [(1,)]),
    "su2adjk": (lambda: su2_adj(1, kappa=Fr(1, 3)), [Fr(1, 2)], [(1,)]),
    "gl2f": (lambda: gl_theory(2, 1, 0), [Fr(1, 2)], [(1,)]),
    "U11": (lambda: u1([1, 1]), [Fr(1, 2)] * 2, [(1, 0), (0, 1)]),
    "U1q2": (lambda: u1([2]), [Fr(1, 2)], [(1,)]),
    "U1q2k": (lambda: u1([2], kappa=1), [Fr(1, 2)], [(1,)]),
    "U1pm": (lambda: u1([1, -1]), [Fr(1, 2)] * 2, [(1, 0), (0, 1)]),
    "U1pmk": (lambda: u1([1, -1], kappa=1), [Fr(1, 2)] * 2, [(1, 0), (0, 1)]),
    "U11k": (lambda: u1([1, 1], kappa=1), [Fr(1, 2)] * 2, [(1, 0), (0, 1)]),
}
mk, rs, fl = SPEC[name]
th = mk()
S = 1
for r in rs:
    S = lcm(S, (r / 2).denominator, (1 - r / 2).denominator, r.denominator)
rq = S  # Q raises R by 1 -> Rw by S
def jw(L):
    n, i, d = L
    if n == 'c': return S * d
    if n == 'B': return S + S * d
    k = int(n[1:]); r = rs[k]
    return int(S * (r / 2 if n[0] == 'p' else 1 - r / 2)) + S * d
def rw(L):
    n = L[0]
    if n == 'c': return S
    if n == 'B': return 0
    k = int(n[1:]); r = rs[k]
    return int(S * (r if n[0] == 'p' else 1 - r))
def charges(L):
    n = L[0]
    if n in ('c', 'B'):
        return tuple(0 for _ in fl[0])
    k = int(n[1:]); s = 1 if n[0] == 'p' else -1
    return tuple(s * x for x in fl[k])
print("R-charges", [str(r) for r in rs], "weights scaled by", S)
E = Engine(th, load_kernels(4))
A = E.A
letters = []
for d in range(Jmax + 1):
    for a in range(th.G):
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
                ch = tuple(sum(x) for x in zip(*[charges(L) for L in m])) if m else None
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

print("theory", th.name, "W", "yes" if th.W else "no")
t0 = time.time()
for J in range(1, Jmax + 1):
    secs = enum(J)
    for ch in sorted(secs, key=lambda x: (x is None, x)):
        bas = secs[ch]
        idx = {R: {m: i for i, m in enumerate(b)} for R, b in bas.items()}
        Z, Bim = {}, {}
        for R, b in bas.items():
            imgs = []
            for m in b:
                P = E.Q0({m: Fr(1)})
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
        ranks = {}
        for R, h in Hn.items():
            if not H.get(R + rq):
                continue
            ech = Echelon(); base = Bim.get(R + rq, Echelon()); ech.rows = dict(base.rows); r0 = ech.rank()
            reps = []
            for zv in Z[R]:
                P = vec_to_poly(zv, bas[R])
                Q2P = E.Q2(P)
                if Q2P and ech.add(poly_to_vec(Q2P, idx[R + rq])):
                    reps.append((P, Q2P))
            if ech.rank() > r0:
                ranks[R] = ech.rank() - r0
                for P, Q2P in reps[:2]:
                    print("    lifted: O =", pstr(P, 6), "\n       Q2 O =", pstr(Q2P, 6))
        print("Jw=%d charges=%s: dim C = %s ; dim H(Q0) = %s%s" % (J, ch, {R: len(b) for R, b in sorted(bas.items())}, Hn,
              (" ; rank d2 (R -> R+%d) = %s" % (rq, ranks)) if ranks else ""))
        sys.stdout.flush()
print("time %.0fs" % (time.time() - t0))
