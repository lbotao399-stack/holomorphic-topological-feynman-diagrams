"""H^3, H^4 of the G-invariant CE complex of sl_N[[z1,z2]] (mod sum p) at weights (2,2), (3,3) [and (4,2), (5,1),
(6,0) for H^4], for N = 3, 4, 5 (faithful fixed-N coordinates, vfix) and N = 7, 11 (trace basis, faithful since
k <= 5 < N).  C-odd / C-even subcomplexes, SL2-invariant subcomplex.  Classes of theta', {D,iota}, o_3' =
{D,iota} + theta' (from v3, formal), and of a generator h of H^4_(3,3); C-parity of H^4.
All ranks exact over Q (sympy sparse rref over QQ), plus the same ranks mod p = 2147483629 printed as a cross-check."""
import pickle, time
from fractions import Fraction as Fr
from vlib import *
import vla
from vfix import Faithful

t0 = time.time()
SL = True
A = pickle.load(open('data/assembled_sl.pkl', 'rb'))
o3p = lin((1, A['th']['d']), (1, A['DI']))
thp = A['th']['d']
DI = A['DI']

NS = [3, 4, 5, 7, 11]
FA = {}
for N in NS:
    for k in (2, 3, 4, 5):
        mode = 'all' if (N == 3 and k <= 5) else 'random'
        FA[(N, k)] = Faithful(N, SL, k, mode=mode)
print("faithful ranks:", {key: F.rank for key, F in FA.items()}, flush=True)


def rows(cochains, N, k):
    F = FA[(N, k)]
    return [F.proj(w, N) for w in cochains]


def formal_basis(cs, k):
    """independent subset at N = 11 (trace basis faithful)."""
    idx = vla.indep_rows(rows(cs, 11, k))
    return [cs[i] for i in idx]


def complex_dims(deg, parity=None, inv=False, tag=""):
    res = {}
    B = {}
    for k in (2, 3, 4):
        if inv and k == 2:
            cs = span_set(k, deg, SL, sl2inv=True, parity=parity) if deg[0] == deg[1] else []
        elif inv:
            cs = span_set(k, deg, SL, sl2inv=True, parity=parity)
        else:
            cs = span_set(k, deg, SL, parity=parity)
        B[k] = formal_basis(cs, k) if cs else []
    dB = {k: [d(w, k, SL) for w in B[k]] for k in B}
    for N in NS:
        dimC = {k: (vla.rankQ(rows(B[k], N, k)) if B[k] else 0) for k in (2, 3, 4)}
        rk = {k: (vla.rankQ(rows(dB[k], N, k + 1)) if B[k] else 0) for k in (2, 3, 4)}
        rkp = {k: (vla.rankp(rows(dB[k], N, k + 1)) if B[k] else 0) for k in (2, 3, 4)}
        H3 = dimC[3] - rk[3] - rk[2]
        H4 = dimC[4] - rk[4] - rk[3]
        res[N] = (dimC, rk, H3, H4)
        print("  %s deg %s N=%d: dim C^2,C^3,C^4 = %d,%d,%d; rank d2,d3,d4 = %d,%d,%d (mod p: %d,%d,%d); H^3 = %d, H^4 = %d"
              % (tag, deg, N, dimC[2], dimC[3], dimC[4], rk[2], rk[3], rk[4], rkp[2], rkp[3], rkp[4], H3, H4), flush=True)
    return res, B, dB


print("== full complex ==")
full = {}
for deg in [(2, 2), (3, 3)]:
    full[deg] = complex_dims(deg, tag="full")
print("== C-odd / C-even subcomplexes, weight (3,3) ==")
for par in (-1, 1):
    complex_dims((3, 3), parity=par, tag="C=%+d" % par)
print("== SL2-invariant subcomplex ==")
inv33 = complex_dims((3, 3), inv=True, tag="SL2-inv")
for par in (-1, 1):
    complex_dims((3, 3), parity=par, inv=True, tag="SL2-inv C=%+d" % par)
print("== H^4 at the other bidegrees of total degree 6 (SL2 structure) ==")
for deg in [(4, 2), (5, 1), (6, 0)]:
    for k in (3, 4):
        pass
    cs3 = span_set(3, deg, SL)
    cs4 = span_set(4, deg, SL)
    B3, B4 = formal_basis(cs3, 3), formal_basis(cs4, 4)
    dB3 = [d(w, 3, SL) for w in B3]
    dB4 = [d(w, 4, SL) for w in B4]
    for N in NS:
        c4 = vla.rankQ(rows(B4, N, 4)) if B4 else 0
        r3 = vla.rankQ(rows(dB3, N, 4)) if B3 else 0
        r4 = vla.rankQ(rows(dB4, N, 5)) if B4 else 0
        print("  deg %s N=%d: dim C^4 %d, rank d3 %d, rank d4 %d, H^4 = %d" % (deg, N, c4, r3, r4, c4 - r3 - r4), flush=True)

print("== classes at (3,3) ==")
res, B, dB = full[(3, 3)]
for N in NS:
    Bd3 = rows(dB[3], N, 4)                    # B^4
    rB = vla.rankQ(Bd3)
    for name, c in [("theta'", thp), ("{D,iota}", DI), ("o_3' = {D,iota}+theta'", o3p), ("-{D,iota}+theta'", lin((1, thp), (-1, DI)))]:
        cr = rows([c], N, 4)[0]
        closed = not rows([d(c, 4, SL)], N, 5)[0]
        r2 = vla.rankQ(Bd3 + [cr])
        print("  N=%d %s: nonzero %s, closed %s, rank[B^4|c] = %d, rank B^4 = %d -> exact: %s" % (N, name, bool(cr), closed, r2, rB, r2 == rB))
    # generator h of H^4: a cocycle not in B^4
    Z = None
    import sympy as sp
    d4rows = rows(dB[4], N, 5)
    # kernel of the map C^4 -> C^5 on the formal basis B[4], coordinates at N
    cols = {}
    for r in d4rows:
        for kk in r:
            cols.setdefault(kk, len(cols))
    M = sp.zeros(len(cols), len(B[4]))
    for i, r in enumerate(d4rows):
        for kk, x in r.items():
            M[cols[kk], i] = sp.Rational(x.numerator, x.denominator)
    from sympy.polys.matrices import DomainMatrix
    from sympy import QQ
    ns = DomainMatrix.from_Matrix(M).convert_to(QQ).nullspace().to_Matrix()
    zs = []
    for j in range(ns.rows):
        vec = ns.row(j)
        z = {}
        for i in range(len(B[4])):
            if vec[i] != 0:
                z = lin((1, z), (Fr(int(vec[i].p), int(vec[i].q)), at_N(B[4][i], N)))
        zs.append(z)
    # exact closedness of the kernel vectors at N (coordinates)
    okc = all(not rows([d(z, 4, SL)], N, 5)[0] for z in zs)
    zrows = rows(zs, N, 4)
    rZ = vla.rankQ(zrows)
    gens = [z for z, zr in zip(zs, zrows) if vla.rankQ(Bd3 + [zr]) > rB]
    h = gens[0]
    hr = rows([h], N, 4)[0]
    Ch = rows([Cop(h, 4)], N, 4)[0]
    diff = {kk: hr.get(kk, 0) - Ch.get(kk, 0) for kk in set(hr) | set(Ch)}
    diff = {kk: v for kk, v in diff.items() if v}
    sumv = {kk: hr.get(kk, 0) + Ch.get(kk, 0) for kk in set(hr) | set(Ch)}
    sumv = {kk: v for kk, v in sumv.items() if v}
    print("  N=%d: dim Z^4 = %d (kernel vectors exactly closed: %s), rank B^4 = %d -> H^4 = %d;  generator h: h - Ch in B^4: %s, h + Ch in B^4: %s"
          % (N, rZ, okc, rB, rZ - rB, vla.rankQ(Bd3 + [diff]) == rB if diff else True,
             (vla.rankQ(Bd3 + [sumv]) == rB) if sumv else True), flush=True)
    # SL2-weight of h: E h and F h exact?
    Eh = rows([sl2_E(h, 4)], N, 4)[0]
    print("      h's momentum polynomials SL2-invariant (E h = 0 as cochain): %s" % (not Eh))
print("time %.1fs" % (time.time() - t0))
