"""(f) Formal sl_N (kappa = tr), weight (3,3): C^3, C^4, B^4 = d C^3, Z^4 = ker d_4 (C^5 coordinates), H^4,
dual functional Phi (Phi(B^4) = 0, Phi(h) = 1), classes of
    theta'  (two-loop graphs with one D vertex),   {D, iota},   o_3' = {D, iota} + theta'
and of candidate basis classes.  O_3 cochain: o_3 = eps{j,i_2} + theta = a pi nu ({D,iota} + theta'), j = aD,
i_2 = eps pi nu iota.  Ranks certified (la.certified_rank: two primes + exact minor + exact kernel).
Saves data/h4w3_sl.pkl."""
import time
from common import *

t0 = time.time()
sl = True
deg = (3, 3)
th = load('theta.pkl')
i2 = load('i2_sl.pkl')
iota = i2['iota']
D = J_cochain('d', sl)
thp = ce.restrict_sl(th['th']['d'])

C3, k3 = spaces.basis(3, deg, sl)
C4, k4 = spaces.basis(4, deg, sl)
k5 = spaces.keys_of(5, deg, sl)
print("formal sl_N weight %s: dim C^2 = %d, dim C^3 = %d, dim C^4 = %d; ambient coordinates %d, %d, %d"
      % (deg, len(spaces.basis(2, deg, sl)[0]), len(C3), len(C4), len(k3), len(k4), len(k5)), flush=True)
B4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 3, sl)) for c in C3]
_, MB = la.matrix_from(B4, k4)
cB = la.certified_rank(MB)
print("rank d3 = %d (minor det %s, mod-p ranks %s)" % (cB['rank'], cB['det'], cB['ranks_modp']), flush=True)
dC4 = [spaces.const_part(ce.d_ce(spaces.to_np(c), 4, sl)) for c in C4]
_, MZ = la.matrix_from(dC4, k5)
cZ = la.certified_rank(MZ)
print("rank d4 = %d (minor det %s, mod-p ranks %s, %d exact kernel vectors)" % (
    cZ['rank'], cZ['det'], cZ['ranks_modp'], len(cZ['kernel'])), flush=True)
Z4 = []
for v in cZ['kernel']:
    z = {}
    for j, cj in enumerate(v):
        if cj:
            for kk, x in C4[j].items():
                z[kk] = z.get(kk, 0) + cj * x
    Z4.append({kk: x for kk, x in z.items() if x})
print("dim Z^4 = %d (closed exactly: %s), dim B^4 = %d, dim H^4 = %d" % (
    len(Z4), all(is_zero(ce.d_ce(spaces.to_np(z), 4, sl)) for z in Z4), cB['rank'], len(Z4) - cB['rank']), flush=True)


def npow_parts(c):
    pw = sorted({p for d in c.values() for p in d})
    return {p: {kk: d[p] for kk, d in c.items() if d.get(p)} for p in pw}


def exact_test(cconst):
    _, M = la.matrix_from(B4 + [cconst], k4)
    return la.certified_rank(M)['rank'] == cB['rank']


# candidate basis classes
G2 = ce.alt_raw(ce.raw_from(ce.poly_mul(ce.bracket_poly(0, 1, 4), ce.bracket_poly(2, 3, 4)),
                            {((0, 1), (2, 3)): {0: Fr(1)}}), 4, sl=True)
obs = spaces.obstruction_cochains()
F1 = {kk: {0: d[1]} for kk, d in ce.restrict_sl(obs[(2, 2, 3)]['F']).items() if d.get(1)}
cands = [("{D, G}, G = alt([p1p2][p3p4] tr(12)tr(34))", ce.bracket(D, 3, G2, 4, sl=True)),
         ("{D, F_1}, F_1 = N^1-part of F_223", ce.bracket(D, 3, F1, 4, sl=True)),
         ("{D, [p1p2]^2 f}", ce.bracket(D, 3, i2['u2f'], 3, sl=True))]
cl = {}
for nm, c in cands:
    parts = npow_parts(c)
    st = {p: (is_zero(ce.d_ce(spaces.to_np(x), 4, sl)), exact_test(x)) for p, x in parts.items()}
    print("candidate %s: N-power parts (closed, exact): %s" % (nm, st), flush=True)
    cl[nm] = (c, parts, st)

# basis class h: N^0 part of {D, [p1p2]^2 f} if not exact, else first non-exact candidate part
h = None
for nm, (c, parts, st) in cl.items():
    for p, (cl_, ex) in sorted(st.items()):
        if cl_ and not ex and h is None:
            h = (nm + " [N^%d part]" % p, parts[p])
if h is None:
    for i, z in enumerate(Z4):
        if not exact_test(z) and h is None:
            h = ("z%d (kernel vector of d4 in the basis C^4)" % i, z)
print("basis class h =", h[0], flush=True)
print("  h =", show(h[1], 4, maxterms=12))
hc = h[1]
# dual functional
Gm = B4 + [hc]
_, MG = la.matrix_from(Gm, k4)
cG = la.certified_rank(MG)
assert cG['rank'] == len(Gm)
rows = cG['pivrows']
sub = [[MG[i][j] for j in range(len(Gm))] for i in rows]
subT = [[sub[i][j] for i in range(len(rows))] for j in range(len(Gm))]
x = la.solve_square(subT, [Fr(0)] * len(B4) + [Fr(1)])
Phi = {k4[rows[i]]: x[i] for i in range(len(rows)) if x[i] != 0}
coord = lambda c: sum(Fr(c.get(kk, 0)) * v for kk, v in Phi.items())
print("Phi(B^4) = 0: %s ; Phi(h) = 1: %s" % (all(coord(b) == 0 for b in B4), coord(hc) == 1))
print("Phi on Z^4 basis:", [str(coord(z)) for z in Z4])


def cls(c, name):
    parts = npow_parts(c)
    out = {}
    for p, x in parts.items():
        assert is_zero(ce.d_ce(spaces.to_np(x), 4, sl)), (name, p)
        out[p] = coord(x)
    print("[%s] = %s [h]   (N-power: coefficient)" % (name, {p: str(v) for p, v in out.items() if v}), flush=True)
    return out


for nm, (c, parts, st) in cl.items():
    cls(c, nm)
r_th = cls(thp, "theta'")
br = ce.bracket(D, 3, iota, 3, sl=True)
r_br = cls(br, "{D, iota}")
o3p = ce.lin((1, br), (1, thp))
r_o3 = cls(o3p, "{D, iota} + theta'")
r_ctrl = cls(ce.lin((-1, br), (1, thp)), "control -{D, iota} + theta'")
save(dict(C3=C3, C4=C4, B4=B4, Z4=Z4, h=h, Phi=Phi, k4=k4, thp=thp, br=br, o3p=o3p,
          classes=dict(theta=r_th, bracket=r_br, o3=r_o3)), 'h4w3_sl.pkl')
print("time %.1fs" % (time.time() - t0))
