"""S2: precise counting  dim g(beta-1^5) - dim coker(beta) = B^1(beta) - B^1(gamma(beta)),  canonical pairing radicals,
flat sources of the U(1) letters, orphans, and a weight-by-weight table.  Exact.

python3 run_S2.py 4 > S2.log
"""
import sys, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la, symmetry as sy

t0 = time.time()
R = int(sys.argv[1]) if len(sys.argv) > 1 else 3
def inR(b):
    return all(x in (0, 1) for x in b[2:]) and min(b[:2]) >= 1 and b != (1,) * 5
def typ(b):
    return "M1" if sy.in_M1(b) else ("1^5" if b == (1,) * 5 else ("R'" if inR(b) else "-"))
def Wn_of(b):
    return (b[0] - 1, b[1] - 1) + tuple(-x - 1 for x in b[2:])
def pf(om):
    out = []
    for (m, S, J, a), c in sorted(om.items()):
        f = []
        for i in range(2):
            if m[i]: f.append("z%d%s" % (i + 1, "^%d" % m[i] if m[i] > 1 else ""))
        f += ["dz%d" % (i + 1) for i in S] + ["t%d" % (I + 1) for I in J]
        for I in range(3):
            if a[I]: f.append("dt%d%s" % (I + 1, "^%d" % a[I] if a[I] > 1 else ""))
        out.append("%s%s" % (("+" if c > 0 else "-") + (str(abs(c)) + "*" if abs(c) != 1 else ""), "*".join(f) or "1"))
    return " ".join(out) or "0"

print("label beta = W_r + 1^5 (regular, C^{2|3} weight gamma = (beta_z - 1, 1 - beta_theta)); Laurent Wn = (beta_z - 1, -beta_theta - 1)")
print("P_D(alpha, zeta) = P_A(alpha, zeta_{p>=2}) (beta_z != 0), K(alpha, zeta_{p>=2}) (beta_z = 0): the flat Koszul pairing with B^1 projected out")
stats = dict(n=0, z=0, form=0, excess=0, rad=0)
cnt = {}
for b in itertools.product(range(R + 1), repeat=5):
    if list(b[:2]) != sorted(b[:2], reverse=True) or list(b[2:]) != sorted(b[2:], reverse=True):
        continue
    Wr = tuple(x - 1 for x in b)
    Wn = Wn_of(b)
    g = len(sy.g_basis(Wr))
    zA0 = len(sy.Z_A0(Wr))
    Z2, im = bx.coker_phi_t0(Wn)
    zA2 = len(Z2)
    cok = zA2 - la.rank(im)
    B1, B1g = sy.B1_dim(b), sy.B1_dim(sy.gamma_of_label(b))
    d = sy.projected_pairing_data(b)
    stats['n'] += 1
    stats['z'] += zA0 != zA2
    stats['form'] += (g != zA0 - B1g) or (cok != zA2 - B1)
    stats['excess'] += (g - cok) != B1 - B1g or (g - cok) != (1 if (sy.in_M1(b) or b == (1,) * 5) else 0)
    stats['rad'] += d['left_rad'] != B1 or d['right_rad'] != (1 if inR(b) else 0) or d['dim_coker'] != cok
    key = typ(b)
    cnt.setdefault(key, [0, 0, 0])
    cnt[key][0] += 1
    cnt[key][1] += d['left_rad']
    cnt[key][2] += d['right_rad']
print("[chk] beta in [0,%d]^5 mod S2xS3: %d weights" % (R, stats['n']))
print("      dim Z(A0)(beta-1^5) != dim Z(A2)(Wn): %d" % stats['z'])
print("      dim g != dim Z(A0) - B^1(gamma(beta))  or  dim coker != dim Z(A2) - B^1(beta): %d" % stats['form'])
print("      dim g - dim coker != B^1(beta) - B^1(gamma(beta)) = [M1] + [1^5]: %d" % stats['excess'])
print("      P_D radicals: left != B^1(beta) or right != [beta in R']: %d" % stats['rad'])
print("      by type: n, sum dim(left radical), sum dim(right radical):", cnt)

# flat sources of the U(1) letters
n = bad = 0
for bz in itertools.product(range(R + 1), repeat=2):
    for bt in itertools.product(range(2), repeat=3):
        b = bz + bt
        if b == (0,) * 5:
            continue
        us = sy.u1_source(b) if bz != (0, 0) else sy.ks_source_z0(b)
        pred = sy.u1_formula(b)
        n += 1
        bad += not (len(us) == 1 and sy.proportional(us[0], pred) and sy.alpha0_is_constant(pred) and not bx.delta(pred))
print("[chk] annihilator of B^{>=2}(beta) in Z(A0)(beta-1^5) (P_A; K at beta_z = 0) = C u_beta, u_beta = theta^{1-beta_theta} H_{z^beta_z}")
print("      (H_1 := d/dz1 ^ d/dz2), and u_beta in g: %d B^1-weights with beta_z in [0,%d]^2, %d failures" % (n, R, bad))

# orphans
print("orphans (right radical of P_D in coker) and the broken shift f_beta = z^{beta_z-1} theta^{1-beta_theta}:")
for b in [(1, 1, 0, 0, 0), (1, 1, 1, 0, 0), (2, 1, 1, 1, 1), (1, 2, 0, 1, 1), (2, 2, 1, 1, 0), (1, 1, 1, 1, 1)]:
    Wr = tuple(x - 1 for x in b)
    Wn = Wn_of(b)
    g = sy.g_basis(Wr)
    if b == (1,) * 5:
        E = la.Echelon(); E.add(sy.UNIT)
        g = [v for v in g if E.add(v)[0]]
    byp = sy.degree_split_Z2(Wn)
    high = [z for p, zs in byp.items() if p >= 2 for z in zs]
    cols = [{i: bx.pair(bx.fischer_res, a, z) for i, a in enumerate(g)} for z in high]
    ker = la.kernel(cols)
    rad = []
    for comb in ker:
        v = {}
        for j, c in comb.items():
            la.vadd(v, high[j], c)
        rad.append(la.clean(v))
    K = tuple(I for I in range(3) if b[2 + I] == 0)
    f = sy.bulk_of_pv((b[0] - 1, b[1] - 1), K, (0, 0, 0), ())
    pairs_f = [bx.pair(bx.fischer_res, f, r) for r in rad]
    print("   beta=%s %s: dim g%s=%d coker=%d; orphan(s): %s ; P_A(f_beta, orphan)=%s ; f_beta=%s" % (
        b, typ(b), "'" if b == (1,) * 5 else "", len(g), len(high), [pf(dr.Phi(r)) for r in rad], pairs_f, sy.pretty_pv(f)))

# weight-by-weight table
basis = sy.psl33_basis()
lab = {}
for (name, M, v) in basis:
    W = bx.weight(next(iter(v))) if name[0] == 'E' else sy.ZERO5
    lab.setdefault(sy.label_of_regular(W), []).append(name)
print("\nweight-by-weight table (beta up to S2xS3, |beta| <= 4, plus the psl(3|3) labels):")
print("   beta             type  U(N)tr  g  coker  B1  B1(gam) exc  psl(3|3) at label        u_beta (U(1)-letter source)         f_beta (Higgsed shift)")
rows = set()
for b in itertools.product(range(4), repeat=5):
    if sum(b) <= 4 and list(b[:2]) == sorted(b[:2], reverse=True) and list(b[2:]) == sorted(b[2:], reverse=True):
        rows.add(b)
rows |= set(lab)
for b in sorted(rows, key=lambda x: (sum(x), x)):
    Wr = tuple(x - 1 for x in b)
    Wn = Wn_of(b)
    g = len(sy.g_basis(Wr))
    Z2, im = bx.coker_phi_t0(Wn)
    cok = len(Z2) - la.rank(im)
    Bl = dr.B_list(b)
    un = sum(Bl[1:]) + (1 if b == (0,) * 5 else 0)
    B1, B1g = sy.B1_dim(b), sy.B1_dim(sy.gamma_of_label(b))
    u = sy.pretty_pv(sy.u1_formula(b)) if B1 else ""
    K = tuple(I for I in range(3) if b[2 + I] == 0)
    fb = sy.pretty_pv(sy.bulk_of_pv((b[0] - 1, b[1] - 1), K, (0, 0, 0), ())) if B1g else ("1 (unbroken)" if b == (1,) * 5 else "")
    print("   %-16s %-4s %6d %3d %5d %3d %5d %4d  %-24s %-35s %s" % (b, typ(b), un, g, cok, B1, B1g, g - cok,
          ",".join(lab.get(b, [])), u, fb))
print("time %.1fs" % (time.time() - t0))
