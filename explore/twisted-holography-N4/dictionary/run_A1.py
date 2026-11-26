"""A1: full bulk complex T = A0[[t]] (+) A2[[t]], D = t del + N{F,-}; exact.

python3 run_A1.py  > A1.log
"""
import sys, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from fractions import Fraction
from holography import bulk_x as bx, derham as dr, linalg as la

t0 = time.time()
M1 = (-1, -1, -1, -1, -1)
M2 = (-2, -2, -2, -2, -2)


def canon(W):
    return tuple(sorted(W[:2], reverse=True)) + tuple(sorted(W[2:], reverse=True))


# ---------------------------------------------------------------- D^2 = 0, N-orders separately
def D_t(v):
    return bx.D(v, N=0)


def D_F(v):
    out = la.vadd(bx.D(v, N=1), D_t(v), -1)
    return out


n = bad = [0, 0, 0]
cnt = 0
bad = [0, 0, 0]
for Wh in itertools.product(range(-2, 3), range(-2, 3), range(-4, 2), range(-4, 2), range(-4, 2)):
    for h in range(-1, 4):
        for b in bx.T_basis(Wh, h):
            v = {b: Fraction(1)}
            cnt += 1
            bad[0] += bool(D_t(D_t(v)))
            bad[1] += bool(la.vadd(D_t(D_F(v)), D_F(D_t(v))))
            bad[2] += bool(D_F(D_F(v)))
print(f"[chk] D^2=0 on {cnt} basis elements of T (What in [-2,2]^2x[-4,1]^3, h in [-1,3]): "
      f"N^0 fails {bad[0]}, N^1 fails {bad[1]}, N^2 fails {bad[2]}")

# ---------------------------------------------------------------- BV bracket = explicit Schouten, all alpha in A0
cnt = badb = 0
for W in itertools.product(range(-1, 2), repeat=5):
    for b in bx.basis(0, W):
        a = {b: Fraction(1)}
        cnt += 1
        badb += bx.bv_bracket(bx.F, a) != bx.schouten(bx.F, a)
print(f"[chk] {{F,alpha}} (BV formula) = explicit Schouten on {cnt} A0 monomials (W in [-1,1]^5): fails {badb}")

# ---------------------------------------------------------------- Phi: (A2, del) -> (Omega(C^{2|3}), d) chain iso
cnt = badp = 0
for W in itertools.product(range(-1, 3), range(-1, 3), range(-4, 0), range(-4, 0), range(-4, 0)):
    for k in bx.basis(2, W):
        cnt += 1
        badp += dr.Phi(bx.delta({k: 1})) != dr.d_form(dr.Phi({k: 1}))
print(f"[chk] Phi o del = d o Phi on {cnt} A2 monomials: fails {badp};  Phi(F) = {dr.Phi(bx.F)}")

# ---------------------------------------------------------------- H(D) per (What, h)
print()
print("H^h(T_What) = coker-part (image of H(A2[[t]],t del)) + ker-part (image in H(A0[[t]],t del))")
print("prediction: coker-part = dim B^{>=2}(beta)+[beta=0] (h=0), [What=-1] (h>=1);"
      " ker-part = kerphi_formula(What+1) (h=-1), [What=-2] (h>=0)")
stats = dict(classes=0, nonzero=0, mis_coker=0, mis_ker=0, les=0, sum=0)
rows = []
seen = set()
Bcache = {}
for Wh in itertools.product(range(-2, 3), range(-2, 3), range(-5, 2), range(-5, 2), range(-5, 2)):
    c = canon(Wh)
    if c in seen:
        continue
    seen.add(c)
    beta = bx.beta_of(c)
    for h in range(-1, 4):
        r = bx.class_cohomology(c, h)
        stats['classes'] += 1
        if r['dimH'] == 0 and r['coker'] == 0 and r['ker'] == 0:
            pass
        else:
            stats['nonzero'] += 1
        if h == 0:
            if min(beta) < 0:
                pc = 0
            else:
                if beta not in Bcache:
                    Bcache[beta] = dr.B_list(beta)
                B = Bcache[beta]
                pc = sum(B[2:]) + (1 if beta == (0, 0, 0, 0, 0) else 0)
        elif h >= 1:
            pc = 1 if c == M1 else 0
        else:
            pc = 0
        pk = bx.ker_phi_formula(tuple(x + 1 for x in c)) if h == -1 else (1 if c == M2 else 0)
        stats['mis_coker'] += r['coker'] != pc
        stats['mis_ker'] += r['ker'] != pk
        stats['les'] += (r['les_coker'] != r['coker']) or (r['les_ker'] != r['ker'])
        stats['sum'] += r['dimH'] != r['coker'] + r['ker']
        if (r['dimH'] or pc or pk) and max(abs(x) for x in beta) <= 2 and min(beta) >= -1:
            rows.append((c, h, beta, r['dimH'], r['coker'], pc, r['ker'], pk))
print(f"[chk] classes (What up to S2xS3, What_z in [-2,2], What_w in [-5,1], h in [-1,3]): {stats['classes']}, "
      f"nonzero {stats['nonzero']}; coker-part mismatches {stats['mis_coker']}, ker-part mismatches {stats['mis_ker']}, "
      f"LES (phi_*) inconsistencies {stats['les']}, dimH != coker+ker {stats['sum']}")
print()
print("slot table (What, h, beta(What), dimH, coker, pred, ker, pred) for beta in [-1,2]^5, nonzero rows:")
for row in rows:
    print("  ", row)

# ---------------------------------------------------------------- ker phi_* by weight
print()
print("ker phi_*|t^0 by regular weight W_r (up to S2xS3, entries in [-1,2]); formula 2^(s-1) - [W_r_z>=0, W_r_w in {0,-1}^3, W_r!=0]:")
seen = set()
tab = []
for Wr in itertools.product(range(-1, 3), repeat=5):
    c = canon(Wr)
    if c in seen:
        continue
    seen.add(c)
    d = len(bx.ker_phi_t0(c))
    _, Z0 = bx.ker_del(0, c)
    tab.append((c, d, len(Z0), bx.ker_phi_formula(c)))
bad = sum(1 for x in tab if x[1] != x[3])
for x in tab:
    print("   W_r=%s  dim ker phi_*=%d  (dim ker del|A0=%d, formula %d)" % x)
print(f"[chk] ker phi_*|t^0 formula mismatches: {bad} of {len(tab)}")
print("t^k tower (k>=1): ker phi_* = C xi_top t^k at W_r=-1;  coker phi_* = C F t^k at W_n=-1  (h>=1 rows above)")
print(f"time {time.time()-t0:.1f}s")
