"""A3: cyclic homology of A = C[z1,z2] (x) Lambda[theta1,theta2,theta3] from the Connes complex, the
cyclic-word CE model of planar single-trace cochains, and LLC sector labels.  Exact.

python3 run_A3.py > A3.log
"""
import sys, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from holography import cyclic as cy, derham as dr, charges as ch

t0 = time.time()


def canon_betas(total_max, theta_max=None):
    out = []
    for b in itertools.product(range(total_max + 1), repeat=5):
        if sum(b) > total_max:
            continue
        if list(b[:2]) != sorted(b[:2], reverse=True) or list(b[2:]) != sorted(b[2:], reverse=True):
            continue
        out.append(b)
    return out


# ------------------------------------------------ Connes complex: absolute and reduced HC_n, n <= 4
print("HC_n(A)_beta (Connes complex) vs HKR: absolute = B^{n+1}(beta) + [n even][beta=0]; reduced = B^{n+1}(beta)")
mis_a = mis_r = cnt_a = cnt_r = 0
rows = []
for beta in canon_betas(5):
    B = dr.B_list(beta)
    Bn = lambda n: B[n + 1] if n + 1 < len(B) else 0
    rec = []
    for n in range(5):
        hr, dr_ = cy.HC(beta, n, reduced=True)
        cnt_r += 1
        mis_r += hr != Bn(n)
        if sum(beta) <= 5:
            ha, da = cy.HC(beta, n)
            pa = Bn(n) + (1 if (beta == (0,) * 5 and n % 2 == 0) else 0)
            cnt_a += 1
            mis_a += ha != pa
        else:
            ha = None
        rec.append((ha, hr, Bn(n)))
    rows.append((beta, rec))
for beta, rec in rows:
    print("  beta=%s  [(HC_n, HCbar_n, B^{n+1}) n=0..4] = %s" % (beta, rec))
print(f"[chk] absolute: {cnt_a} (beta,n) slots, |beta|<=5, mismatches {mis_a};  reduced: {cnt_r} slots, mismatches {mis_r}")
# HC_n = Omega^n/dOmega^{n-1} (+) C[n even, n>=2, beta=0];  the form  Omega^n/dOmega^{n-1} (+) C[n even]  double counts at n=0
mis_om = mis_lit = cnt_om = 0
for beta, rec in rows:
    for n in range(5):
        quot = len(dr.forms(beta, n)) - dr.B_dim(beta, n)
        ha = rec[n][0]
        cnt_om += 1
        mis_om += ha != quot + (1 if (beta == (0,) * 5 and n % 2 == 0 and n >= 2) else 0)
        mis_lit += ha != quot + (1 if (beta == (0,) * 5 and n % 2 == 0) else 0)
print(f"[chk] HC_n vs dim Omega^n/dOmega^(n-1) + [n even, n>=2, beta=0]: {cnt_om} slots, mismatches {mis_om};"
      f"  vs Omega^n/dOmega^(n-1) + [n even, beta=0]: mismatches {mis_lit} (n=0, beta=0)")

# reduced at |beta| = 6 for n <= 4
mis6 = cnt6 = 0
for beta in canon_betas(6):
    if sum(beta) != 6:
        continue
    B = dr.B_list(beta)
    for n in range(5):
        hr, _ = cy.HC(beta, n, reduced=True)
        cnt6 += 1
        mis6 += hr != (B[n + 1] if n + 1 < len(B) else 0)
print(f"[chk] reduced, |beta|=6, n<=4: {cnt6} slots, mismatches {mis6}")

# ------------------------------------------------ cyclic-word CE model, LLC sectors
print()
print("planar single-trace CE cohomology (cyclic words, Q C = C.C), sector (theta1,theta2,theta3; n+,n-)_Y,")
print("charges (E; R1,R2,R3; J_L, J_R); columns: H_abs, H_rel, B^Y(beta), HC_{Y-1}, HCbar_{Y-1}")
mis = cnt = 0
for beta in canon_betas(5):
    B = dr.B_list(beta)
    for Y in range(1, 6):
        hrel, nw = cy.CE_single_trace(beta, Y, relative=True)
        habs, _ = cy.CE_single_trace(beta, Y, relative=False)
        BY = B[Y] if Y < len(B) else 0
        hc, _ = cy.HC(beta, Y - 1)
        hcr, _ = cy.HC(beta, Y - 1, reduced=True)
        cnt += 1
        bad = (hrel != BY) or (hcr != BY) or (habs != BY + (1 if beta == (0,) * 5 and Y % 2 == 1 else 0)) or (hc != habs)
        mis += bad
        if habs or hrel or BY:
            theta, npm = beta[2:], beta[:2]
            E, R, JL, JR = ch.llc_charges(theta, npm, Y)
            lab = "(%d,%d,%d;%d,%d)_%d" % (theta + npm + (Y,))
            chg = "(%s; %s,%s,%s; %s,%s)" % (E, R[0], R[1], R[2], JL, JR)
            print(f"  {lab:18s} {chg:34s} H_abs={habs} H_rel={hrel} B^Y={BY} HC={hc} HCbar={hcr}"
                  + ("  MISMATCH" if bad else ""))
print(f"[chk] CE single-trace sectors (|beta|<=5, Y<=5): {cnt}, mismatches {mis}")
print(f"time {time.time()-t0:.1f}s")
