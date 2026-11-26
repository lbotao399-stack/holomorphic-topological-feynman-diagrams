"""A2: regular modes (ker phi_*) vs Laurent modes (coker phi_*): dimensions and pairings.  Exact.

python3 run_A2.py > A2.log
"""
import sys, itertools, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from holography import bulk_x as bx, derham as dr

t0 = time.time()
R = int(sys.argv[1]) if len(sys.argv) > 1 else 3


def in_M1(b):
    return b != (0,) * 5 and all(x in (0, 1) for x in b[2:]) and min(b[:2]) == 0


def in_R1(b):  # B^1 weights with both z-entries >= 1, excluding (1,1,1,1,1)
    return all(x in (0, 1) for x in b[2:]) and min(b[:2]) >= 1 and b != (1,) * 5


print("rule W_r = (W_n,z, -2-W_n,w) = beta(W_n) - (1,1,1,1,1); t^0;  dimK = dim ker phi_*(W_r), dimC = dim coker phi_*(W_n)")
print("pairing P = Taylor in (z, xi_z) [z^m <-> z^m, weight m!, xi_z equal] (x) Res_w [w^b <-> w^(-b-1), xi_w complementary]")
print("columns: beta, dimK, dimC, dimK-dimC, pred, dimB>=1, direct(descends, rank on K x coker), KS(indep, descends, rank),")
print("         closed(dim Z(A0)_Wr, dim Z(A2)_Wn, rank P), dim B^1, dim{k in K: P(k,B^1)=0}")
stats = dict(n=0, dim_mis=0, pred_mis=0, unm_mis=0)
types = {}
for b in itertools.product(range(R + 1), repeat=5):
    if list(b[:2]) != sorted(b[:2], reverse=True) or list(b[2:]) != sorted(b[2:], reverse=True):
        continue
    Wn = (b[0] - 1, b[1] - 1, -b[2] - 1, -b[3] - 1, -b[4] - 1)
    d = bx.pairing_data(Wn)
    B = dr.B_list(b)
    Bge1 = sum(B[1:]) + (1 if b == (0,) * 5 else 0)
    pred = (1 if in_M1(b) else 0) + (1 if b == (1,) * 5 else 0)
    stats['n'] += 1
    stats['dim_mis'] += d['dimK'] != d['dimC']
    stats['pred_mis'] += (d['dimK'] - d['dimC']) != pred
    stats['unm_mis'] += (d['dimK'] - Bge1) != (-1 if in_R1(b) else 0)
    zt = ('z0' if b[0] == 0 else 'z+') + ('z0' if b[1] == 0 else 'z+')
    tt = 'th01' if all(x in (0, 1) for x in b[2:]) else 'th2'
    key = (zt, tt, b == (0,) * 5, b == (1,) * 5)
    typ = types.setdefault(key, dict(n=0, eq=0, dir_ok=0, ks_ok=0, closed_perfect=0, K_perp_B1=0))
    typ['n'] += 1
    typ['eq'] += d['dimK'] == d['dimC']
    typ['dir_ok'] += d['dir_desc'] and d['dir_rank'] == d['dimK'] == d['dimC']
    typ['ks_ok'] += d['ks_desc'] and d['ks_rank'] == d['dimK'] == d['dimC']
    typ['closed_perfect'] += d['closed_rank'] == d['dimZ0'] == d['dimZ2']
    typ['K_perp_B1'] += d['orth_dim'] == d['dimK']
    print("  %s  %2d %2d %+d %d  %2d  dir(%s,%d)  KS(%s,%s,%d)  closed(%d,%d,%d)  B1=%d  Kperp=%d" % (
        b, d['dimK'], d['dimC'], d['dimK'] - d['dimC'], pred, Bge1, d['dir_desc'], d['dir_rank'],
        d['ks_indep'], d['ks_desc'], d['ks_rank'], d['dimZ0'], d['dimZ2'], d['closed_rank'], d['dimB1'], d['orth_dim']))
print(f"[chk] beta in [0,{R}]^5 up to S2xS3: {stats['n']} weights; dimK != dimC at {stats['dim_mis']}; "
      f"dimK-dimC != [beta in M1]+[beta=1^5] at {stats['pred_mis']};  "
      f"dimK-dim(B^>=1 + [beta=0]) != -[beta in R'] at {stats['unm_mis']}")
print("M1 = {beta != 0: beta_theta in {0,1}^3, min(beta_z) = 0};  R' = {beta_theta in {0,1}^3, beta_z >= 1, beta != 1^5}")
print("by type (z-pattern, theta-pattern, beta=0, beta=1^5): n, #(dimK=dimC), #(direct descends & nondegenerate),")
print("   #(KS descends & nondegenerate), #(P perfect on Z(A0) x Z(A2)), #(K orthogonal to B^1 = im phi_*):")
for k in sorted(types):
    print("  ", k, types[k])
print(f"time {time.time()-t0:.1f}s")
