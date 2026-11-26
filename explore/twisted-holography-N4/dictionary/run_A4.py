"""A4: supercharacter of B^{>=2}(C^{2|3}) with the LLC charge map vs the KMMR single-particle index.

python3 run_A4.py [tmax] > A4.log
"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + '/code')
from fractions import Fraction
from holography import derham as dr, charges as ch
from holography.linalg import vadd

t0 = time.time()
tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 18

chi = {p: {} for p in range(0, 4 * tmax)}
for beta in ch.sector_list(tmax):
    npm, theta = beta[:2], beta[2:]
    for p, dim in enumerate(dr.B_list(beta)):
        if not dim:
            continue
        sign, mono = ch.kmmr_monomial(theta, npm, p)
        assert mono[0] <= tmax
        vadd(chi[p], {mono: Fraction(sign * dim)})

Bge1, Bge2 = {}, {}
for p, s in chi.items():
    if p >= 1:
        vadd(Bge1, s)
    if p >= 2:
        vadd(Bge2, s)
Isp = ch.kmmr_Isp(tmax)
f = ch.kmmr_f(tmax)
Isp_minus_f = vadd(dict(Isp), f, -1)


def first_diff(a, b):
    d = vadd(dict(a), b, -1)
    if not d:
        return None
    k = min(d, key=lambda m: (m[0], m))
    return k, d[k], len(d)


print(f"series through t^{tmax}; monomial (t,y,v,w) = t^(2Theta+3n) y^(n+-n-) v^(th1-th2) w^(th2-th3), sign (-1)^(Y+Theta)")
print("[chk] sch(B^1) - f (KMMR 4.2):", first_diff(chi[1], f))
print("[chk] sum_{p>=1} sch(B^p) - I_sp (KMMR 4.14):", first_diff(Bge1, Isp))
print("[chk] sum_{p>=2} sch(B^p) - (I_sp - f):", first_diff(Bge2, Isp_minus_f))
print("[chk] sum_{p>=2} sch(B^p) - I_sp: first difference", first_diff(Bge2, Isp))
nterms = sum(1 for _ in Bge2)
print(f"number of monomials in sum_(p>=2) sch(B^p) through t^{tmax}: {nterms}")
# low-order display of sum_{p>=2} sch B^p (t-degree <= 6)
low = sorted((k, v) for k, v in Bge2.items() if k[0] <= 6)
print("sum_{p>=2} sch(B^p), t-degree <= 6:")
for k, v in low:
    print("   ", v, k)
print(f"time {time.time()-t0:.1f}s")
