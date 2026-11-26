"""B2/B3: Budzik-Gaiotto k = 1 meson complex for the 1/16-BPS word algebra vs the finite-N tree-level
cohomology H(Q_0) of U(N), sector by sector.

usage: python3 b2_bg_vs_UN.py N Jmin Jmax [field: Q | p (default 2^31-1) | p2 (2147483629)]
Sectors lam = (theta1,theta2,theta3;n+,n-) with level 2*Theta + 3(n+ + n-) = j, representatives
theta1 >= theta2 >= theta3, n+ >= n-.  Ghost number g of the meson complex <-> Y = g + Theta.
"""
import sys
import time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_complex as bc, uN_q0 as uq, rank_mod

N = int(sys.argv[1])
J0, J1 = int(sys.argv[2]), int(sys.argv[3])
FIELD = {'Q': 'Q', 'p': rank_mod.P1, 'p2': rank_mod.P2}[sys.argv[4] if len(sys.argv) > 4 else 'p']
alg, keys = bw.n4_letters(J1)
assert not bw.check_Q2(alg)
bg = bc.BG1(alg)


def reps(j):
    out = []
    for n in range(j // 3 + 1):
        r = j - 3 * n
        if r % 2:
            continue
        Th = r // 2
        for npl in range(n, (n + 1) // 2 - 1, -1):
            nmi = n - npl
            if npl < nmi:
                continue
            for t1 in range(Th, -1, -1):
                for t2 in range(min(t1, Th - t1), -1, -1):
                    t3 = Th - t1 - t2
                    if t3 > t2 or t3 < 0:
                        continue
                    out.append((t1, t2, t3, npl, nmi))
    return out


print(f"# N = {N}: k=1 Budzik-Gaiotto homology H_g vs U({N}) H^Y(Q_0), Y = g + Theta; field = {FIELD}")
print("# columns: level | sector | U(N) {Y: dim H^Y} | BG {Y: dim H_g} | agree | time")
tot_agree = tot = 0
for j in range(J0, J1 + 1):
    for lam in reps(j):
        t0 = time.time()
        Th = sum(lam[:3])
        hu = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, N, field=FIELD).items() if v[2]}
        hb_raw = bc.homology(bg, lam, N, field=FIELD)
        hb = {g + Th: v[3] for g, v in hb_raw.items() if v[3]}
        ok = (hu == hb)
        tot += 1
        tot_agree += ok
        sizes = sum(v[0] for v in hb_raw.values())
        print(f"{j:3d} | {lam} | {hu} | {hb} | {'Y' if ok else 'NO'} | dimC'={sizes} | {time.time()-t0:.1f}s",
              flush=True)
print(f"# agree {tot_agree}/{tot}")
