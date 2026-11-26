"""k = 2 meson complex (bg_complex_k) for the 1/16-BPS word algebra vs U(N) H(Q_0) in low sectors.
usage: python3 b2_bg_k2_small.py N Jmax"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_complex_k as bk, uN_q0 as uq
N = int(sys.argv[1]); J = int(sys.argv[2])
alg, keys = bw.n4_letters(J)
bg = bk.BGk(alg, 2)
exec(open("b2_bg_vs_UN.py").read().split("print(f\"# N")[0].split("bg = bc.BG1(alg)")[1])   # reps()
print(f"# N = {N}, k = 2: sector | U(N) {{Y: dim}} | BG k=2 {{Y: dim}} | agree | dim V_T | time")
for j in range(2, J + 1):
    for lam in reps(j):
        t0 = time.time()
        Th = sum(lam[:3])
        hu = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, N).items() if v[2]}
        res = bk.homology(bg, lam, N)
        hb = {g + Th: v[3] for g, v in res.items() if v[3]}
        print(f"{j:3d} | {lam} | {hu} | {hb} | {'Y' if hu == hb else 'NO'} | {sum(v[0] for v in res.values())} | {time.time()-t0:.1f}s", flush=True)
