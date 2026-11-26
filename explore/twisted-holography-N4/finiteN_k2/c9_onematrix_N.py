"""One-matrix sector (n,0,0;0,0) (= A = C[theta] of B-G App. C) for U(N), k = 2, all ghost numbers.
usage: python3 c9_onematrix_N.py N n gmin gmax [p]"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf
N = int(sys.argv[1]); n = int(sys.argv[2]); g0, g1 = int(sys.argv[3]), int(sys.argv[4])
p = int(sys.argv[5]) if len(sys.argv) > 5 else bf.P1
lam = (n, 0, 0, 0, 0)
alg, _ = bw.n4_letters(2 * n)
t0 = time.time()
h = bf.homology(alg, lam, 2, N, p=p, ghs=list(range(g0, g1 + 1)), mode=1, log=f"c9_N{N}_n{n}_p{p}.raw")
print(f"# k=2 N={N} x^{n} p={p} ({time.time()-t0:.1f}s)  H_g = { {g: v[3] for g, v in h.items() if v[3]} }")
for g in sorted(h):
    print(f"  g={g}: dimV_T={h[g][0]} dimR={h[g][0]-h[g][1]} rankd={h[g][2]} H={h[g][3]}", flush=True)
