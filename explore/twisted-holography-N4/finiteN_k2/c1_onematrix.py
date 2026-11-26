"""k-convergence test in the one-matrix sectors (n,0,0;0,0) of U(2): k = 1, 2 meson complexes.
usage: python3 c1_onematrix.py k n gmin gmax [p]"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf
k = int(sys.argv[1]); n = int(sys.argv[2]); g0, g1 = int(sys.argv[3]), int(sys.argv[4])
p = int(sys.argv[5]) if len(sys.argv) > 5 else bf.P1
lam = (n, 0, 0, 0, 0)
alg, _ = bw.n4_letters(2 * n)
t0 = time.time()
h = bf.homology(alg, lam, k, 2, p=p, ghs=list(range(g0, g1 + 1)), mode=1, log=f"c1_onematrix_k{k}_n{n}_p{p}.raw")
print(f"# k={k} N=2 sector {lam} p={p}  ({time.time()-t0:.1f}s)")
for g in sorted(h):
    print(f"  g={g} Y={g+n}: dimV_T={h[g][0]} rankS={h[g][1]} rankd={h[g][2]} H={h[g][3]}")
