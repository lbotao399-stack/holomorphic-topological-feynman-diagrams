"""k = 2 meson complex for U(2) in one sector, ghost numbers gmin..gmax, at prime p (bg_fast kernel, block mode).
usage: python3 c5_k2_sector.py 't1,t2,t3,n+,n-' gmin gmax [p]"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf
lam = tuple(int(x) for x in sys.argv[1].split(','))
g0, g1 = int(sys.argv[2]), int(sys.argv[3])
p = int(sys.argv[4]) if len(sys.argv) > 4 else bf.P1
L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
Th = sum(lam[:3])
alg, _ = bw.n4_letters(L)
tag = "_".join(map(str, lam))
t0 = time.time()
h = bf.homology(alg, lam, 2, 2, p=p, ghs=list(range(g0, g1 + 1)), mode=1, log=f"c5_k2_{tag}_p{p}.raw")
print(f"# k=2 N=2 sector {lam} level {L} p={p} ({time.time()-t0:.1f}s)")
for g in sorted(h):
    print(f"  g={g} Y={g+Th}: dimV_T={h[g][0]} rankS={h[g][1]} dimR={h[g][0]-h[g][1]} rankd={h[g][2]} H={h[g][3]}", flush=True)
