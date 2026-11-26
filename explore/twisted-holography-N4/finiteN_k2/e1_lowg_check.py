"""Claude check: k=2, N=2 one-matrix sectors (n,0,0;0,0) over all g = -2n..0 (n=5,6,7), and (8,0,0;0,0) at g = -16..-8
(the c5 run covered g = -8..0 only).  usage: python3 e1_lowg_check.py"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf, bg_sizes as bs
def run(lam, k, N, g0, g1, p=bf.P1):
    L = 2*sum(lam[:3]) + 3*sum(lam[3:]); Th = sum(lam[:3])
    alg, _ = bw.n4_letters(L)
    t0 = time.time()
    h = bf.homology(alg, lam, k, N, p=p, ghs=list(range(g0, g1 + 1)), mode=1, log=f"e1_low_{'_'.join(map(str,lam))}_k{k}_N{N}_p{p}.raw")
    print(lam, "k", k, "N", N, f"g={g0}..{g1}", f"p={p}", f"{time.time()-t0:.1f}s", flush=True)
    for g in sorted(h):
        print(f"  g={g} Y={g+Th}: dimV_T={h[g][0]} dimR={h[g][0]-h[g][1]} rankd={h[g][2]} H={h[g][3]}", flush=True)
for n in (5, 6, 7):
    run((n,0,0,0,0), 2, 2, -2*n, 0)
run((8,0,0,0,0), 2, 2, -16, -8)
