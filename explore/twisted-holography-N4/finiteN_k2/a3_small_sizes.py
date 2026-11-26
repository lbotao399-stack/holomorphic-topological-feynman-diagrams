"""Chain-space dimensions of the HA^! (small) model, n+ = n- = 0 sectors.
usage: python3 a3_small_sizes.py k N 't1,t2,t3,0,0' gmin gmax"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_sizes as bs
k = int(sys.argv[1]); N = int(sys.argv[2])
lam = tuple(int(x) for x in sys.argv[3].split(','))
g0, g1 = int(sys.argv[4]), int(sys.argv[5])
ghs = list(range(g0, g1 + 1))
t0 = time.time()
tcl, wcl = bs.small_model_classes(lam)
if k == 2:
    extra = [((N+1, N-1), (N, N)), ((N, N), (N+1, N-1)), ((N+1, N-1), (N+1, N-1)), ((N+1, N), (N+1, N))]
else:
    extra = [((N + 1,), (N + 1,))]
d = bs.chain_dims_classes(lam, k, ghs, N, tcl, wcl, extra)
T = ((N,) * k, (N,) * k)
print(f"# small model k={k} N={N} sector {lam} ({time.time()-t0:.1f}s)")
for key in [T] + extra:
    print(f"  Lc={key[0]} Rc={key[1]}: " + "  ".join(f"g{g}:{int(d[(key[0], key[1], g)])}" for g in ghs))
if k == 2:
    aL, aR, aLR = extra[0], extra[1], extra[2]
    inv = [int(d[(T[0],T[1],g)] - d[(aL[0],aL[1],g)] - d[(aR[0],aR[1],g)] + d[(aLR[0],aLR[1],g)]) for g in ghs]
    print("  invariants per g: " + "  ".join(f"g{g}:{v}" for g, v in zip(ghs, inv)))
