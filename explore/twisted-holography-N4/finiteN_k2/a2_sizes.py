"""Chain-space dimensions dim V_{T,g} (and V_{T - ch(beta), g} for the b-images) of the k-meson complex.
usage: python3 a2_sizes.py k N 't1,t2,t3,n+,n-' gmin gmax"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_sizes as bs
k = int(sys.argv[1]); N = int(sys.argv[2])
lam = tuple(int(x) for x in sys.argv[3].split(','))
g0, g1 = int(sys.argv[4]), int(sys.argv[5])
L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
alg, _ = bw.n4_letters(L)
t0 = time.time()
ghs = list(range(g0, g1 + 1))
if k == 2:
    # T, T+alpha_L, T+alpha_R, T+alpha_L+alpha_R, and sources of b-images: T - ch(P0_xy) = T -+ alpha_L, etc.
    extra = [((N+1, N-1), (N, N)), ((N, N), (N+1, N-1)), ((N+1, N-1), (N+1, N-1)),
             ((N+1, N), (N+1, N))]
else:
    extra = [((N + 1,), (N + 1,))]
d = bs.chain_dims(alg, lam, k, ghs, N, extra_charges=extra)
T = ((N,) * k, (N,) * k)
print(f"# k={k} N={N} sector {lam} level {L}  ({time.time()-t0:.1f}s)")
for key in [T] + extra:
    row = [int(d[(key[0], key[1], g)]) for g in ghs]
    print(f"  Lc={key[0]} Rc={key[1]}: " + "  ".join(f"g{g}:{v}" for g, v in zip(ghs, row)))
if k == 2:
    # dim of the det^N-isotypic part (relative invariants) per g: Weyl formula for gl2 x gl2
    aL = ((N+1, N-1), (N, N)); aR = ((N, N), (N+1, N-1)); aLR = ((N+1, N-1), (N+1, N-1))
    inv = [int(d[(T[0],T[1],g)] - d[(aL[0],aL[1],g)] - d[(aR[0],aR[1],g)] + d[(aLR[0],aLR[1],g)]) for g in ghs]
    print("  invariants (det^N isotypic) per g: " + "  ".join(f"g{g}:{v}" for g, v in zip(ghs, inv)))
