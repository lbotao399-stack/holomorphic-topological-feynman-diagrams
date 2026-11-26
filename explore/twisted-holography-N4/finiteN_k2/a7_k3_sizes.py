"""k = 3 chain dims (exact sparse DP, bg_sizes.chain_dims_general) vs k = 1, 2 in BMN sectors (n,n,n;0,0),
for the ratio dim V^{k=3}/dim V^{k=2} at Y = Theta - 5 .. Theta - 3 ... (all g printed)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_sizes as bs
for lam in [(1,1,1,0,0), (2,1,1,0,0), (2,2,1,0,0), (2,2,2,0,0), (3,2,2,0,0)]:
    L = 2 * sum(lam[:3]); Th = sum(lam[:3])
    alg, _ = bw.n4_letters(L)
    ghs = list(range(-2 * Th, 1))
    t0 = time.time()
    out = {}
    for k in (1, 2, 3):
        if k <= 2:
            d = bs.chain_dims(alg, lam, k, ghs, 2)
            out[k] = {g: int(round(d[((2,)*k, (2,)*k, g)])) for g in ghs}
        else:
            d = bs.chain_dims_general(alg, lam, 3, ghs, 2)
            out[k] = {g: d[((2,)*3, (2,)*3, g)] for g in ghs}
    for g in ghs:
        if out[3][g]:
            print(f"{lam} g={g} Y={g+Th}: k1 {out[1][g]}  k2 {out[2][g]}  k3 {out[3][g]}  k3/k2 {out[3][g]/max(out[2][g],1):.1f}", flush=True)
    print(f"# {lam} {time.time()-t0:.1f}s", flush=True)
