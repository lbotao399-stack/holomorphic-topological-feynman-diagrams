"""k = 3: kernel vs Python bg_complex_k (tiny sectors), and chain_dims_general vs kernel enumeration."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_complex_k as bk, bg_fast as bf, bg_sizes as bs
alg, _ = bw.n4_letters(6)
for N in (0, 1):
    bg = bk.BGk(alg, 3)
    for lam in [(1,0,0,0,0), (0,0,0,1,0), (2,0,0,0,0), (1,1,0,0,0)]:
        t0 = time.time()
        ref = bk.homology(bg, lam, N)
        t1 = time.time()
        fast = bf.homology(alg, lam, 3, N, ghs=sorted(ref), mode=0)
        fast1 = bf.homology(alg, lam, 3, N, ghs=sorted(ref), mode=1)
        ok = all(tuple(fast[g]) == tuple(ref[g]) for g in ref) and all(fast1[g][3] == ref[g][3] for g in ref)
        gd = bs.chain_dims_general(alg, lam, 3, sorted(ref), N)
        okd = all(gd[((N,)*3, (N,)*3, g)] == ref[g][0] for g in ref)
        print(N, lam, 'homology', 'OK' if ok else ('MISMATCH', ref, fast, fast1), 'dims', 'OK' if okd else 'BAD', {g: ref[g] for g in ref}, '%.1fs' % (t1 - t0), flush=True)
