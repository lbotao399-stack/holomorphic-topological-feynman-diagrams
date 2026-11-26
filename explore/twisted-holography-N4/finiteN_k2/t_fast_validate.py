import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_complex_k as bk, bg_fast as bf
alg, keys = bw.n4_letters(8)
bg = bk.BGk(alg, 2)
for N in (1, 2):
    for lam in [(1,0,0,0,0),(2,0,0,0,0),(1,1,0,0,0),(1,0,0,1,0),(3,0,0,0,0),(2,1,0,0,0),(0,0,0,2,0),(0,0,0,1,1)]:
        t0 = time.time()
        ref = bk.homology(bg, lam, N)
        t1 = time.time()
        for mode in (0, 1):
            fast = bf.homology(alg, lam, 2, N, ghs=sorted(ref), mode=mode)
            ok = all(tuple(fast[g]) == tuple(ref[g]) if mode == 0 else fast[g][3] == ref[g][3] for g in ref)
            print(N, lam, 'mode', mode, 'OK' if ok else 'MISMATCH', {g: ref[g] for g in ref} if not ok else '', {g: fast[g] for g in fast} if not ok else '', '%.2fs/%.2fs' % (t1-t0, time.time()-t1), flush=True)
