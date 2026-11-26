"""k = 1 (both primes) and U(2) E_1 (uN_q0, p) for a list of sectors.
usage: python3 c7_k1_list.py 's1' 's2' ..."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf, uN_q0 as uq, rank_mod
for s in sys.argv[1:]:
    lam = tuple(int(x) for x in s.split(','))
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:]); Th = sum(lam[:3])
    alg, _ = bw.n4_letters(L)
    t0 = time.time()
    tb = bf.sector_tables(alg, lam)
    hk = {}
    for p in (bf.P1, bf.P2):
        h = bf.homology(alg, lam, 1, 2, p=p, tables=tb)
        hk[p] = {g + Th: v[3] for g, v in h.items() if v[3]}
    t1 = time.time()
    hu = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, 2).items() if v[2]}
    print(f"{L:3d} | {lam} | k=1 {hk[bf.P1]} (p' same: {hk[bf.P1] == hk[bf.P2]}) | U(2) {hu} | k=1 = U(2): {hk[bf.P1] == hu} | {t1-t0:.1f}s + {time.time()-t1:.1f}s", flush=True)
