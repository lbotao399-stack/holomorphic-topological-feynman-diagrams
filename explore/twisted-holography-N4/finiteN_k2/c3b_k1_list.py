"""k = 1 meson complex for U(2), both primes, list of sectors (resume of c3_k1_scan).
usage: python3 c3b_k1_list.py nproc 's1' 's2' ..."""
import sys, time
from multiprocessing import Pool
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf
def work(s):
    lam = tuple(int(x) for x in s.split(','))
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:]); Th = sum(lam[:3])
    alg, _ = bw.n4_letters(L)
    t0 = time.time()
    tb = bf.sector_tables(alg, lam)
    out = {}
    for p in (bf.P1, bf.P2):
        h = bf.homology(alg, lam, 1, 2, p=p, tables=tb)
        out[p] = {g + Th: v[3] for g, v in h.items() if v[3]}
    chi = sum((-1) ** ((Y - Th) % 2) * d for Y, d in out[bf.P1].items())
    return f"{L:3d} | {lam} | {out[bf.P1]} | {out[bf.P1] == out[bf.P2]} | {chi} | {time.time()-t0:.1f}s"
if __name__ == "__main__":
    with Pool(int(sys.argv[1])) as pool:
        for line in pool.imap_unordered(work, sys.argv[2:]):
            print(line, flush=True)
