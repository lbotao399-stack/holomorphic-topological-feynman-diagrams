"""Validation of bg_fast (C++ kernel) against every sector of the previous report:
k = 1 for U(1), U(2), U(3), U(4) (logs bg_vs_U*_k1*.log) and k = 2 for U(2) (bg_vs_U2_k2.log),
at both primes p = 2^31-1 and p' = 2147483629, mode 1 (reduced b-operators); k = 2 also mode 0."""
import sys, re, ast, time
from multiprocessing import Pool
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf
D = (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/explore/twisted-holography-N4/finiteN/"
jobs = []
for fn, N, k in [("bg_vs_U1_k1.log", 1, 1), ("bg_vs_U2_k1.log", 2, 1), ("bg_vs_U2_k1_L11_12.log", 2, 1),
                 ("bg_vs_U3_k1.log", 3, 1), ("bg_vs_U4_k1.log", 4, 1), ("bg_vs_U2_k2.log", 2, 2)]:
    for line in open(D + fn):
        if line.startswith("#") or "|" not in line:
            continue
        parts = [x.strip() for x in line.split("|")]
        lev = int(parts[0]); lam = ast.literal_eval(parts[1]); hu = ast.literal_eval(parts[2]); hb = ast.literal_eval(parts[3])
        jobs.append((fn, N, k, lev, lam, hu, hb))

def work(job):
    fn, N, k, lev, lam, hu, hb = job
    alg, _ = bw.n4_letters(lev)
    t0 = time.time()
    Th = sum(lam[:3])
    res = {}
    tb = bf.sector_tables(alg, lam)
    for p in (bf.P1, bf.P2):
        modes = (1, 0) if k == 2 else (1,)
        for mode in modes:
            h = bf.homology(alg, lam, k, N, p=p, mode=mode, tables=tb)
            res[(p, mode)] = {g + Th: v[3] for g, v in h.items() if v[3]}
    ok = all(r == hb for r in res.values())
    return fn, N, k, lev, lam, hb, hu, ok, res, time.time() - t0

if __name__ == "__main__":
    bf.build()
    nok = 0
    with Pool(12) as pool:
        for fn, N, k, lev, lam, hb, hu, ok, res, dt in pool.imap(work, jobs):
            nok += ok
            print(f"{fn} N={N} k={k} {lev:3d} {lam} prev={hb} U(N)={hu} fast={'same' if ok else res} {'OK' if ok else 'MISMATCH'} {dt:.1f}s", flush=True)
    print(f"# {nok}/{len(jobs)} sectors reproduced at p and p' (k=2 also with all b-operators)")
