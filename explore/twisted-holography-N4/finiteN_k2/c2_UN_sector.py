"""E_1^{U(2)} = H(Q_0) of the relative twisted U(2) complex per sector (uN_q0, evaluation on 2x2 matrices),
both primes.  usage: python3 c2_UN_sector.py N 't1,t2,t3,n+,n-' [...]"""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, uN_q0 as uq, rank_mod
N = int(sys.argv[1])
for s in sys.argv[2:]:
    lam = tuple(int(x) for x in s.split(','))
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
    alg, _ = bw.n4_letters(L)
    t0 = time.time()
    out = {}
    for p in (rank_mod.P1, rank_mod.P2):
        r = uq.cohomology_sector(alg, lam, N, field=p)
        out[p] = {Y: v[2] for Y, v in r.items() if v[2]}
        dims = {Y: v[0] for Y, v in r.items()}
    same = out[rank_mod.P1] == out[rank_mod.P2]
    print(f"{L:3d} | {lam} | U({N}) H^Y = {out[rank_mod.P1]} | p,p' agree: {same} | dim Inv^Y = {dims} | {time.time()-t0:.1f}s", flush=True)
