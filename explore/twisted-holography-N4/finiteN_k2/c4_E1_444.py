"""E_1 = H(Q_0) for U(2) and SU(2) in the O_0 sector (4,4,4;0,0) (uN_q0, both primes), plus the stable
(U(inf)) cohomology and the U(1) factor check E_1^{U(2)} = E_1^{SU(2)} (x) Sym(V_1)."""
import sys, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, uN_q0 as uq, rank_mod
lam = (4, 4, 4, 0, 0)
alg, _ = bw.n4_letters(24)
for special, name in ((True, "SU(2)"), (False, "U(2)")):
    for p in (rank_mod.P1, rank_mod.P2):
        t0 = time.time()
        r = uq.cohomology_sector(alg, lam, 2, special=special, field=p)
        print(f"{name} p={p}: H^Y = { {Y: v[2] for Y, v in r.items() if v[2]} }  dim Inv^Y = { {Y: v[0] for Y, v in r.items()} }  ({time.time()-t0:.1f}s)", flush=True)
