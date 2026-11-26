"""For every sector of the k=1 scan logs: U(2) E_1 (uN_q0, p) and comparison with H^BG_{k=1}.
usage: python3 c8_compare_k1_UN.py log1 [log2 ...]"""
import sys, ast, time
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, uN_q0 as uq
done = set()
for fn in sys.argv[1:]:
    for line in open(fn):
        if line.startswith("#") or "|" not in line:
            continue
        p = [x.strip() for x in line.split("|")]
        lam = ast.literal_eval(p[1]); hb = ast.literal_eval(p[2])
        if lam in done:
            continue
        done.add(lam)
        L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
        alg, _ = bw.n4_letters(L)
        t0 = time.time()
        hu = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, 2).items() if v[2]}
        diff = {Y: hb.get(Y, 0) - hu.get(Y, 0) for Y in set(hb) | set(hu) if hb.get(Y, 0) != hu.get(Y, 0)}
        print(f"{L:3d} | {lam} | k=1 {hb} | U(2) {hu} | equal {hb == hu} | k=1 - U(2) = {diff} | {time.time()-t0:.1f}s", flush=True)
