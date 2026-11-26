"""k = 1 meson complex (bg_fast kernel), U(2), all sector representatives at the given levels, both primes.
Records H_g^BG (Y = g + Theta), the Euler characteristic and the index target [Z_inf (1 + G_2^(1))]_lam.
usage: python3 c3_k1_scan.py Jmin Jmax maxsize"""
import sys, ast, time
from multiprocessing import Pool
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_fast as bf, bg_sizes as bs
res = {}
for line in open("a1_residual.log"):
    if line.startswith("#"): continue
    p = [x.strip() for x in line.split("|")]
    res[ast.literal_eval(p[1])] = (int(p[2]), int(p[3]), int(p[4]))

def work(lam):
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
    alg, _ = bw.n4_letters(L)
    Th = sum(lam[:3])
    t0 = time.time()
    tb = bf.sector_tables(alg, lam)
    out = {}
    for p in (bf.P1, bf.P2):
        h = bf.homology(alg, lam, 1, 2, p=p, tables=tb)
        out[p] = {g + Th: v[3] for g, v in h.items() if v[3]}
    chi = sum((-1) ** ((Y - Th) % 2) * d for Y, d in out[bf.P1].items())
    return lam, L, out, chi, time.time() - t0

if __name__ == "__main__":
    J0, J1, maxsize = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    sizes = {}
    for line in open("a4_size_table.log"):
        parts = line.split("|")
        lam = ast.literal_eval(parts[0].split(None, 1)[1].split(" res=")[0])
        sizes[lam] = int(parts[1].split("max")[1].split("@")[0])
    secs = [lam for lam in res if J0 <= 2 * sum(lam[:3]) + 3 * sum(lam[3:]) <= J1 and sizes.get(lam, 0) <= maxsize]
    secs.sort(key=lambda l: sizes.get(l, 0))
    print(f"# k=1 U(2) scan levels {J0}-{J1}, {len(secs)} sectors (max chain dim <= {maxsize:g})", flush=True)
    print("# level | sector | H^BG_k=1 {Y: dim} at p | p' agrees | chi(H^BG) | [Z_inf(1+G1)]_lam | Z_U2_lam | time", flush=True)
    with Pool(10) as pool:
        for lam, L, out, chi, dt in pool.imap_unordered(work, secs):
            zu, z1, r = res[lam]
            print(f"{L:3d} | {lam} | {out[bf.P1]} | {out[bf.P1] == out[bf.P2]} | {chi} | {z1} | {zu} | {dt:.1f}s", flush=True)
