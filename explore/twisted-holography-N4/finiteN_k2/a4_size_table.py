"""Size table: max_g dim V_{T,g} (cobar, small model) for k = 1, 2 in sectors with nonzero k=1 index residual.
usage: python3 a4_size_table.py Jmin Jmax"""
import sys, ast
from multiprocessing import Pool
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_sizes as bs
res = {}
for line in open("a1_residual.log"):
    if line.startswith("#"): continue
    p = [x.strip() for x in line.split("|")]
    res[ast.literal_eval(p[1])] = int(p[4])
def work(lam):
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
    alg, _ = bw.n4_letters(L)
    gmin = -2 * sum(lam[:3]) - 2; gmax = sum(lam[3:])
    ghs = list(range(gmin, gmax + 1))
    out = {}
    for k in (1, 2):
        d = bs.chain_dims(alg, lam, k, ghs, 2)
        out[('cobar', k)] = {g: int(d[((2,)*k, (2,)*k, g)]) for g in ghs}
        if lam[3] == 0 and lam[4] == 0:
            tcl, wcl = bs.small_model_classes(lam)
            d = bs.chain_dims_classes(lam, k, ghs, 2, tcl, wcl)
            out[('small', k)] = {g: int(d[((2,)*k, (2,)*k, g)]) for g in ghs}
    return lam, L, out
if __name__ == "__main__":
    J0, J1 = int(sys.argv[1]), int(sys.argv[2])
    secs = [lam for lam in res if J0 <= 2*sum(lam[:3]) + 3*sum(lam[3:]) <= J1]
    with Pool(12) as pool:
        for lam, L, out in pool.imap(work, secs):
            s = f"{L:3d} {lam} res={res[lam]:6d}"
            for key in sorted(out):
                v = out[key]; gm = max(v, key=lambda g: v[g])
                s += f" | {key[0]} k{key[1]}: max {v[gm]} @g{gm}"
            print(s, flush=True)
