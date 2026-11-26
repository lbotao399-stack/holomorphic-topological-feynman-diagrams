"""dim V_{T,g} (k = 1, 2, cobar and HA^! model) and gauge-theory cochain dims for the O_0 sector and
neighbouring BMN sectors, at the ghost numbers g = Y - Theta for Y in Ywin."""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, bg_sizes as bs
secs = [(4,4,4,0,0), (5,4,3,0,0), (6,3,3,0,0), (4,4,3,0,0), (4,4,2,0,0), (3,3,3,0,0), (4,4,4,1,0)]
for lam in secs:
    L = 2 * sum(lam[:3]) + 3 * sum(lam[3:])
    Th = sum(lam[:3])
    alg, _ = bw.n4_letters(L)
    Ys = [6, 7, 8]
    ghs = [Y - Th for Y in Ys]
    line = f"{L:3d} {lam}"
    for k in (1, 2):
        d = bs.chain_dims(alg, lam, k, ghs, 2)
        line += f" | cobar k{k}: " + " ".join(f"Y{Y}:{d[((2,)*k,(2,)*k,Y-Th)]:.3g}" for Y in Ys)
        if lam[3] == 0:
            tcl, wcl = bs.small_model_classes(lam)
            d = bs.chain_dims_classes(lam, k, ghs, 2, tcl, wcl)
            line += f" | HA! k{k}: " + " ".join(f"Y{Y}:{d[((2,)*k,(2,)*k,Y-Th)]:.3g}" for Y in Ys)
    print(line, flush=True)
