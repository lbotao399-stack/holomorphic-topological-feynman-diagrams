"""Per-sector tree-level cohomology: stable (U(inf), no trace relations), U(2), SU(2), U(3).
usage: python3 b3_sector_table.py Jmax"""
import sys
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import bg_words as bw, uN_q0 as uq
J = int(sys.argv[1])
alg, keys = bw.n4_letters(J)
exec(open("b2_bg_vs_UN.py").read().split("print(f\"# N")[0].split("bg = bc.BG1(alg)")[1])   # reps()
def fH(h):
    return " ".join(f"{h[Y]}@Y{Y}" for Y in sorted(h) if h[Y]) or "0"
print("level | sector | U(inf) stable | U(2) | SU(2) | U(3)   [dim H^Y(Q_0), entries n@Y]")
for j in range(2, J + 1):
    for lam in reps(j):
        st = uq.cohomology_sector_stable(alg, lam)
        u2 = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, 2).items()}
        su2 = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, 2, special=True).items()}
        u3 = {Y: v[2] for Y, v in uq.cohomology_sector(alg, lam, 3).items()}
        print(f"{j:3d} | ({lam[0]},{lam[1]},{lam[2]};{lam[3]},{lam[4]}) | {fH(st)} | {fH(u2)} | {fH(su2)} | {fH(u3)}", flush=True)
