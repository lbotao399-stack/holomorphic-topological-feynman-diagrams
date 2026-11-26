"""Collect the k = 1 Budzik-Gaiotto vs U(N) sector logs: per-level totals with S_3 x Z_2 orbit
multiplicities, comparison with the Chang-Lin U(N) tables, and the per-sector Euler characteristic
against the refined Murthy truncation Z_inf (1 + G_N^(1)) (Budzik-Gaiotto (10) at k = 1).
usage: python3 summarize.py N log1 [log2 ...]"""
import ast
import sys
from itertools import permutations
sys.path.insert(0, (__import__('os').environ.get('HTP_ROOT') or next(str(p) for p in __import__('pathlib').Path(__file__).resolve().parents if (p / '.htp-root').exists())) + "/code")
from holography import giant as g, census_io as cio

N = int(sys.argv[1])
rows = []
for path in sys.argv[2:]:
    for line in open(path):
        if line.startswith("#") or "|" not in line:
            continue
        parts = [p.strip() for p in line.split("|")]
        j = int(parts[0])
        lam = ast.literal_eval(parts[1])
        hu = ast.literal_eval(parts[2])
        hb = ast.literal_eval(parts[3])
        rows.append((j, lam, hu, hb, parts[4]))


def mult(lam):
    return len(set(permutations(lam[:3]))) * len(set(permutations(lam[3:])))


CL = cio.chang_lin_tables()
levels = sorted(set(r[0] for r in rows))
print(f"U({N}), k = 1: per-level totals (sum over sectors with orbit multiplicity)")
print("level | sectors(reps) | all agree | dim H^BG | dim H^U(N) | Chang-Lin d_u | index BG | Chang-Lin i_u")
for j in levels:
    rs = [r for r in rows if r[0] == j]
    dimB = dimU = indB = 0
    for (_, lam, hu, hb, ok) in rs:
        m = mult(lam)
        Th = sum(lam[:3])
        dimB += m * sum(hb.values())
        dimU += m * sum(hu.values())
        indB += m * sum(((-1) ** ((Y + Th) % 2)) * h for Y, h in hb.items())
    cl = CL.get(N)
    dcl = cl["d_u"][j] if cl and j < len(cl["d_u"]) else None
    icl = cl["i_u"][j] if cl and j < len(cl["i_u"]) else None
    allok = all(r[4] == "Y" for r in rs)
    print(f"{j:3d} | {len(rs)} | {allok} | {dimB} | {dimU} | {dcl} | {indB} | {icl}")

print("\nper-sector Euler characteristic of the k=1 complex vs [Z_inf (1 + G_N^(1))]_lam:")
bad = 0
for (j, lam, hu, hb, ok) in rows:
    Th = sum(lam[:3])
    chi = sum(((-1) ** ((Y + Th) % 2)) * h for Y, h in hb.items())
    box = lam
    Zi = g.Zinf_box(box)
    G1 = g.G1_box(N, box)
    S = dict(G1)
    S[(0, 0, 0, 0, 0)] = S.get((0, 0, 0, 0, 0), 0) + 1
    val = 0
    for a, ca in Zi.items():
        b = tuple(y - x for x, y in zip(a, lam))
        if all(x >= 0 for x in b):
            val += ca * S.get(b, 0)
    flag = "" if val == chi else "  <-- mismatch"
    if val != chi:
        bad += 1
    print(f"{j:3d} {lam}: chi(H^BG) = {chi}, [Z_inf(1+G1)] = {val}{flag}")
print("mismatches:", bad)
