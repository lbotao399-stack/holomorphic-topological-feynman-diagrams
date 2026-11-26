"""dim H^n(C_g^0)_{(a,b)} for sl2 (ad(h)-charge 0) and abelian1, n = 1..5 (exact ranks, cohom.py)."""
import sys
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import jets as J, cohom as C
for th, ch in ((J.sl2(), 0), (J.abelian(1), None)):
    for ab in [(0, 0), (1, 0), (2, 0), (1, 1), (3, 0), (2, 1), (2, 2), (3, 1), (4, 0)]:
        pc = C.Piece(th, None, ab, holo=True, charge=ch)
        res = C.cohomology(pc, range(1, 6))
        print(th.name, ab, "H^n:", {n: r["H"][2] for n, r in res.items()},
              " modp agree:", all(r["H"][0] == r["H"][1] == r["H"][2] for r in res.values()),
              " dim C^n:", {n: r["dimM"] - r["rkD"][2] for n, r in res.items()}, flush=True)
