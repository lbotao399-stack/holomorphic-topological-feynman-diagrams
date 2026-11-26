"""Compact sector tables from the k = 1 logs (one line per sector representative)."""
import ast
import sys


def load(paths):
    rows = {}
    for path in paths:
        try:
            fh = open(path)
        except FileNotFoundError:
            continue
        for line in fh:
            if line.startswith("#") or "|" not in line:
                continue
            p = [x.strip() for x in line.split("|")]
            j = int(p[0])
            lam = ast.literal_eval(p[1])
            rows[(j, lam)] = (ast.literal_eval(p[2]), ast.literal_eval(p[3]), p[4])
    return rows


def fmt(lam):
    return "({},{},{};{},{})".format(*lam)


def fmtH(h):
    return " ".join(f"{h[Y]}@Y{Y}" for Y in sorted(h)) or "0"


out = []
for N, paths in ((1, ["bg_vs_U1_k1.log"]), (2, ["bg_vs_U2_k1.log", "bg_vs_U2_k1_L11_12.log"]),
                 (3, ["bg_vs_U3_k1.log"]), (4, ["bg_vs_U4_k1.log"])):
    rows = load(paths)
    out.append(f"U({N}), k = 1: sector | dim H^Y(Q_0)^U({N}) (= dim H_g^BG, g = Y - Theta) | agree")
    for (j, lam) in sorted(rows, key=lambda x: (x[0], [-v for v in x[1]])):
        hu, hb, ok = rows[(j, lam)]
        tag = "=" if ok == "Y" else "MISMATCH BG: " + fmtH(hb)
        out.append(f"  {j:2d} {fmt(lam):16s} {fmtH(hu):40s} {tag}")
    out.append(f"  sectors: {len(rows)}, mismatches: {sum(1 for v in rows.values() if v[2] != 'Y')}")
    out.append("")
open("tables_k1.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
