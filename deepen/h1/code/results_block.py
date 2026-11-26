"""Text block for REPORT.txt Sec. 4 from logs/summary.log: per-log counts, and the nonzero entries for sl2, gl2,
abelian2 at (2,1), (2,2), and the eta-invariant complexes at w = 0."""
import os
LOGD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")
lines = open(os.path.join(LOGD, "summary.log")).read().splitlines()
out = []
for l in lines:
    if "compared entries" in l:
        out.append("  " + l)
out.append("  selected nonzero entries (full list: logs/summary.log):")
for l in lines:
    if not l.startswith("  check_"):
        continue
    keep = (" sl2 " in l or " gl2 " in l or ("abelian2" in l and ("(2, 2)" in l or "(2, 1)" in l) and ("w=(0, 0, 0)" in l or "w=(0, 0, -1)" in l))
            or (" Veta " in l and "w=(0, 0, 0)" in l and "abelian1" not in l))
    if keep:
        out.append("  " + l.strip())
print("\n".join(out))
