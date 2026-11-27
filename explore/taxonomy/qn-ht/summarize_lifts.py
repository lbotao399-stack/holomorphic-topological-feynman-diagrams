"""Parse lift_*.log: per theory, completed sectors (all sectors with sum(w) <= W_done, plus partial),
sectors with H(Q0) != 0, total dim H, cocycle basis vectors with Q2^hat z != 0 (target H != 0), total rank d2.
Log: lift_summary.log"""
import glob, re
from fractions import Fraction as Fr

out = open("lift_summary.log", "w")
for f in sorted(glob.glob("lift_*.log")):
    if f == "lift_summary.log":
        continue
    lines = open(f).read().splitlines()
    head = lines[0] if lines else ""
    done = any(l.startswith("time") for l in lines)
    secs = [l for l in lines if l.startswith("w=")]
    nH = totH = nq2 = rk = 0
    sums = []
    for l in secs:
        w = re.findall(r"Fraction\((\d+), (\d+)\)", l.split(" charges")[0])
        if w:
            s = sum(Fr(int(a), int(b)) for a, b in w)
        else:
            s = sum(int(x) for x in re.findall(r"\d+", l.split(" charges")[0]))
        sums.append(s)
        hpart = l.split("dim H(Q0) = ")[1].split(" ;")[0]
        hd = eval(hpart)
        nH += 1
        totH += sum(hd.values())
        m = re.search(r"Q2\^hat z != 0 on (\d+)", l)
        if m:
            nq2 += int(m.group(1))
        m = re.search(r"rank d2 = (\{[^}]*\})", l)
        if m:
            rk += sum(eval(m.group(1)).values())
    smax = max(sums) if sums else 0
    out.write("%s: %s\n   finished: %s; sectors listed: %d (sum w <= %s%s); dim H(Q0) total %d; "
              "cocycle basis vectors with Q2^hat z != 0: %d; total rank d2: %d\n"
              % (f, head.split(";")[0], done, nH, smax, "" if done else ", last weight partial", totH, nq2, rk))
out.close()
print(open("lift_summary.log").read())
