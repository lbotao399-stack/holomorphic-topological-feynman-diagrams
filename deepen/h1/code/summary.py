"""Collect the check logs: number of compared entries, mismatches, nonzero entries with the representative test."""
import re, glob, os

LOGD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "logs")
pat = re.compile(r"^\s+(\S+)\s+w=\((.*?)\)\s+\(a,b\)=\((.*?)\) g=\s*(-?\d+)\s+dim M=\s*(\d+).*H\(V\)=(\d+) \[modp (\d+) (\d+)\]\s+pred H\^(\S+)\(Q\)=(\d+)\s+(OK|MISMATCH)")
pate = re.compile(r"^\s+(\S+)\s+V\^eta(\(zb1,zb2 only\))?\s+w=\((.*?)\)\s+\(a,b\)=\((.*?)\) g=\s*(-?\d+)\s+dim V\^eta=\s*(\d+)\s+H\(V\^eta\)=(\d+)\s+pred=(\d+)(.*?)\s+(OK|MISMATCH)$")
patp = re.compile(r"psi_S, S=(\d+): dim Z\^(\d+)\(Q\)=(\d+), psi\(Z\) closed in V: (\w+), rank increase over B\^(-?\d+)\(V\) = (\d+) \(= H\): (\w+)")

out = []
tot = mis = 0
nz = []
for f in sorted(glob.glob(os.path.join(LOGD, "check_*.log"))):
    name = os.path.basename(f)
    n_ok = n_mis = 0
    last = None
    for line in open(f):
        m = pat.match(line)
        if m:
            alg, w, ab, g, dimM, H, h1, h2, n, pred, st = m.groups()
            n_ok += st == "OK"; n_mis += st == "MISMATCH"
            last = (alg, w, ab, g, H, n, pred)
            if int(H) or int(pred):
                nz.append([name, alg, "V", w, ab, g, dimM, H, n, pred, st, ""])
            continue
        m = patp.search(line)
        if m and nz and nz[-1][0] == name:
            S, n, dz, closed, gb, inc, ok = m.groups()
            nz[-1][-1] = "psi_S S=%s: closed %s, rank incr %s, %s" % (S, closed, inc, ok)
            continue
        m = pate.match(line)
        if m:
            alg, only, w, ab, g, dimK, H, pred, extra, st = m.groups()
            n_ok += st == "OK"; n_mis += st == "MISMATCH"
            if int(H) or int(pred):
                nz.append([name, alg, "Veta" + ("12" if only else ""), w, ab, g, dimK, H, "-", pred, st, extra.strip()])
    tot += n_ok + n_mis; mis += n_mis
    out.append("%-24s compared entries %4d   mismatches %d" % (name, n_ok + n_mis, n_mis))
out.append("total compared entries %d, mismatches %d" % (tot, mis))
out.append("")
out.append("nonzero entries (log, algebra, complex, w, (a,b), g, dim slice, computed H, CE degree, predicted, status, representative test):")
for r in nz:
    out.append("  %-20s %-8s %-6s w=(%s) (a,b)=(%s) g=%s dim=%s H=%s n=%s pred=%s %s %s" % tuple(r))
txt = "\n".join(out)
print(txt)
open(os.path.join(LOGD, "summary.log"), "w").write(txt + "\n")
