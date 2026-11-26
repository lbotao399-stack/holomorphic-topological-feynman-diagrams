"""Tally of all comparison lines in the verifier's logs."""
import glob, re
tot = bad = rep = repok = 0
rows = []
for f in sorted(glob.glob(__file__.rsplit('/code/', 1)[0] + '/logs/*.log')):
    L = [l for l in open(f) if (' OK' in l or 'MISMATCH' in l) and ' g=' in l]
    if not L:
        continue
    n = len(L); m = sum('MISMATCH' in l for l in L)
    r = sum('psi_S rep' in l for l in L); ro = sum('REP-OK' in l for l in L)
    nz = sum(bool(re.search(r'pred=[1-9]', l)) for l in L)
    exact = sum(',Q)=' in l for l in L)
    rows.append((f.rsplit('/', 1)[1], n, m, nz, r, ro, exact))
    tot += n; bad += m; rep += r; repok += ro
print(f"{'log':24s} {'entries':>8s} {'mismatch':>8s} {'pred>0':>7s} {'rep':>5s} {'repOK':>6s} {'exactQ':>7s}")
for row in rows:
    print(f"{row[0]:24s} {row[1]:8d} {row[2]:8d} {row[3]:7d} {row[4]:5d} {row[5]:6d} {row[6]:7d}")
print(f"TOTAL entries {tot}, mismatches {bad}, representative tests {rep}, passed {repok}")
