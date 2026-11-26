"""Slot sizes (character formula of t7_sizes.py) of C_quot^{w=1} for n = 3 at the slots skipped in the exact
rank computations, with their SL_3-invariant parts (Weyl formula, validated in check_weyl.log).  Writes sizes_n3.log."""
from t7_sizes import dimC_N, dimV_N, weyl_mult
out = open('sizes_n3.log', 'w')
def log(s):
    print(s, flush=True); out.write(s + '\n')
for t in [(0, 0, -1), (0, -1, -1), (-1, -1, -1)]:
    G = -1 - sum(t)
    log('n=3 w=1 G=%d t=%s:' % (G, t))
    for N in range(0, 20):
        dC = sum(dimC_N(3, 1, t, N).values())
        if not dC:
            continue
        dV = sum(dimV_N(3, 1, t, N).values())
        iC = weyl_mult(3, 1, t, N, 'C') if len(set(t)) == 1 else None
        log('  N=%2d dim C=%8d dim V=%8d %s' % (N, dC, dV, ('SL3-inv C=%d' % iC) if iC is not None else ''))
