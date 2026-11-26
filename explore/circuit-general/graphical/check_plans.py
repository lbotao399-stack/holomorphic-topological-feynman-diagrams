"""Plan independence: one period computed through several distinct sequential plans
(different decompletion vertices, hub pairs, path orientations) must give the same value.
usage: python3 check_plans.py 'EDGES' 'ALPHAS' Q [max_plans]"""
import ast, os, sys, time
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import seq, mzv
E = ast.literal_eval(sys.argv[1]); A = ast.literal_eval(sys.argv[2]); q = int(sys.argv[3])
mx = int(sys.argv[4]) if len(sys.argv) > 4 else 6
cands = seq.plans_for(E, A, q)
seen = set(); vals = []
pref = F(1, seq.gamma_int(q))
for a in A:
    pref *= seq.gamma_int(q - a)
for pl in cands:
    key = (pl['first'], tuple((o[0], str(o[1])) for o in pl['ops']))
    if key in seen:
        continue
    seen.add(key)
    t = time.time()
    v = mzv.scal(seq.run(pl), pref)
    vals.append(v)
    print('  %-14s hubs %-8s path %-16s appends %d: %s (%.1fs)' % (pl['route'], pl['hubs'], pl['path'], pl['appends'], mzv.tostr(v), time.time() - t), flush=True)
    if len(vals) >= mx:
        break
print('%d plans, all equal: %s' % (len(vals), all(v == vals[0] for v in vals)))
